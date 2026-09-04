# load.jl
#
# Public `load(outdir)` for reading the C++ side's HDF5 output. The on-disk
# hierarchy is:
#
#   <outdir>/Data/*.h5
#     └── g4run_<RunId>/
#          ├── input/
#          ├── particle_hits/
#          ├── optical_photons/
#          ├── sipm_hits/
#          ├── sipm_digi_summary/
#          └── sipm_voltage_trace/
#
# A single call to `load` walks every .h5 file, every g4run_* group, and
# concatenates rows. Unit suffixes on column names are stripped and the
# corresponding `Unitful` units attached (see columns.jl).

"""Group names emitted by the C++ HDF5 writer (see g4scintkit/src/DataProcessing/HDF5Writer.cc).
Anything appearing here can be passed to `load(...; groups=...)`."""
const ALL_GROUPS = (:input, :particle_hits, :optical_photons,
                    :sipm_hits, :sipm_digi_summary, :sipm_voltage_trace)

"""
    SimulationOutput

Container returned by [`load`](@ref). Wraps a NamedTuple of per-group
StructArrays so the REPL display stays compact for real-size runs while
keeping `data.input`, `iterate`, etc. ergonomic.
"""
struct SimulationOutput
    groups::NamedTuple
end

Base.getproperty(o::SimulationOutput, s::Symbol) =
    s === :groups ? getfield(o, :groups) : getfield(getfield(o, :groups), s)
Base.propertynames(o::SimulationOutput) = propertynames(getfield(o, :groups))
Base.keys(o::SimulationOutput)          = propertynames(getfield(o, :groups))
Base.values(o::SimulationOutput)        = values(getfield(o, :groups))
Base.length(o::SimulationOutput)        = length(getfield(o, :groups))
Base.iterate(o::SimulationOutput, s...) = iterate(pairs(getfield(o, :groups)), s...)
Base.haskey(o::SimulationOutput, k::Symbol) = haskey(getfield(o, :groups), k)
Base.getindex(o::SimulationOutput, k::Symbol) = getfield(getfield(o, :groups), k)

function Base.show(io::IO, ::MIME"text/plain", o::SimulationOutput)
    ks = propertynames(o)
    println(io, "SimulationOutput with $(length(ks)) group", length(ks) == 1 ? "" : "s", ":")
    namew = isempty(ks) ? 0 : maximum(length ∘ string, ks)
    for k in ks
        sa = getproperty(o, k)
        ncols = length(propertynames(sa))
        println(io, "  ", rpad(string(k), namew + 2),
                lpad(string(length(sa)), 7), " rows × ",
                lpad(string(ncols), 2), " cols")
    end
end

# Resolve a path argument to a list of .h5 files. Accepts either an outdir
# (looks in <outdir>/Data/*.h5) or a single .h5 file.
function _resolve_files(path::AbstractString)
    if isfile(path)
        endswith(path, ".h5") || error("load: $path is not an .h5 file")
        return [path]
    elseif isdir(path)
        datadir = isdir(joinpath(path, "Data")) ? joinpath(path, "Data") : path
        files = sort(filter(f -> endswith(f, ".h5"),
                           readdir(datadir; join=true)))
        isempty(files) && error("load: no .h5 files found under $datadir")
        return files
    else
        error("load: path does not exist: $path")
    end
end

# Extract the integer RunId from a "g4run_<N>" group name. Errors on unexpected
# group names so we don't silently mis-parse.
function _parse_run_id(name::AbstractString)
    m = match(r"^g4run_(\d+)$", name)
    m === nothing && error("load: unexpected run group name: $name")
    parse(Int, m.captures[1])
end

# List the "g4run_*" group names in an open HDF5 file, in numeric run-id order.
function _run_groups(f::HDF5.File)
    names = filter(n -> startswith(n, "g4run_"), keys(f))
    sort!(names; by=_parse_run_id)
    names
end

# Per-event waveform columns arrive as an (n_samples x n_events) matrix while
# every sibling column is a length-n_events vector. StructArray requires all
# components to share a shape, so split the matrix into one view per event.
# Only sipm_voltage_trace.voltages is shaped this way today, but the check is
# on the shape rather than the name, so a future waveform column needs no
# change here.
function _split_matrix_columns(nt::NamedTuple)
    any(v -> v isa AbstractMatrix, values(nt)) || return nt
    vecs = Iterators.filter(v -> v isa AbstractVector, values(nt))
    isempty(collect(Iterators.take(vecs, 1))) &&
        error("load: group has a matrix column but no vector column to give " *
              "the event count")
    n = length(first(vecs))
    return map(nt) do v
        v isa AbstractMatrix || return v
        size(v, 2) == n || error(
            "load: waveform column has $(size(v, 2)) columns but the group " *
            "has $n events; expected (n_samples x n_events)")
        return [view(v, :, i) for i in 1:n]
    end
end

# Build the cleaned StructArray for one group under one run, adding a
# `run_id` column. Returns `nothing` if the group is absent in this run.
function _read_one_run_group(run_g::HDF5.Group, gname::Symbol, run_id::Int)
    skey = String(gname)
    haskey(run_g, skey) || return nothing
    raw = _group_to_namedtuple(run_g[skey])
    cleaned = _split_matrix_columns(_clean_columns(raw))
    n = length(first(values(cleaned)))
    sa = StructArray(merge(cleaned, (run_id = fill(run_id, n),)))
    return sa
end

# Vertically concat a sequence of StructArrays sharing a schema. Empty
# fallback returns `nothing` to signal "no data anywhere".
function _vcat_sas(sas::Vector)
    isempty(sas) && return nothing
    out = copy(sas[1])
    for i in 2:length(sas)
        append!(out, sas[i])
    end
    return out
end


# ---------------------------------------------------------------------------
#  Voltage-trace units
# ---------------------------------------------------------------------------
#
# `sipm_voltage_trace` is the one group whose units depend on which g4sipm model
# produced it, so it gets handled here rather than inside the generic reader.
#
# Every g4sipm model but two builds its trace from G4SipmGenericVoltageTraceModel,
# whose amplitude/baseline/noise are real voltages (50 mV / 25 mV / 1 mV, each
# scaled by CLHEP::volt). HamamatsuS12573100C and ...100X instead define their own
# nested VoltageTraceModel returning bare numbers -- amplitude 1.0, v0 0.0, noise
# sigma 0.01 -- normalised to one fired cell, with no voltage scale at all. That is
# upstream g4sipm behaviour, present since its initial commit.
#
# The C++ writer divides the trace by CLHEP::volt exactly as it does every other
# quantity, which is right for the first group and, for the other two, inflates a
# per-photoelectron amplitude by 1/CLHEP::volt = 1e6. Nothing is lost: the stored
# value is exactly photoelectrons * 1e6, so the scaling below is exact rather than
# approximate.
#
# The model is not recorded in the HDF5, but RunSimulation.cc dumps the manifest it
# actually built from to <Data>/geometry.manifest on every run, so the SIPM lines
# there are the authority. When that cannot be resolved we drop the unit rather
# than guess -- a number labelled volts that is not a voltage is worse than a bare
# number.

const _PE_PER_STORED_UNIT = 1e-6

# Path to the manifest the C++ dumped beside the data, or `nothing`.
function _manifest_beside(path::AbstractString)
    dir = if isfile(path)
        dirname(path)
    elseif isdir(path)
        isdir(joinpath(path, "Data")) ? joinpath(path, "Data") : path
    else
        return nothing
    end
    mf = joinpath(dir, "geometry.manifest")
    return isfile(mf) ? mf : nothing
end

# Model aliases named on the manifest's SIPM lines. A SIPM line with no `model=`
# is a GODDeSS photon detector, which produces no g4sipm trace at all.
function _manifest_sipm_models(manifest::AbstractString)
    models = String[]
    for line in eachline(manifest)
        startswith(line, "SIPM") || continue
        m = match(r"(?:^|\s)model=(\S+)", line)
        m === nothing || push!(models, String(m.captures[1]))
    end
    return models
end

# :volt, :photoelectron, :mixed (SiPMs disagree) or :unknown (no manifest, no
# g4sipm SiPM in it, or an alias this package does not know).
function _trace_units(path::AbstractString)
    mf = _manifest_beside(path)
    mf === nothing && return :unknown
    models = _manifest_sipm_models(mf)
    isempty(models) && return :unknown
    infos = [get(SIPM_MODELS, a, nothing) for a in models]
    any(i -> i === nothing, infos) && return :unknown
    units = unique(i.trace_units for i in infos)
    return length(units) == 1 ? only(units) : :mixed
end

# Rescale/relabel the `voltages` column of a sipm_voltage_trace StructArray.
# Values arrive carrying u"V" from _clean_columns; that is correct only for
# :volt models.
function _rescale_voltage_trace(sa, path::AbstractString)
    units = _trace_units(path)
    units === :volt && return sa
    cols = NamedTuple(k => getproperty(sa, k) for k in propertynames(sa))
    haskey(cols, :voltages) || return sa
    bare = [ustrip.(u"V", t) for t in cols.voltages]
    if units === :photoelectron
        @info "load: SiPM model reports its voltage trace per photoelectron, " *
              "not as a voltage; `voltages` is in photoelectrons (dimensionless)."
        bare = [t .* _PE_PER_STORED_UNIT for t in bare]
    else
        @warn "load: could not determine the SiPM model for this output " *
              "($(units)), so the voltage trace unit is unknown; returning " *
              "`voltages` as stored, without a unit."
    end
    return StructArray(merge(cols, (voltages = bare,)))
end

"""
    load(path; groups = ALL_GROUPS) -> SimulationOutput

Read the HDF5 output produced by `run_simulation` / `g4scint`. The returned
[`SimulationOutput`](@ref) supports dot-access (`data.input`),
indexing (`data[:input]`), `haskey`, `keys`/`values`, and iteration over
`(group_name, struct_array)` pairs.

`path` may be either an output directory (in which case `<path>/Data/*.h5`
is read; bare `*.h5` is also accepted if `Data/` is absent) or a single
`.h5` file. Rows are concatenated across files and across `g4run_*` groups
within each file; a synthesized `run_id::Int` column identifies the source
run. The HDF5 column `g4event_id` is renamed to `event_id`.

Numeric columns whose HDF5 name ends in a unit suffix (`_mm`, `_MeV`,
`_ns`, `_V`, …) have the suffix stripped and the corresponding `Unitful`
unit attached. Suffix-less columns (ints, dimensionless doubles) pass
through unchanged.

`sipm_voltage_trace.voltages` is a special case. Most g4sipm models report a
real voltage there and it carries `u"V"`, but `hamamatsu-s12573-100c`/`-100x`
normalise their trace to one photoelectron instead, so for those it is rescaled
and returned dimensionless, in photoelectrons. The model is read from the
`geometry.manifest` the C++ writes beside the data; if that is missing the unit
is dropped rather than guessed, with a warning.

`groups` selects which subgroups to read; defaults to all known names
(see `ALL_GROUPS`). Groups that don't appear in any run are simply omitted
from the returned `SimulationOutput` — test with `haskey(data, :name)`.

```julia
data = load("output/B3_test")
data.input.init_kinetic_energy             # Vector{Quantity{Float64, …, u"MeV"}}
data.particle_hits.event_id                # Vector{Int}
data.sipm_voltage_trace.voltages           # Matrix or Vector{Vector}, with u"V"
```
"""
function load(path::AbstractString; groups = ALL_GROUPS)
    files = _resolve_files(path)

    # accumulator: gname => Vector{StructArray}
    bins = Dict{Symbol, Vector{Any}}(g => [] for g in groups)

    for file in files
        h5open(file, "r") do f
            for run_name in _run_groups(f)
                run_id = _parse_run_id(run_name)
                run_g  = f[run_name]
                for g in groups
                    sa = _read_one_run_group(run_g, g, run_id)
                    sa === nothing || push!(bins[g], sa)
                end
            end
        end
    end

    # Build the final NamedTuple: only groups that produced ≥1 row are
    # included. Use `haskey(data, :group)` to test presence.
    pairs = Pair{Symbol, Any}[]
    for g in groups
        v = _vcat_sas(bins[g])
        v === nothing && continue
        g === :sipm_voltage_trace && (v = _rescale_voltage_trace(v, path))
        push!(pairs, g => v)
    end
    SimulationOutput(NamedTuple{Tuple(first.(pairs))}(Tuple(last.(pairs))))
end

"""
    available_groups(path) -> Vector{Symbol}

Return the union of `g4run_*` subgroup names present anywhere under `path`
(an outdir or a single .h5 file). Useful for discovering what's in a file
whose schema you don't know up front.
"""
function available_groups(path::AbstractString)
    files = _resolve_files(path)
    found = Set{Symbol}()
    for file in files
        h5open(file, "r") do f
            for run_name in _run_groups(f)
                for k in keys(f[run_name])
                    push!(found, Symbol(k))
                end
            end
        end
    end
    sort!(collect(found))
end
