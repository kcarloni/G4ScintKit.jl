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

# Build the cleaned StructArray for one group under one run, adding a
# `run_id` column. Returns `nothing` if the group is absent in this run.
function _read_one_run_group(run_g::HDF5.Group, gname::Symbol, run_id::Int)
    skey = String(gname)
    haskey(run_g, skey) || return nothing
    raw = _group_to_namedtuple(run_g[skey])
    cleaned = _clean_columns(raw)
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
        v === nothing || push!(pairs, g => v)
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
