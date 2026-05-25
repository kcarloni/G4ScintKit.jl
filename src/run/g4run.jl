# g4run.jl
#
# Julia wrappers around the bash entry points that drive the g4scint
# executable. They source `bash_scripts/setup_paths.sh` (Geant4 + GODDESS
# environment) and then `g4scintkit/run.sh`, forwarding keyword arguments
# as `--flag value` pairs.
#
# A geometry can be supplied directly as a `GeometryManifest` (see manifest.jl)
# via the `manifest` keyword; it is serialised to disk and passed to run.sh as
# `--manifest <file>` (geometry input is now mandatory).

"""Default G4ScintKit project root (the parent of the G4ScintKit.jl package)."""
# this file lives at src/run/ — three levels below the G4ScintKit project root.
_default_projectdir() = dirname(dirname(dirname(@__DIR__)))

# Resolve a `manifest` argument to a file path.
#   - `GeometryManifest` -> serialised to `dest`, returns `dest`
#   - `AbstractString`   -> treated as an existing manifest file, returned as-is
function _resolve_manifest(manifest::GeometryManifest, dest::String)
    mkpath(dirname(dest))
    return write_manifest(dest, manifest)
end

function _resolve_manifest(manifest::AbstractString, ::String)
    isfile(manifest) || error("run: manifest file not found: $manifest")
    return String(manifest)
end

# Format a keyword value for run.sh's CLI. Unitful values become
# "<num> <unit>" strings in Geant4-internal units (mm/MeV/ns/rad), so the
# rest of the pipeline doesn't need to know about Unitful. Plain values
# pass through `string`. Used with `--key=value` so values may contain
# spaces.
_format_g4cli(v) = string(v)
_format_g4cli(v::Bool) = v ? "true" : "false"
function _format_g4cli(q::Quantity)
    d = dimension(q)
    if     d === dimension(u"m");   return string(ustrip(u"mm",  q), " mm")
    elseif d === dimension(u"J");   return string(ustrip(u"MeV", q), " MeV")
    elseif d === dimension(u"s");   return string(ustrip(u"ns",  q), " ns")
    elseif d === dimension(u"V");   return string(ustrip(u"V",   q), " V")
    elseif d === dimension(u"rad"); return string(ustrip(u"rad", q), " rad")
    else
        error("_format_g4cli: unsupported unit dimension $(d) for value $(q)")
    end
end
function _format_g4cli(t::Tuple)
    isempty(t) && return ""
    if all(x -> x isa Quantity && dimension(x) === dimension(u"m"), t)
        return string(join((string(ustrip(u"mm", x)) for x in t), " "), " mm")
    elseif all(x -> x isa Real, t)
        return join((string(x) for x in t), " ")
    else
        error("_format_g4cli: tuple must be all lengths or all reals; got $(typeof(t))")
    end
end

# Append `--key=value` pairs for each keyword argument. We use `=` rather
# than a separate token so values may contain spaces (e.g. "3 GeV") without
# extra quoting through bash positional parameters.
function _append_kwargs!(args::Vector{String}, kwargs)
    for (k, v) in kwargs
        push!(args, string("--", k, "=", _format_g4cli(v)))
    end
    return args
end

"""
    run_simulation(manifest; outdir, projectdir=_default_projectdir(), kwargs...) -> outdir
    run_simulation(spec::DetectorSpec; outdir, kwargs...) -> outdir

Run a batch g4scint simulation, writing results into `outdir`.

`manifest` is either a `GeometryManifest` (serialised to
`joinpath(outdir, "geometry.manifest")` as provenance) or an `AbstractString`
path to an existing manifest file. The `DetectorSpec` form first calls
[`build_manifest`](@ref).

`projectdir` defaults to the G4ScintKit repository root (the parent of the
G4ScintKit.jl package), so callers normally need not set it. Remaining
keyword arguments are forwarded to `run.sh` as `--key=value` pairs
(e.g. `nevents=10`, `injparticle="mu-"`, `injenergy=3u"GeV"`).
"""
function run_simulation(manifest::Union{GeometryManifest,AbstractString};
                          outdir::String,
                          projectdir::String = _default_projectdir(),
                          kwargs...)

    mkpath(outdir)
    manifest_file = _resolve_manifest(manifest,
                                      joinpath(outdir, "geometry.manifest"))

    args = String["--outdir", outdir, "--manifest", manifest_file]
    _append_kwargs!(args, kwargs)

    setup = joinpath(projectdir, "bash_scripts", "setup_paths.sh")
    # Pass run.sh arguments as bash positional parameters ($2…) rather than
    # interpolating them into the -c script, so values may contain spaces.
    # If the user hasn't created a local setup_paths.sh, fall back to the
    # committed .example (see bash_scripts/setup_paths.sh.example).
    script = raw"""
        PATHS="$1"
        [[ -f "$PATHS" ]] || PATHS="${PATHS}.example"
        source "$PATHS" && source "$SIMDIR/run.sh" "${@:2}"
    """
    run(`bash -c $script g4scint $setup $args`)

    return outdir
end

run_simulation(spec::DetectorSpec; outdir::String, kwargs...) =
    run_simulation(build_manifest(spec); outdir, kwargs...)

"""
    run_simulation(manifest, particles::AbstractVector{ParticleListEntry};
                   outdir, kwargs...) -> outdir
    run_simulation(spec::DetectorSpec, particles::AbstractVector{ParticleListEntry};
                   outdir, kwargs...) -> outdir

Particle-list overload: serialise `particles` to `<outdir>/particle_list.csv`
via [`write_particle_list`](@ref) and forward it as `--particlelist`. An
explicit `particlelist=` kwarg takes precedence (the CSV is still written,
but the explicit path wins).
"""
function run_simulation(manifest::Union{GeometryManifest,AbstractString},
                          particles::AbstractVector{ParticleListEntry};
                          outdir::String, kwargs...)
    mkpath(outdir)
    csv = joinpath(outdir, "particle_list.csv")
    write_particle_list(csv, particles)
    if haskey(kwargs, :particlelist)
        return run_simulation(manifest; outdir, kwargs...)
    end
    return run_simulation(manifest; outdir, particlelist = csv, kwargs...)
end

run_simulation(spec::DetectorSpec,
               particles::AbstractVector{ParticleListEntry};
               outdir::String, kwargs...) =
    run_simulation(build_manifest(spec), particles; outdir, kwargs...)

"""
    run_visu(manifest; projectdir=_default_projectdir(), kwargs...) -> outdir
    run_visu(spec::DetectorSpec; kwargs...) -> outdir

Launch g4scint in visualization mode via `bash_scripts/run_visu.sh`.

Output goes to `<projectdir>/output/visu`, which `run_visu.sh` wipes on each
call; the path is returned. `manifest` behaves as in [`run_simulation`](@ref),
except a `GeometryManifest` is serialised to a temporary file (the visu output
directory is recreated each run). Remaining keyword arguments are forwarded to
`run.sh` as `--key=value` pairs.
"""
function run_visu(manifest::Union{GeometryManifest,AbstractString};
                    projectdir::String = _default_projectdir(),
                    kwargs...)

    manifest_file = _resolve_manifest(manifest, tempname() * ".manifest")

    args = String["--manifest", manifest_file]
    _append_kwargs!(args, kwargs)

    script = joinpath(projectdir, "bash_scripts", "run_visu.sh")
    run(`bash $script $args`)

    return joinpath(projectdir, "output", "visu")
end

run_visu(spec::DetectorSpec; kwargs...) = run_visu(build_manifest(spec); kwargs...)
