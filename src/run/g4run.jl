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

# Resolve the `manifest` keyword to a file path (or `nothing`).
#   - `GeometryManifest` -> serialised to `dest`, returns `dest`
#   - `AbstractString`   -> treated as an existing manifest file, returned as-is
#   - `nothing`          -> returns `nothing`
function _resolve_manifest(manifest, dest::String)
    manifest === nothing && return nothing
    if manifest isa GeometryManifest
        mkpath(dirname(dest))
        return write_manifest(dest, manifest)
    elseif manifest isa AbstractString
        isfile(manifest) || error("run: manifest file not found: $manifest")
        return String(manifest)
    else
        error("run: `manifest` must be a GeometryManifest, a path String, or " *
              "nothing; got $(typeof(manifest))")
    end
end

# Append `--key value` pairs for each keyword argument.
function _append_kwargs!(args::Vector{String}, kwargs)
    for (k, v) in kwargs
        push!(args, "--$(k)", string(v))
    end
    return args
end

"""
    run_simulation(; outdir, projectdir=_default_projectdir(),
                     manifest=nothing, kwargs...) -> outdir

Run a batch g4scint simulation, writing results into `outdir`.

`projectdir` is the G4ScintKit repository root (it defaults to the parent of
the G4ScintKit.jl package, so callers normally need not set it).

`manifest` selects the detector geometry:
  - a `GeometryManifest` is serialised to `joinpath(outdir, "geometry.manifest")`
    — persisted alongside the output as provenance — and passed to run.sh;
  - an `AbstractString` is treated as an existing manifest file path;
  - `nothing` is no longer accepted: pass a `GeometryManifest` or a path. If `nothing`, run.sh will fail with a clear error.

Any remaining keyword arguments are forwarded to `run.sh` as `--key value`
pairs (e.g. `nevents=10`, `injparticle="mu-"`, `setup="B2"`).
"""
function run_simulation(; outdir::String,
                          projectdir::String = _default_projectdir(),
                          manifest = nothing,
                          kwargs...)

    mkpath(outdir)
    manifest_file = _resolve_manifest(manifest,
                                      joinpath(outdir, "geometry.manifest"))

    args = String["--outdir", outdir]
    manifest_file === nothing || append!(args, ["--manifest", manifest_file])
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

run_simulation(manifest::GeometryManifest; outdir::String, kwargs...) =
    run_simulation(; outdir, manifest, kwargs...)

run_simulation(spec::DetectorSpec; outdir::String, kwargs...) =
    run_simulation(build_manifest(spec); outdir, kwargs...)

"""
    run_visu(; projectdir=_default_projectdir(), manifest=nothing, kwargs...) -> outdir

Launch g4scint in visualization mode via `bash_scripts/run_visu.sh`.

Output goes to `<projectdir>/output/visu`, which `run_visu.sh` wipes on each
call; the path is returned. `manifest` behaves as in [`run_simulation`](@ref),
except a `GeometryManifest` is serialised to a temporary file (the visu output
directory is recreated each run). Remaining keyword arguments are forwarded to
`run.sh` as `--key value` pairs.
"""
function run_visu(; projectdir::String = _default_projectdir(),
                    manifest = nothing,
                    kwargs...)

    manifest_file = _resolve_manifest(manifest, tempname() * ".manifest")

    args = String[]
    manifest_file === nothing || append!(args, ["--manifest", manifest_file])
    _append_kwargs!(args, kwargs)

    script = joinpath(projectdir, "bash_scripts", "run_visu.sh")
    run(`bash $script $args`)

    return joinpath(projectdir, "output", "visu")
end

run_visu(manifest::GeometryManifest; kwargs...) = run_visu(; manifest, kwargs...)

run_visu(spec::DetectorSpec; kwargs...) = run_visu(build_manifest(spec); kwargs...)
