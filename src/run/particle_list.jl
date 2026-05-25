# particle_list.jl
#
# Julia mirror of the C++ `ParticleListSource` CSV input — the primary-particle
# list consumed via `--particlelist` (run.sh) / `--useParticleList` +
# `--particleSourceInput` (g4scint binary). See
# goddess-package/source/G4BasicObjects/G4GeneralParticleSource/src/ParticleListSource.cc
# for the reader.
#
# Schema (11 columns, comma-separated, optional `#` comment lines and one
# optional header row): run_id, event_id, pdg, t_ns, x_mm, y_mm, z_mm,
# ekin_MeV, dx, dy, dz. Rows sharing (run_id, event_id) become multiple
# primaries in the same event. Direction is renormalised by the C++ side, so
# `dir` here is a plain Float64 triple (not unitful).

"""
    ParticleListEntry(; run_id=0, event_id, pdg, t, pos, ekin, dir)

One primary particle in a `ParticleListSource` CSV row. Multiple entries with
the same `(run_id, event_id)` are emitted together as primaries of one event.

- `t`   :: `TimeQ`             — emission time (any time unit; mm-equivalent ns)
- `pos` :: `NTuple{3,LengthQ}` — emission position (any length unit; → mm)
- `ekin`:: `EnergyQ`           — kinetic energy (any energy unit; → MeV)
- `dir` :: `NTuple{3,Float64}` — momentum direction (renormalised C++-side)
"""
Base.@kwdef struct ParticleListEntry
    run_id::Int = 0
    event_id::Int
    pdg::Int
    t::TimeQ
    pos::NTuple{3,LengthQ}
    ekin::EnergyQ
    dir::NTuple{3,Float64}
end

# Standard 11-column header. Mirrors the column order ParticleListSource.cc
# parses in the hasRunId branch. Emitted as the first row (no `#` prefix);
# the C++ reader detects header rows by their leading non-digit/non-`-`.
const _PARTICLE_LIST_HEADER =
    "run_id,event_id,pdg,t_ns,x_mm,y_mm,z_mm,ekin_MeV,dx,dy,dz"

# %.17g for floating-point fields, exactly as manifest.jl does — preserves a
# round-trip through C `strtod`. Integers use `string` (no formatting).
_pl_num(v::Float64) = @sprintf("%.17g", v)
_pl_num(v::Integer) = string(v)

function _write_entry(io::IO, e::ParticleListEntry)
    print(io,
        e.run_id, ",",
        e.event_id, ",",
        e.pdg, ",",
        _pl_num(ustrip(u"ns",  e.t)), ",",
        _pl_num(ustrip(u"mm",  e.pos[1])), ",",
        _pl_num(ustrip(u"mm",  e.pos[2])), ",",
        _pl_num(ustrip(u"mm",  e.pos[3])), ",",
        _pl_num(ustrip(u"MeV", e.ekin)), ",",
        _pl_num(e.dir[1]), ",",
        _pl_num(e.dir[2]), ",",
        _pl_num(e.dir[3]), "\n")
end

"""
    write_particle_list(path, entries; header=true) -> path
    write_particle_list(io,   entries; header=true)

Serialise `entries` (any iterable of `ParticleListEntry`) in the
11-column CSV format read by `ParticleListSource::LoadFile`. Lengths
convert to mm, times to ns, energies to MeV (Geant4-internal units); a
`%.17g` round-trip preserves the values exactly. With `header=true`
(default), the first row names the columns — recognised and skipped by
the C++ reader.

The path-based form returns `path` for chaining; the `IO` form returns
`nothing`.
"""
function write_particle_list(path::AbstractString, entries; header::Bool=true)
    open(io -> write_particle_list(io, entries; header), path, "w")
    return path
end

function write_particle_list(io::IO, entries; header::Bool=true)
    header && print(io, _PARTICLE_LIST_HEADER, "\n")
    for e in entries
        _write_entry(io, e)
    end
    return nothing
end
