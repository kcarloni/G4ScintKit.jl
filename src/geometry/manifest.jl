# manifest.jl
#
# Julia mirror of the C++ geometry manifest — the low-level, fully-explicit
# detector description that the C++ `DetectorConstruction::PlaceManifest`
# interpreter consumes (via `--manifest`).
#
# The structs mirror g4scintkit/include/Preparation/GeometryManifest.hh and the
# text serializer (`write_manifest`) reproduces g4scintkit/src/Preparation/
# ManifestFile.cc byte-for-byte (C `%.17g` number formatting). All values are
# plain Float64 in Geant4 internal units (mm, radian) — this is the machine
# representation. The user-facing, unitful layer is `DetectorSpec`
# (see detector_spec.jl); `build_manifest(::DetectorSpec)` produces a
# `GeometryManifest`.

# ---------------------------------------------------------------------------
#  Manifest data structures (mirror GeometryManifest.hh)
# ---------------------------------------------------------------------------

"""Abstract supertype for the four manifest placement-entry kinds."""
abstract type PlacementEntry end

"""A scintillator tile placement. Mirrors C++ `ScintEntry`."""
Base.@kwdef struct ScintEntry <: PlacementEntry
    name::String
    g4name::String = ""                                  # "" = don't call SetScintillatorName
    dims::NTuple{3,Float64}                              # full dimensions (x,y,z)
    pos::NTuple{3,Float64} = (0.0, 0.0, 0.0)             # translation in the mother volume
    rot::NTuple{9,Float64} = (1.0, 0.0, 0.0,             # row-major rotation matrix
                              0.0, 1.0, 0.0,
                              0.0, 0.0, 1.0)
    mother::String = "world"
    material_file::String
    sensitive::Bool = false
end

"""A single fibre segment (straight or circular arc). Mirrors C++ `FiberEntry`.
`stop` is the segment end point (`end` is reserved in Julia)."""
Base.@kwdef struct FiberEntry <: PlacementEntry
    name::String
    kind::String = "straight"                            # "straight" or "bent"
    mother::String = "world"
    start::NTuple{3,Float64}
    stop::NTuple{3,Float64}
    bend_angle::Float64 = 0.0                            # bent only
    bend_axis::NTuple{3,Float64} = (0.0, 0.0, 0.0)       # bent only
    material_file::String
    reference::String = ""                               # "" = SetFibreReferenceVolume not called
    glued::Bool = false
    glue_file::String = ""
    glue_profile::String = ""
    start_reflectivity::Float64 = NaN                    # NaN = SetFibreStartPointReflectivity not called
    end_reflectivity::Float64 = NaN                      # NaN = SetFibreEndPointReflectivity not called
    loop_id::Int = -1                                    # -1 = print length individually
end

"""A reflective wrapping around a scintillator tile. Mirrors C++ `WrapEntry`.
An empty `cut` means `PlaceManifest` auto-derives the cut candidates."""
Base.@kwdef struct WrapEntry <: PlacementEntry
    scint::String
    g4name::String = ""
    material_file::String
    cut::Vector{String} = String[]
end

"""A SiPM photon detector plus its optical coupling to a fibre. Mirrors C++ `SipmEntry`.
`model` names a g4sipm model alias (see [`SIPM_MODELS`](@ref)); an empty
`model` means "no g4sipm digitization — use the GODDESS photodetector"."""
Base.@kwdef struct SipmEntry <: PlacementEntry
    name::String
    ref_volume::String
    face_dir::NTuple{3,Float64}
    rel_pos::NTuple{3,Float64}
    edge_length::Float64
    fiber::String
    coupling_normal::NTuple{3,Float64}
    coupling_pos::NTuple{3,Float64}
    coupling_width::Float64
    fiber_is_base::Bool = true
    model::String = ""
end

"""The outer aluminum box + lead sheet, derived from the module bounding box.
Mirrors C++ `CasingSpec`."""
Base.@kwdef struct CasingSpec
    module_half_x::Float64 = 0.0
    module_min_y::Float64 = 0.0
    module_max_y::Float64 = 0.0
    module_half_z::Float64 = 0.0
    aluminum_thickness::Float64 = 0.0                    # <= 0 disables the aluminum box
    lead_thickness::Float64 = 0.0                        # <= 0 disables the lead sheet
    num_bars::Int = 0
    bar_width::Float64 = 0.0
    scinti_z::Float64 = 0.0
end

"""A complete, flat detector geometry: one ordered placement list + casing.
Mirrors C++ `GeometryManifest`."""
Base.@kwdef struct GeometryManifest
    setup_label::String = ""
    placements::Vector{PlacementEntry} = PlacementEntry[]
    casing::CasingSpec = CasingSpec()
end

# ---------------------------------------------------------------------------
#  Typed placement accessors
# ---------------------------------------------------------------------------

"""
    scintillators(m::GeometryManifest) -> Vector{ScintEntry}
    fibers(m::GeometryManifest)        -> Vector{FiberEntry}
    wraps(m::GeometryManifest)         -> Vector{WrapEntry}
    sipms(m::GeometryManifest)         -> Vector{SipmEntry}

The placements of one kind, in manifest (placement) order.
"""
scintillators(m::GeometryManifest) = ScintEntry[p for p in m.placements if p isa ScintEntry]
fibers(m::GeometryManifest)        = FiberEntry[p for p in m.placements if p isa FiberEntry]
wraps(m::GeometryManifest)         = WrapEntry[p  for p in m.placements if p isa WrapEntry]
sipms(m::GeometryManifest)         = SipmEntry[p  for p in m.placements if p isa SipmEntry]

# ---------------------------------------------------------------------------
#  Fibre lengths
# ---------------------------------------------------------------------------

"""
    fiber_length(f::FiberEntry) -> Quantity

Length of one fibre segment as a unitful quantity (mm). Straight segments
return the chord `norm(stop − start)`; bent (circular-arc) segments return
the arc length `r · bend_angle`, with the radius derived from chord and bend
angle as `r = chord / (2·sin(bend_angle / 2))`. Length is invariant under
the entry's reference frame, so no frame resolution is needed.
"""
function fiber_length(f::FiberEntry)
    chord = norm(SVector{3,Float64}(f.stop) - SVector{3,Float64}(f.start))
    L_mm = if f.kind == "straight"
        chord
    elseif f.kind == "bent"
        θ = f.bend_angle
        # Degenerate guard: a zero-angle bent segment is straight; avoid 0/0.
        abs(θ) < eps(Float64) ? chord : (chord / (2 * sin(θ / 2))) * θ
    else
        error("fiber_length: unknown fibre kind '$(f.kind)' " *
              "(expected 'straight' or 'bent')")
    end
    return L_mm * u"mm"
end

"""
    fiber_lengths(m::GeometryManifest) -> NamedTuple

Total length of each *physical* fibre in the manifest, as a `NamedTuple` of
unitful quantities (mm) in manifest order (first occurrence of each fibre
in the placement list). Segments sharing a non-negative `loop_id` are
summed under the field `:loop_<id>` (matching the C++ `PlaceManifest` print
format); single-segment fibres (`loop_id < 0`) appear under their own
`name`. Pair with [`fiber_length`](@ref) for the per-segment value.
"""
function fiber_lengths(m::GeometryManifest)
    L = typeof(1.0u"mm")
    keys = Symbol[]
    idx  = Dict{Symbol,Int}()
    vals = L[]
    for f in fibers(m)
        key = Symbol(f.loop_id < 0 ? f.name : "loop_$(f.loop_id)")
        if haskey(idx, key)
            vals[idx[key]] += fiber_length(f)
        else
            push!(keys, key); push!(vals, fiber_length(f))
            idx[key] = length(keys)
        end
    end
    return NamedTuple{Tuple(keys)}(Tuple(vals))
end

# ---------------------------------------------------------------------------
#  Text serializer (reproduces ManifestFile::Write byte-for-byte)
# ---------------------------------------------------------------------------

# C `std::snprintf(buf, sizeof(buf), "%.17g", v)` — full precision for an exact
# strtod round-trip.
function _numstr(v::Float64)
    isnan(v) && return "nan"                             # C prints "nan"; Julia "NaN"
    isinf(v) && return v < 0 ? "-inf" : "inf"
    return @sprintf("%.17g", v)
end

_vecstr(v::NTuple{3,Float64}) =
    string(_numstr(v[1]), ",", _numstr(v[2]), ",", _numstr(v[3]))

_rotstr(r::NTuple{9,Float64}) = join((_numstr(x) for x in r), ",")

_liststr(xs::Vector{String}) = join(xs, ",")

_boolstr(b::Bool) = b ? "1" : "0"

# The flat key=value format is space-delimited, so a value must not contain a
# space. Mirrors ManifestFile's requireNoSpace: fail fast rather than emit a
# line that the C++ reader cannot parse back.
function _require_no_space(value::AbstractString, field::AbstractString)
    if occursin(' ', value)
        error("manifest: '$field' value contains a space, which the flat " *
              "manifest format cannot represent: '$value'")
    end
end

function _write_entry(io::IO, s::ScintEntry)
    _require_no_space(s.name, "name")
    _require_no_space(s.g4name, "g4name")
    _require_no_space(s.mother, "mother")
    _require_no_space(s.material_file, "material")
    print(io, "SCINT",
        " name=", s.name,
        " g4name=", s.g4name,
        " dims=", _vecstr(s.dims),
        " pos=", _vecstr(s.pos),
        " rot=", _rotstr(s.rot),
        " mother=", s.mother,
        " material=", s.material_file,
        " sensitive=", _boolstr(s.sensitive),
        "\n")
end

function _write_entry(io::IO, f::FiberEntry)
    _require_no_space(f.name, "name")
    _require_no_space(f.mother, "mother")
    _require_no_space(f.reference, "reference")
    _require_no_space(f.material_file, "material")
    _require_no_space(f.glue_file, "glue_file")
    _require_no_space(f.glue_profile, "glue_profile")
    print(io, "FIBER",
        " name=", f.name,
        " kind=", f.kind,
        " mother=", f.mother,
        " start=", _vecstr(f.start),
        " end=", _vecstr(f.stop),
        " bend_angle=", _numstr(f.bend_angle),
        " bend_axis=", _vecstr(f.bend_axis),
        " material=", f.material_file,
        " reference=", f.reference,
        " glued=", _boolstr(f.glued),
        " glue_file=", f.glue_file,
        " glue_profile=", f.glue_profile)
    # Conditionally emitted: defaults are NaN ("not called"); omitting at
    # default keeps existing baselines stable while letting new designs
    # surface the knob. The C++ parser tolerates a missing `start_refl`
    # (optNum -> NaN). `end_refl` stays unconditional for backward
    # compatibility with existing files that always carry it.
    if !isnan(f.start_reflectivity)
        print(io, " start_refl=", _numstr(f.start_reflectivity))
    end
    print(io,
        " end_refl=", _numstr(f.end_reflectivity),
        " loop_id=", string(f.loop_id),
        "\n")
end

function _write_entry(io::IO, w::WrapEntry)
    _require_no_space(w.scint, "scint")
    _require_no_space(w.g4name, "g4name")
    _require_no_space(w.material_file, "material")
    for c in w.cut
        _require_no_space(c, "cut")
    end
    print(io, "WRAP",
        " scint=", w.scint,
        " g4name=", w.g4name,
        " material=", w.material_file,
        " cut=", _liststr(w.cut),
        "\n")
end

function _write_entry(io::IO, sp::SipmEntry)
    _require_no_space(sp.name, "name")
    _require_no_space(sp.ref_volume, "ref_volume")
    _require_no_space(sp.fiber, "fiber")
    _require_no_space(sp.model, "model")
    print(io, "SIPM",
        " name=", sp.name,
        " ref_volume=", sp.ref_volume,
        " face_dir=", _vecstr(sp.face_dir),
        " rel_pos=", _vecstr(sp.rel_pos),
        " edge_length=", _numstr(sp.edge_length),
        " fiber=", sp.fiber,
        " coupling_normal=", _vecstr(sp.coupling_normal),
        " coupling_pos=", _vecstr(sp.coupling_pos),
        " coupling_width=", _numstr(sp.coupling_width),
        " fiber_is_base=", _boolstr(sp.fiber_is_base))
    # `model` is optional: emitted only when non-empty so existing manifests
    # (and the C++ reader pre-C2) stay backward-compatible. Empty = no g4sipm.
    isempty(sp.model) || print(io, " model=", sp.model)
    print(io, "\n")
end

function _write_entry(io::IO, c::CasingSpec)
    print(io, "CASING",
        " module_half_x=", _numstr(c.module_half_x),
        " module_min_y=", _numstr(c.module_min_y),
        " module_max_y=", _numstr(c.module_max_y),
        " module_half_z=", _numstr(c.module_half_z),
        " aluminum_thickness=", _numstr(c.aluminum_thickness),
        " lead_thickness=", _numstr(c.lead_thickness),
        " num_bars=", string(c.num_bars),
        " bar_width=", _numstr(c.bar_width),
        " scinti_z=", _numstr(c.scinti_z),
        "\n")
end

"""
    write_manifest(path::AbstractString, manifest::GeometryManifest) -> path
    write_manifest(io::IO, manifest::GeometryManifest)

Serialise `manifest` in the exact `ManifestFile` text format that the C++
`DetectorConstruction` reads via `--manifest`. The output is byte-for-byte
identical to a C++ `ManifestFile::Write` of the same geometry.

The path-based form returns `path` (convenient for chaining); the `IO` form
returns `nothing`.
"""
function write_manifest(path::AbstractString, m::GeometryManifest)
    open(io -> write_manifest(io, m), path, "w")
    return path
end

function write_manifest(io::IO, m::GeometryManifest)
    print(io, "# G4ScintKit geometry manifest\n")
    print(io, "# units: Geant4 internal (mm, radian)\n")
    print(io, "# SCINT/FIBER/WRAP/SIPM lines, in file order, are the placement order\n")
    isempty(m.setup_label) || print(io, "SETUP ", m.setup_label, "\n")
    for p in m.placements
        _write_entry(io, p)
    end
    _write_entry(io, m.casing)
    return nothing
end
