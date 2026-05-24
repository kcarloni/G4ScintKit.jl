# detector_spec.jl
#
# The user-facing, high-level detector-design layer. A `DetectorSpec` is a
# typed, unitful struct whose fields are the design knobs (the "type-2 knob
# schema" of the architecture plan). `build_manifest(::DetectorSpec)` runs the
# layout trig and emits a low-level `GeometryManifest` (see manifest.jl).
#
# `build_manifest` mirrors the C++ `DetectorConstruction::BuildManifest_B1/B2`
# expression-by-expression. User-facing values are unitful; the trig is done in
# plain Float64 in Geant4 internal units (mm, radian) so the result is
# byte-compatible with the C++ reference manifest. Material abbreviations are
# resolved by materials.jl.

"""Length quantity type for `DetectorSpec` fields — forces inputs to a length
dimension (any length unit; converted to mm on construction)."""
const LengthQ  = typeof(1.0u"mm")

"""Energy quantity (e.g. column types on the readback side)."""
const EnergyQ  = typeof(1.0u"MeV")

"""Time quantity (e.g. column types on the readback side)."""
const TimeQ    = typeof(1.0u"ns")

"""Voltage quantity (e.g. column types on the readback side)."""
const VoltageQ = typeof(1.0u"V")

# ---------------------------------------------------------------------------
#  Detector designs
# ---------------------------------------------------------------------------

"""Abstract supertype for typed detector designs. Each concrete subtype's
fields are the design's knob schema; `build_manifest` turns one into a
`GeometryManifest`."""
abstract type DetectorSpec end

# ---------------------------------------------------------------------------
#  build_manifest — layout trig (mirrors C++ BuildManifest_B1 / BuildManifest_B2)
# ---------------------------------------------------------------------------

"""
    build_manifest(spec::DetectorSpec) -> GeometryManifest

Compute the flat geometry manifest for a detector design. The layout trig runs
in plain Float64 in Geant4 internal units (mm, radian), mirroring the C++
`BuildManifest_*` expression-by-expression so the result is byte-compatible
with the C++ reference manifest. Pair with [`write_manifest`](@ref) to emit a
file the C++ side can consume via `--manifest`.
"""
function build_manifest end

"""
    strip_units(spec::DetectorSpec) -> NamedTuple

A unit-free view of `spec` for use inside `build_manifest`: every length field
is converted to plain `Float64` millimetres (Geant4 internal units); all other
fields (counts, material names, paths) pass through unchanged. Collapses the
per-field `ustrip(u"mm", spec.x)` boilerplate into one call:

    g = strip_units(spec)
    g.scint_width        # Float64, in mm
"""
function strip_units(spec::DetectorSpec)
    names = fieldnames(typeof(spec))
    return NamedTuple{names}(map(n -> _strip_field(getfield(spec, n)), names))
end

# A length quantity becomes Float64 millimetres; any other field passes through.
_strip_field(v) = v isa LengthQ ? ustrip(u"mm", v) : v

"""
    write_manifest(path::AbstractString, spec::DetectorSpec) -> path

Build the geometry manifest for `spec` and serialise it to `path` in one step —
equivalent to `write_manifest(path, build_manifest(spec))`. `build_manifest`
already runs the pre-flight geometry check.
"""
write_manifest(path::AbstractString, spec::DetectorSpec) =
    write_manifest(path, build_manifest(spec))
