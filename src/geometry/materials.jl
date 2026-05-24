# materials.jl
#
# Material-file resolution for detector designs: the abbreviation maps (ported
# from run.sh's resolve_* tables), the lookup that turns an abbreviation into an
# absolute GODDESS `.properties` path, the per-design `ResolvedMaterials` bundle,
# and the fibre-geometry reader (`fiber_cross_section`).
#
# A `DetectorSpec` carries material abbreviations (`scint_material`, …); this
# layer turns them into the concrete paths the C++ side reads.

# ---------------------------------------------------------------------------
#  Material-file abbreviation maps (port run.sh's resolve_* tables)
# ---------------------------------------------------------------------------

const _SCINT_FILES = Dict(
    "fermilab"   => "Fermilab_scintillator.properties",
    "bc404"      => "Saint-Gobain_BC-404.properties",
    "bc408"      => "Saint-Gobain_BC-408.properties",
    "bc452-2pb"  => "Saint-Gobain_BC-452_2perCentLead.properties",
    "bc452-5pb"  => "Saint-Gobain_BC-452_5perCentLead.properties",
    "bc452-10pb" => "Saint-Gobain_BC-452_10perCentLead.properties",
)

const _WRAP_FILES = Dict(
    "tio2"   => "Wrapping_TiO2.properties",
    "teflon" => "Wrapping_Teflon.properties",
    "alu"    => "Wrapping_Aluminum.properties",
    "bc620"  => "Wrapping_Saint-Gobain_BC-620.properties",
    "tyvek"  => "Wrapping_Tyvek.properties",
)

const _FIBER_FILES = Dict(
    "y11-200-r1"    => "Kuraray_Y11-200_round_1mm.properties",
    "y11-200-r1-sc" => "Kuraray_Y11-200_round_1mm_singleClad.properties",
    "y11-200-s1-sc" => "Kuraray_Y11-200_square_1mm_singleClad.properties",
    "y11-300-r1"    => "Kuraray_Y11-300_round_1mm.properties",
    "bcf10-ms1"     => "Saint-Gobain_BCF-10_multi_square_1mm.properties",
    "bcf92-r1"      => "Saint-Gobain_BCF-92_round_1mm.properties",
    "bcf92-r1-sc"   => "Saint-Gobain_BCF-92_round_1mm_singleClad.properties",
    "bcf92-r2"      => "Saint-Gobain_BCF-92_round_2mm.properties",
    "bcf92-r2-sc"   => "Saint-Gobain_BCF-92_round_2mm_singleClad.properties",
    "bcf92-q1"      => "Saint-Gobain_BCF-92_quadratic_1mm.properties",
    "bcf92-q1-sc"   => "Saint-Gobain_BCF-92_quadratic_1mm_singleClad.properties",
    "bcf92-q2"      => "Saint-Gobain_BCF-92_quadratic_2mm.properties",
    "bcf92-q2-sc"   => "Saint-Gobain_BCF-92_quadratic_2mm_singleClad.properties",
    "bcf98-r1"      => "Saint-Gobain_BCF-98_round_1mm.properties",
    "bcf98-r1-sc"   => "Saint-Gobain_BCF-98_round_1mm_singleClad.properties",
    "bcf98-r2"      => "Saint-Gobain_BCF-98_round_2mm.properties",
    "bcf98-r2-sc"   => "Saint-Gobain_BCF-98_round_2mm_singleClad.properties",
    "bcf98-q1"      => "Saint-Gobain_BCF-98_quadratic_1mm.properties",
    "bcf98-q1-sc"   => "Saint-Gobain_BCF-98_quadratic_1mm_singleClad.properties",
    "bcf98-q2"      => "Saint-Gobain_BCF-98_quadratic_2mm.properties",
    "bcf98-q2-sc"   => "Saint-Gobain_BCF-98_quadratic_2mm_singleClad.properties",
    "eo534b"        => "EO-534B.properties",
)

const _CEMENT_FILES = Dict(
    "bc600"  => "Saint-Gobain_BC-600.properties",
    "air"    => "air.properties",
    "air1mm" => "air_1mm_fiber.properties",
)

"""Default GODDESS package root — `goddess-package/` sits beside `G4ScintKit.jl/`."""
# this file lives at src/geometry/ — three levels below the G4ScintKit project
# root, where goddess-package/ sits beside the G4ScintKit.jl package.
_default_goddess_root() =
    normpath(joinpath(@__DIR__, "..", "..", "..", "goddess-package"))

# Resolve a material abbreviation to the absolute .properties path, matching
# the path run.sh builds: <goddess>/source/MaterialProperties/<subdir>/<file>.
function _resolve_material(goddess_root::AbstractString, subdir::AbstractString,
                           table::Dict{String,String}, key::AbstractString)
    haskey(table, key) || error("unknown material abbreviation '$key' (known: " *
        join(sort(collect(keys(table))), ", ") * ")")
    return joinpath(goddess_root, "source", "MaterialProperties", subdir, table[key])
end

# ---------------------------------------------------------------------------
#  Resolved materials
# ---------------------------------------------------------------------------

"""
    ResolvedMaterials(scint, wrap, wls, cement)
    ResolvedMaterials(spec)

The four `.properties` file paths a detector design needs, resolved once. The
`spec` form reads the standard abbreviation fields (`goddess_root`,
`scint_material`, `wrap_material`, `wls_material`, `cement_material`) and
resolves each against the GODDESS package — the same paths `run.sh` builds.
"""
struct ResolvedMaterials
    scint::String
    wrap::String
    wls::String
    cement::String
end

function ResolvedMaterials(spec)
    root = spec.goddess_root
    return ResolvedMaterials(
        _resolve_material(root, "Scintillator",  _SCINT_FILES,  spec.scint_material),
        _resolve_material(root, "Scintillator",  _WRAP_FILES,   spec.wrap_material),
        _resolve_material(root, "Fibre",         _FIBER_FILES,  spec.wls_material),
        _resolve_material(root, "OpticalCement", _CEMENT_FILES, spec.cement_material),
    )
end

# ---------------------------------------------------------------------------
#  Fibre geometry from the material file
# ---------------------------------------------------------------------------

# Read the raw value of `key` from a GODDESS `.properties` file. Lines are
# `key: value` or `key = value`; `#` starts a comment.
function _read_property(file::AbstractString, key::AbstractString)
    for line in eachline(file)
        ci = findfirst(==('#'), line)
        content = ci === nothing ? line : line[1:prevind(line, ci)]
        si = findfirst(c -> c == ':' || c == '=', content)
        si === nothing && continue
        strip(content[1:prevind(content, si)]) == key || continue
        return strip(content[nextind(content, si):end])
    end
    error("property '$key' not found in $file")
end

# Parse a GODDESS length value (e.g. "0.50 * mm") to plain Float64 millimetres.
function _parse_length_mm(s::AbstractString)
    parts = split(s, '*')
    value = parse(Float64, strip(parts[1]))
    length(parts) == 1 && return value                   # bare number: assume mm
    unit = strip(parts[2])
    factor = unit == "mm" ? 1.0   :
             unit == "cm" ? 10.0  :
             unit == "m"  ? 1000.0 :
             error("unsupported length unit '$unit' (expected mm, cm or m)")
    return value * factor
end

"""
    fiber_cross_section(properties_file) -> (; profile::Symbol, width::Float64)

Read a fibre's transverse geometry from its GODDESS `.properties` file.
`profile` is `:round` or `:quadratic`; `width` is the full transverse size in
millimetres — `2 × radius` for a round fibre, `edge_length` for a quadratic
one. Lets a `build_manifest` derive the fibre diameter from the chosen
`wls_material` rather than hard-coding it.
"""
function fiber_cross_section(properties_file::AbstractString)
    isfile(properties_file) ||
        error("fiber_cross_section: no such file: $properties_file")
    profile = Symbol(_read_property(properties_file, "profile"))
    if profile === :round
        width = 2 * _parse_length_mm(_read_property(properties_file, "radius")) * u"mm"
    elseif profile === :quadratic
        width = _parse_length_mm(_read_property(properties_file, "edge_length")) * u"mm"
    else
        error("fiber_cross_section: unknown profile '$profile' in " *
              "$properties_file (expected round or quadratic)")
    end
    return (; profile = profile, width = width)
end
