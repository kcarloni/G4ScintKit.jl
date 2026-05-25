# sipms.jl
#
# SiPM model registry: the table of model aliases that may appear in a
# `SipmEntry.model` field, the geometric dimensions (cell pitch and cell
# count → die edge length) for each, and helpers to look them up.
#
# Two backing routes (recorded as `factory`):
#   - `:builtin`     — g4sipm ships a hardcoded C++ subclass with bespoke
#                      physics (PDE curve, voltage trace). C++ side will
#                      construct it via `G4SipmModelFactory::create<…>`.
#   - `:config_file` — a generic `G4SipmConfigFileModel` constructed from a
#                      `.properties` file under `<g4sipm>/sample/resources/`.
#                      The `properties_file` field names the file basename.
#
# Geometric dimensions live here as a hand-maintained table because g4sipm's
# `:builtin` models bake their numbers into C++ source (e.g.
# `g4sipm/g4sipm/src/model/impl/HamamatsuS12573100C.cc`) rather than into any
# parseable config file. The values are stable for a vendored g4sipm submodule;
# update this table only if g4sipm changes its model classes. The two models
# that *do* have matching `.properties` files (s10362-11-100c, s10362-33-050c)
# are cross-checked by the test suite.
#
# `edge_length = cell_pitch * sqrt(n_cells)` (cells are arranged in a square
# grid; cellPitch is the centre-to-centre cell spacing).

"""
    SipmModelInfo(; cell_pitch, n_cells, edge_length, factory, properties_file="")

Geometry + dispatch info for one g4sipm SiPM model.

- `cell_pitch`     :: `LengthQ` — microcell centre-to-centre spacing
- `n_cells`        :: `Int`     — total number of microcells in the die
- `edge_length`    :: `LengthQ` — die edge length (`cell_pitch * sqrt(n_cells)`)
- `factory`        :: `Symbol`  — `:builtin` or `:config_file`
- `properties_file`:: `String`  — `.properties` basename for `:config_file`
                                  models (empty for `:builtin`)
"""
Base.@kwdef struct SipmModelInfo
    cell_pitch::LengthQ
    n_cells::Int
    edge_length::LengthQ
    factory::Symbol
    properties_file::String = ""
end

# Construct a `:builtin` entry, deriving edge_length from cell_pitch * sqrt(n_cells).
_builtin(cell_pitch::LengthQ, n_cells::Integer) = SipmModelInfo(
    cell_pitch  = cell_pitch,
    n_cells     = Int(n_cells),
    edge_length = cell_pitch * sqrt(n_cells),
    factory     = :builtin)

"""
    SIPM_MODELS :: Dict{String, SipmModelInfo}

The known g4sipm SiPM model aliases. Mirrors run.sh's old
`resolve_sipmmodel` switch and the C++ `G4SipmModelFactory::create<…>`
entrypoints. Use `sipm_model_info(alias)` rather than indexing directly
so unknown aliases produce a clear error.
"""
const SIPM_MODELS = Dict{String, SipmModelInfo}(
    # g4sipm generic default — implemented as G4SipmGenericSipmModel
    # (g4sipm/g4sipm/src/model/impl/G4SipmGenericSipmModel.cc).
    "generic"                  => _builtin(0.1u"mm",   100),

    # Hardcoded Hamamatsu factory classes in g4sipm. Cell pitch and count are
    # baked into their respective .cc files; see header of this file.
    "hamamatsu-s10362-11-100c" => _builtin(0.1u"mm",    100),  # 1×1 mm
    "hamamatsu-s10362-33-100c" => _builtin(0.1u"mm",    900),  # 3×3 mm
    "hamamatsu-s10362-33-050c" => _builtin(0.05u"mm", 3600),   # 3×3 mm
    "hamamatsu-s12651-050"     => _builtin(0.05u"mm",  400),   # 1×1 mm
    "hamamatsu-s12573-100c"    => _builtin(0.1u"mm",  3600),   # 6×6 mm
    "hamamatsu-s12573-100x"    => _builtin(0.1u"mm",  3600),   # 6×6 mm
)

"""
    sipm_model_info(alias) -> SipmModelInfo

Look up `alias` in [`SIPM_MODELS`](@ref); error with a sorted list of valid
aliases otherwise. Use the empty string to mean "no g4sipm model" (GODDESS
photodetector); that case should be handled by the caller, not passed in.
"""
function sipm_model_info(alias::AbstractString)
    isempty(alias) && error("sipm_model_info: empty model alias " *
                            "(empty means \"no g4sipm model\" — handle that " *
                            "case at the call site)")
    haskey(SIPM_MODELS, alias) && return SIPM_MODELS[alias]
    error("sipm_model_info: unknown SiPM model '$alias' (known: " *
          join(sort(collect(keys(SIPM_MODELS))), ", ") * ")")
end

"""
    sipm_edge_length(alias) -> LengthQ

Die edge length of the SiPM model `alias`. Shorthand for
`sipm_model_info(alias).edge_length`.
"""
sipm_edge_length(alias::AbstractString) = sipm_model_info(alias).edge_length
