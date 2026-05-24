# columns.jl
#
# Convention for HDF5 column names produced by the C++ side: a trailing
# `_<unit>` suffix annotates the unit of a numeric column. This file converts
# raw HDF5 columns into "cleaned" form: the suffix is stripped from the symbol
# and the corresponding Unitful unit is attached to the values.

# Suffix → unit mapping. Order matters: longer suffixes are tried first so
# that "_MeV" wins over "_V".
const _UNIT_SUFFIXES = (
    # energies
    ("_MeV", u"MeV"),
    ("_GeV", u"GeV"),
    ("_keV", u"keV"),
    ("_eV",  u"eV"),
    # times
    ("_ns",  u"ns"),
    ("_us",  u"μs"),
    ("_μs",  u"μs"),
    ("_ms",  u"ms"),
    # voltages
    ("_mV",  u"mV"),
    ("_V",   u"V"),
    # lengths
    ("_mm",  u"mm"),
    ("_cm",  u"cm"),
    ("_m",   u"m"),
    # bare seconds last (don't shadow "_ns" / "_ms" / "_us")
    ("_s",   u"s"),
)

# (sym, value) -> (clean_sym, value_with_units). Pass-through if no suffix matches.
function _strip_and_attach(sym::Symbol, val)
    s = String(sym)
    for (suf, unit) in _UNIT_SUFFIXES
        if endswith(s, suf)
            clean = Symbol(s[1:end-length(suf)])
            return clean, val .* unit
        end
    end
    return sym, val
end

# Rename HDF5's `g4event_id` to the more idiomatic `event_id`. Other suffix-less
# columns are passed through unchanged.
_rename(sym::Symbol) = sym === :g4event_id ? :event_id : sym

"""
    _clean_columns(nt::NamedTuple) -> NamedTuple

Apply the kit's column-cleanup conventions:
- strip a trailing unit suffix (`_mm`, `_MeV`, `_ns`, `_V`, …) from each name
  and attach the corresponding `Unitful` unit to its values;
- rename `g4event_id` → `event_id`.

Non-suffixed columns (ints, dimensionless doubles) pass through unchanged.
"""
function _clean_columns(nt::NamedTuple)
    ks = keys(nt)
    vs = values(nt)
    new_pairs = Pair{Symbol, Any}[]
    for (k, v) in zip(ks, vs)
        ck, cv = _strip_and_attach(k, v)
        push!(new_pairs, _rename(ck) => cv)
    end
    NamedTuple{Tuple(first.(new_pairs))}(Tuple(last.(new_pairs)))
end
