# g4coordinate.jl
#
# Three frame-bound geometric primitives, each carrying the name of the
# reference volume it is expressed in:
#
#   * `G4Coordinate` — a *position*. Two positions cannot be added (positions
#     have no meaningful sum); subtracting two positions yields a `G4Vector`,
#     the displacement between them.
#   * `G4Vector`     — a *free displacement* (frame-bound 3-vector, any
#     magnitude). Vectors add to vectors, and translate positions
#     (Coord + Vec → Coord). The result of `scalar * G4Direction` is a
#     `G4Vector`.
#   * `G4Direction`  — a *unit-magnitude direction*. Auto-normalised at
#     construction. Use for tangents (`G4DirectedPoint.direction`), face
#     normals, working-plane normals — anywhere the magnitude is meaningless.
#
# This mirrors the standard affine-space split: points + free vectors, with a
# unit-vector subtype where the magnitude is fixed. Same-frame arithmetic is
# always frame-checked.
#
# A bare 3-tuple is ambiguous: every position in the geometry manifest is
# relative to *some* parent volume. `G4Coordinate` carries that frame with the
# position so a coordinate is self-describing; `G4Vector` and `G4Direction` do
# the same for displacements and tangents.

# ---------------------------------------------------------------------------
#  Shared coercion helpers
# ---------------------------------------------------------------------------

# Coerce one length value to plain Float64 millimetres: a real number is taken
# as-is (already mm); a unitful length is converted; a non-length quantity
# makes `ustrip(u"mm", …)` throw a dimension-mismatch error. Shared across
# every user-facing length surface.
_to_mm(v::Real) = Float64(v)
_to_mm(v)       = Float64(ustrip(u"mm", v))   # force Float (ustrip keeps Int when exact)

# Coerce a 3-vector-like to NTuple{3,Float64} in mm.
_to_mm3(v) = (_to_mm(v[1]), _to_mm(v[2]), _to_mm(v[3]))

# ---------------------------------------------------------------------------
#  Volume-reference name coercion
# ---------------------------------------------------------------------------
#
# A `ref` may be supplied as a name string, a `Symbol`, or a placed entry —
# the entry is coerced to its `.name`. Used by all three primitive types and
# by `add_*!` kwargs that name another placement (`mother`, `scint`, `fiber`).
# The `_VolRef` Union is defined *after* the placement entry types so it can
# include them; the entries themselves are declared in manifest.jl.

const _VolRef = Union{AbstractString, Symbol, ScintEntry, FiberEntry, SipmEntry}

_volname(s::AbstractString) = String(s)
_volname(s::Symbol)         = String(s)
_volname(p::ScintEntry)     = p.name
_volname(p::FiberEntry)     = p.name
_volname(p::SipmEntry)      = p.name

# ---------------------------------------------------------------------------
#  G4Coordinate — a position bound to its reference frame
# ---------------------------------------------------------------------------

"""
    G4Coordinate(pos, ref="world")
    G4Coordinate(x, y, z, ref="world")

A 3D *position* `pos` (mm) together with `ref`, the name of the reference
volume it is expressed relative to — either `"world"` or the name of a placed
volume. `pos` may be any indexable 3-element value; each element may be a
plain real number (taken as mm) or a unitful length (converted to mm). `ref`
may be a string, a `Symbol`, or a placed entry — coerced to its `.name`.

A `G4Coordinate` is a position, not a free vector: two positions cannot be
added. Use [`G4Vector`](@ref) for a displacement (`Coord ± Vec → Coord`;
`Coord - Coord → Vec`).

See `docs/placement_rules.md` for which volume each placement kind measures
its coordinates against.
"""
struct G4Coordinate
    pos::NTuple{3,Float64}
    ref::String
end

# Typed outer (not `(::Any, ::Any)`) so the auto-generated outer constructor —
# Julia's fallback `convert`-based one — is not overwritten.
G4Coordinate(pos, ref::_VolRef) = G4Coordinate(_to_mm3(pos), _volname(ref))
G4Coordinate(pos)               = G4Coordinate(pos, "world")
G4Coordinate(x, y, z, ref::_VolRef = "world") =
    G4Coordinate((_to_mm(x), _to_mm(y), _to_mm(z)), _volname(ref))

Base.show(io::IO, c::G4Coordinate) =
    print(io, "G4Coordinate(", c.pos, ", \"", c.ref, "\")")

# ---------------------------------------------------------------------------
#  G4Vector — a free displacement bound to its reference frame
# ---------------------------------------------------------------------------

"""
    G4Vector(vec, ref="world")
    G4Vector(x, y, z, ref="world")

A 3D *free displacement* `vec` (mm) together with `ref`, the name of the
reference frame it is expressed in. Construction mirrors [`G4Coordinate`](@ref)
(tuple, three scalars, mixed plain/unitful lengths, name/symbol/entry `ref`).

The vector is a displacement, not a position: `Coord + Vec → Coord`,
`Vec + Vec → Vec`, and `Coord - Coord → Vec`. Multiplying a [`G4Direction`](@ref)
by a length scalar also yields a `G4Vector`.
"""
struct G4Vector
    vec::NTuple{3,Float64}
    ref::String
end

G4Vector(vec, ref::_VolRef) = G4Vector(_to_mm3(vec), _volname(ref))
G4Vector(vec)               = G4Vector(vec, "world")
G4Vector(x, y, z, ref::_VolRef = "world") =
    G4Vector((_to_mm(x), _to_mm(y), _to_mm(z)), _volname(ref))

Base.show(io::IO, v::G4Vector) =
    print(io, "G4Vector(", v.vec, ", \"", v.ref, "\")")

# ---------------------------------------------------------------------------
#  G4Direction — a unit tangent bound to its reference frame
# ---------------------------------------------------------------------------

"""
    G4Direction(vec, ref="world")
    G4Direction(x, y, z, ref="world")
    G4Direction(v::G4Vector)

A 3D *unit tangent direction* `vec` together with its reference frame. The
vector is auto-normalised; a zero vector is an error. Construction mirrors
[`G4Vector`](@ref) (tuple, three scalars; lengths are unitless so do not pass
unitful values). Passing a `G4Vector` through normalises it and carries its
frame.

Use `G4Direction` wherever a tangent is meaningful only up to magnitude —
`G4DirectedPoint.direction`, working-plane normals, SiPM face directions.
"""
struct G4Direction
    vec::SVector{3,Float64}
    ref::String
    function G4Direction(vec, ref::_VolRef)
        v = SVector{3,Float64}(Float64(vec[1]), Float64(vec[2]), Float64(vec[3]))
        n = norm(v)
        n > 0 || error("G4Direction: direction must be a nonzero vector")
        return new(v / n, _volname(ref))
    end
end

G4Direction(vec) = G4Direction(vec, "world")
G4Direction(x::Real, y::Real, z::Real, ref::_VolRef = "world") =
    G4Direction((x, y, z), ref)
# idempotent: passing a G4Direction through returns it unchanged
G4Direction(d::G4Direction) = d
function G4Direction(d::G4Direction, ref::_VolRef)
    target = _volname(ref)
    d.ref == target || error(
        "G4Direction: direction is in frame '$(d.ref)' but ref='$target' " *
        "was requested — both must match")
    return d
end
# A G4Vector is the natural displacement form; normalise it.
G4Direction(v::G4Vector) = G4Direction(v.vec, v.ref)
function G4Direction(v::G4Vector, ref::_VolRef)
    target = _volname(ref)
    v.ref == target || error(
        "G4Direction: vector is in frame '$(v.ref)' but ref='$target' was " *
        "requested — both must match")
    return G4Direction(v.vec, v.ref)
end

Base.show(io::IO, d::G4Direction) =
    print(io, "G4Direction(", Tuple(d.vec), ", \"", d.ref, "\")")

# ---------------------------------------------------------------------------
#  Affine arithmetic
# ---------------------------------------------------------------------------
#
# Same-frame is required for every binary operation. A bare 3-vector-like
# (tuple, SVector, AbstractVector) is treated as a `G4Vector` in the operand's
# own frame — the common convenience for hand-coded offsets.

const _BareVec = Union{Tuple,AbstractVector}

function _check_same_ref(a_ref::AbstractString, b_ref::AbstractString, op::AbstractString)
    a_ref == b_ref || error(
        "$op: operands are in different frames " *
        "('$a_ref' and '$b_ref') — frame-bound arithmetic requires a " *
        "shared reference frame")
end

# --- G4Vector with G4Vector ------------------------------------------------

function Base.:+(a::G4Vector, b::G4Vector)
    _check_same_ref(a.ref, b.ref, "G4Vector +")
    return G4Vector(a.vec .+ b.vec, a.ref)
end
function Base.:-(a::G4Vector, b::G4Vector)
    _check_same_ref(a.ref, b.ref, "G4Vector -")
    return G4Vector(a.vec .- b.vec, a.ref)
end
Base.:-(a::G4Vector) = G4Vector(.- a.vec, a.ref)

# Bare-vector convenience: treat as a G4Vector in the operand's frame.
Base.:+(a::G4Vector, v::_BareVec) = G4Vector(a.vec .+ _to_mm3(v), a.ref)
Base.:+(v::_BareVec, a::G4Vector) = a + v
Base.:-(a::G4Vector, v::_BareVec) = G4Vector(a.vec .- _to_mm3(v), a.ref)

# --- G4Coordinate with G4Vector / G4Coordinate -----------------------------

function Base.:+(c::G4Coordinate, v::G4Vector)
    _check_same_ref(c.ref, v.ref, "G4Coordinate + G4Vector")
    return G4Coordinate(c.pos .+ v.vec, c.ref)
end
Base.:+(v::G4Vector, c::G4Coordinate) = c + v

function Base.:-(c::G4Coordinate, v::G4Vector)
    _check_same_ref(c.ref, v.ref, "G4Coordinate - G4Vector")
    return G4Coordinate(c.pos .- v.vec, c.ref)
end

# Difference of two positions is the displacement between them.
function Base.:-(a::G4Coordinate, b::G4Coordinate)
    _check_same_ref(a.ref, b.ref, "G4Coordinate - G4Coordinate")
    return G4Vector(a.pos .- b.pos, a.ref)
end

# Bare-vector convenience: treat as a G4Vector in the coord's frame.
Base.:+(c::G4Coordinate, v::_BareVec) = G4Coordinate(c.pos .+ _to_mm3(v), c.ref)
Base.:+(v::_BareVec, c::G4Coordinate) = c + v
Base.:-(c::G4Coordinate, v::_BareVec) = G4Coordinate(c.pos .- _to_mm3(v), c.ref)

# Positions don't add and don't negate. Catch the common mistakes with a
# clear error so failures point at the conceptual issue rather than a missing
# method.
Base.:+(::G4Coordinate, ::G4Coordinate) = error(
    "G4Coordinate + G4Coordinate is not defined — positions don't add. " *
    "Use a G4Vector for the displacement, or take `b - a` if you want the " *
    "vector from `a` to `b`.")
Base.:-(::G4Coordinate) = error(
    "-G4Coordinate is not defined — a position has no sign. " *
    "Negate a G4Vector instead.")

"""
    midpoint(a::G4Coordinate, b::G4Coordinate) -> G4Coordinate

The affine midpoint of two positions, computed as `a + 0.5 * (b - a)` so that
the position/displacement split is respected (positions don't add, but their
half-difference is a `G4Vector` that can translate `a`). Both coordinates must
share a reference frame; the result is in that frame.
"""
midpoint(a::G4Coordinate, b::G4Coordinate) = a + 0.5 * (b - a)

# --- scalar × direction / vector -------------------------------------------

# A length scalar (Real mm or unitful length) times a unit direction produces a
# displacement carrying the direction's frame. Lets natural idioms like
# `sipm_pos + sipm_coupling_width * sipm_facing_dir` work.
Base.:*(s::Number, d::G4Direction) =
    G4Vector(Tuple(_to_mm(s) .* d.vec), d.ref)
Base.:*(d::G4Direction, s::Number) = s * d

Base.:*(s::Number, v::G4Vector) = G4Vector(_to_mm(s) .* v.vec, v.ref)
Base.:*(v::G4Vector, s::Number) = s * v

# --- named axis accessors --------------------------------------------------
#
# `c.x` / `c.y` / `c.z` pull the corresponding component out of the underlying
# `pos` / `vec` triple. The real fields (`pos`/`vec` and `ref`) still resolve
# via `getfield` as usual. Deliberately scoped to these three symbols — no
# iteration, broadcasting, or `getindex` — so the affine-space split (no
# `coord + coord`, no silent frame-stripping) is preserved.
function Base.getproperty(c::G4Coordinate, s::Symbol)
    s === :x ? getfield(c, :pos)[1] :
    s === :y ? getfield(c, :pos)[2] :
    s === :z ? getfield(c, :pos)[3] :
    getfield(c, s)
end
function Base.getproperty(v::G4Vector, s::Symbol)
    s === :x ? getfield(v, :vec)[1] :
    s === :y ? getfield(v, :vec)[2] :
    s === :z ? getfield(v, :vec)[3] :
    getfield(v, s)
end
function Base.getproperty(d::G4Direction, s::Symbol)
    s === :x ? getfield(d, :vec)[1] :
    s === :y ? getfield(d, :vec)[2] :
    s === :z ? getfield(d, :vec)[3] :
    getfield(d, s)
end
