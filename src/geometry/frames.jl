# frames.jl
#
# Reference-frame algebra shared by the geometry checker and the routing
# builders: the vector / rotation aliases, the row-major rotation lift, the
# rigid `Transform`, and `_world_transform` — the recursive resolution of a
# named volume's placement into world coordinates.
#
# A manifest stores geometry as plain NTuple{3} (lengths) and NTuple{9}
# (row-major rotation). These are lifted to StaticArrays / Rotations types at
# the boundary; `Transform` is the lifted, composable placement.
#
# Assumes the manifest structs (manifest.jl) are already in module scope.

# Internal vector / rotation types.
const _Vec = SVector{3,Float64}
const _Rot = RotMatrix3{Float64}

_vec(t::NTuple{3,Float64}) = _Vec(t)

# Lift a row-major NTuple{9} rotation to the true rotation matrix. The manifest
# stores rotations row-major; an SMatrix fills column-major, so the elements are
# reordered here (column k of the result is row-major row k).
_rotmatrix(r::NTuple{9,Float64}) = SMatrix{3,3,Float64,9}(
    r[1], r[4], r[7],
    r[2], r[5], r[8],
    r[3], r[6], r[9])

_rot(r::NTuple{9,Float64}) = _Rot(_rotmatrix(r))

# ---------------------------------------------------------------------------
#  Transform — a rigid placement
# ---------------------------------------------------------------------------

"""
    Transform(translation, rotation)
    Transform()                       # identity

A rigid placement of one volume in another's frame: a `translation` (mm) and a
`rotation`. Applying it to a point — `tf(p)` — gives `rotation*p + translation`;
composing — `a ∘ b` — gives the transform equivalent to applying `b` then `a`,
so a child volume's world transform is `parent_world ∘ local`.
"""
struct Transform
    translation::_Vec
    rotation::_Rot
end

Transform() = Transform(zero(_Vec), one(_Rot))

# apply the transform to a point
(tf::Transform)(p) = tf.rotation * _Vec(p) + tf.translation

# compose: (a ∘ b)(p) == a(b(p))
Base.:∘(a::Transform, b::Transform) =
    Transform(a.rotation * b.translation + a.translation, a.rotation * b.rotation)

# inverse: `inv(tf)(p)` lifts a point from `tf`'s frame back to the frame `tf`
# was applied in (i.e. `inv(tf) ∘ tf == identity`).
Base.inv(tf::Transform) = Transform(-(tf.rotation' * tf.translation), tf.rotation')

# ---------------------------------------------------------------------------
#  World-frame resolution
# ---------------------------------------------------------------------------
#
# A FiberEntry's start/stop are NOT necessarily world coordinates: they are
# expressed in the frame of the fibre's `reference` volume if one is set
# (GODDESS SetFibreReferenceVolume), otherwise in its `mother` volume. B2's
# u-turn segments have mother="world" but reference="scint_N", so each bar's
# routing is written in that bar's local frame. Two volumes must therefore be
# lifted into a common (world) frame before they are compared.

# A reusable world-frame resolver: a name -> ScintEntry index (built once, so
# resolution is not an O(placements) scan per lookup) plus a memo of resolved
# world transforms.
struct _FrameCache
    scints::Dict{String,ScintEntry}
    memo::Dict{String,Transform}
end

_FrameCache(m::GeometryManifest) = _FrameCache(
    Dict{String,ScintEntry}(p.name => p for p in m.placements if p isa ScintEntry),
    Dict{String,Transform}())

"""
    _world_transform(fc::_FrameCache, name) -> Transform

World-frame placement of a named volume. "world" is the identity; any other
name must resolve to a `ScintEntry`, whose own `mother` chain is composed
recursively. Results are memoised in `fc`.
"""
function _world_transform(fc::_FrameCache, name::AbstractString)
    name == "world" && return Transform()
    haskey(fc.memo, name) && return fc.memo[name]
    s = get(fc.scints, name, nothing)
    s === nothing && error("geometry: volume '$name' not found as a " *
        "scintillator — a fibre reference/mother frame must resolve to a " *
        "ScintEntry or 'world'")
    tf = _world_transform(fc, s.mother) ∘ Transform(_vec(s.pos), _rot(s.rot))
    fc.memo[name] = tf
    return tf
end
