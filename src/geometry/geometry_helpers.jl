# geometry_helpers.jl
#
# Reusable geometric primitives — oriented bounding boxes, corner fillets, and
# the bent-fibre arc reconstruction — shared by the geometry checker and the
# routing builders. Pure geometry on StaticArrays / Rotations: the only
# manifest coupling is the `OBB(::ScintEntry, ::Transform)` convenience
# constructor.

"""
    OBB(center, axes, half)

An oriented bounding box. `center` is the world-frame centre; `axes` is a 3×3
matrix whose *columns* are the box's unit axes in world coordinates; `half`
holds the half-extents along those axes.
"""
struct OBB
    center::SVector{3,Float64}
    axes::SMatrix{3,3,Float64,9}
    half::SVector{3,Float64}
end

"""
    OBB(scint::ScintEntry, tf::Transform)

The oriented bounding box of a scintillator tile placed by world transform
`tf`: `tf`'s translation is the box centre, its rotation columns the box axes,
and half the tile's `dims` the extents.
"""
OBB(scint::ScintEntry, tf::Transform) =
    OBB(tf.translation, SMatrix{3,3,Float64,9}(tf.rotation),
        SVector{3,Float64}(scint.dims) ./ 2)

# A world point expressed in the box's local (axis-aligned) coordinates.
_obb_local(box::OBB, p::SVector{3,Float64}) = box.axes' * (p - box.center)

"""
    point_in_obb(p, box; tol=0.0) -> Bool

True if world point `p` lies inside oriented box `box` (within `tol`).
"""
function point_in_obb(p::SVector{3,Float64}, box::OBB; tol::Real=0.0)
    q = _obb_local(box, p)
    return abs(q[1]) <= box.half[1] + tol &&
           abs(q[2]) <= box.half[2] + tol &&
           abs(q[3]) <= box.half[3] + tol
end

"""
    segment_obb_interval(p0, p1, box; tol=1e-9)
        -> Union{Nothing,Tuple{Float64,Float64}}

The sub-interval `(t_in, t_out) ⊆ [0,1]` of the segment `p0 + t·(p1-p0)` that
lies inside `box`, or `nothing` if the segment never enters it (or only grazes
it at a point). Slab method, evaluated in the box's local frame.
"""
function segment_obb_interval(p0::SVector{3,Float64}, p1::SVector{3,Float64},
                              box::OBB; tol::Real=1e-9)
    a = _obb_local(box, p0)
    d = _obb_local(box, p1) - a
    tmin, tmax = 0.0, 1.0
    for k in 1:3
        h = box.half[k]
        if abs(d[k]) < tol
            # segment runs parallel to this slab — must already lie within it
            (-h - tol <= a[k] <= h + tol) || return nothing
        else
            lo, hi = minmax((-h - a[k]) / d[k], (h - a[k]) / d[k])
            tmin = max(tmin, lo)
            tmax = min(tmax, hi)
            tmin > tmax && return nothing
        end
    end
    tmax - tmin <= tol && return nothing          # grazes at a point only
    return (tmin, tmax)
end

"""
    obb_overlap(A::OBB, B::OBB) -> Bool

True iff two oriented boxes intersect. Separating-axis theorem (Ericson,
Real-Time Collision Detection §4.4.1): 3 face axes of `A`, 3 of `B`, and the
9 edge-cross axes `A_i × B_j`.
"""
function obb_overlap(A::OBB, B::OBB)
    # each box's local axes expressed in world coordinates (rotation columns)
    uA = (A.axes[:, 1], A.axes[:, 2], A.axes[:, 3])
    uB = (B.axes[:, 1], B.axes[:, 2], B.axes[:, 3])
    eA, eB = A.half, B.half
    # R[i][j] = uA[i]·uB[j]; AbsR adds an epsilon to stay robust when two edges
    # are parallel (their cross product is ~0 and that axis is degenerate)
    R    = ntuple(i -> ntuple(j -> dot(uA[i], uB[j]), 3), 3)
    AbsR = ntuple(i -> ntuple(j -> abs(R[i][j]) + 1e-9, 3), 3)
    # translation from A to B, in A's frame
    d = B.center - A.center
    t = (dot(d, uA[1]), dot(d, uA[2]), dot(d, uA[3]))

    # 3 face axes of A
    for i in 1:3
        ra = eA[i]
        rb = eB[1]*AbsR[i][1] + eB[2]*AbsR[i][2] + eB[3]*AbsR[i][3]
        abs(t[i]) > ra + rb && return false
    end
    # 3 face axes of B
    for j in 1:3
        ra = eA[1]*AbsR[1][j] + eA[2]*AbsR[2][j] + eA[3]*AbsR[3][j]
        rb = eB[j]
        abs(t[1]*R[1][j] + t[2]*R[2][j] + t[3]*R[3][j]) > ra + rb && return false
    end
    # 9 edge-cross axes A_i × B_j
    abs(t[3]*R[2][1]-t[2]*R[3][1]) > eA[2]*AbsR[3][1]+eA[3]*AbsR[2][1]+eB[2]*AbsR[1][3]+eB[3]*AbsR[1][2] && return false
    abs(t[3]*R[2][2]-t[2]*R[3][2]) > eA[2]*AbsR[3][2]+eA[3]*AbsR[2][2]+eB[1]*AbsR[1][3]+eB[3]*AbsR[1][1] && return false
    abs(t[3]*R[2][3]-t[2]*R[3][3]) > eA[2]*AbsR[3][3]+eA[3]*AbsR[2][3]+eB[1]*AbsR[1][2]+eB[2]*AbsR[1][1] && return false
    abs(t[1]*R[3][1]-t[3]*R[1][1]) > eA[1]*AbsR[3][1]+eA[3]*AbsR[1][1]+eB[2]*AbsR[2][3]+eB[3]*AbsR[2][2] && return false
    abs(t[1]*R[3][2]-t[3]*R[1][2]) > eA[1]*AbsR[3][2]+eA[3]*AbsR[1][2]+eB[1]*AbsR[2][3]+eB[3]*AbsR[2][1] && return false
    abs(t[1]*R[3][3]-t[3]*R[1][3]) > eA[1]*AbsR[3][3]+eA[3]*AbsR[1][3]+eB[1]*AbsR[2][2]+eB[2]*AbsR[2][1] && return false
    abs(t[2]*R[1][1]-t[1]*R[2][1]) > eA[1]*AbsR[2][1]+eA[2]*AbsR[1][1]+eB[2]*AbsR[3][3]+eB[3]*AbsR[3][2] && return false
    abs(t[2]*R[1][2]-t[1]*R[2][2]) > eA[1]*AbsR[2][2]+eA[2]*AbsR[1][2]+eB[1]*AbsR[3][3]+eB[3]*AbsR[3][1] && return false
    abs(t[2]*R[1][3]-t[1]*R[2][3]) > eA[1]*AbsR[2][3]+eA[2]*AbsR[1][3]+eB[1]*AbsR[3][2]+eB[2]*AbsR[3][1] && return false
    return true
end

"""
    arc_polyline(p0, p1, bend_angle, bend_axis; arc_segments=24)
        -> Vector{SVector{3,Float64}}

Sample a manifest bent-fibre arc into a polyline of `arc_segments` chords. The
arc convention (shared by the geometry checker and both fibre routers):

    R      = |p1 - p0| / (2 sin(θ/2))
    centre = midpoint(p0,p1) + R cos(θ/2) (n̂ × ĉ)

where θ = `bend_angle`, n̂ = unit(`bend_axis`), ĉ = unit(p1 - p0); the arc is
swept from `p0` to `p1` by a right-hand rotation about n̂. A degenerate arc
(zero chord or zero angle) collapses to the two endpoints.

Chord sampling sits slightly *inside* the true arc, so overlap detection on the
polyline is marginally optimistic — raise `arc_segments` to tighten it.
"""
function arc_polyline(p0::SVector{3,Float64}, p1::SVector{3,Float64},
                      bend_angle::Real, bend_axis; arc_segments::Int=24)
    θ = bend_angle
    chord = p1 - p0
    L = norm(chord)
    (L == 0 || θ == 0) && return SVector{3,Float64}[p0, p1]
    R      = L / (2 * sin(θ / 2))
    ĉ      = chord / L
    n̂      = normalize(SVector{3,Float64}(bend_axis))
    center = (p0 + p1) / 2 + R * cos(θ / 2) * cross(n̂, ĉ)
    v0     = p0 - center
    return SVector{3,Float64}[center + AngleAxis(θ * i / arc_segments, n̂...) * v0
                              for i in 0:arc_segments]
end

"""
    fillet_corner(p_prev, p_corner, p_next, radius; collinear_tol=1e-9)
        -> Union{Nothing,NamedTuple}

The circular fillet of the corner `p_prev → p_corner → p_next`: an arc of
`radius` tangent to both edges. Returns `nothing` when the corner is collinear
(no fillet needed), otherwise a named tuple

    (tangent_in, tangent_out, angle, axis, trim)

where `tangent_in`/`tangent_out` are the arc endpoints on the two edges,
`angle` is the deflection (radian), `axis` the right-hand bend axis, and `trim`
the tangent length consumed along each edge. Throws if the corner is too sharp
to fillet (a near-reversal).
"""
function fillet_corner(p_prev::SVector{3,Float64}, p_corner::SVector{3,Float64},
                       p_next::SVector{3,Float64}, radius::Real;
                       collinear_tol::Real=1e-9)
    ein  = p_corner - p_prev
    eout = p_next - p_corner
    Lin, Lout = norm(ein), norm(eout)
    (Lin > 0 && Lout > 0) ||
        error("fillet_corner: zero-length edge at corner $(Tuple(p_corner))")
    u_in  = ein  / Lin
    u_out = eout / Lout
    θ = acos(clamp(dot(u_in, u_out), -1.0, 1.0))      # deflection angle
    θ <= collinear_tol && return nothing               # collinear: no corner
    θ >= π - 1e-6 && error("fillet_corner: corner at $(Tuple(p_corner)) is a " *
        "near-reversal ($(round(rad2deg(θ), digits=1))°) and cannot be " *
        "filleted — add an intermediate waypoint")
    axisv = cross(u_in, u_out)
    na = norm(axisv)
    na <= collinear_tol && return nothing
    axis = axisv / na
    trim = radius * tan(θ / 2)
    return (tangent_in  = p_corner - trim * u_in,
            tangent_out = p_corner + trim * u_out,
            angle = θ, axis = axis, trim = trim)
end

# ---------------------------------------------------------------------------
#  Scintillator-tile extent
# ---------------------------------------------------------------------------

const _AXIS_INDEX = (x = 1, y = 2, z = 3)

function _axis_index(axis::Symbol)
    hasproperty(_AXIS_INDEX, axis) ||
        error("extent: axis must be :x, :y or :z, got :$axis")
    return getproperty(_AXIS_INDEX, axis)
end

# (origin, unit direction) — in the tile's mother frame — of the `k`-th axis of
# the frame named by `ref`, against which `extent` projects the tile.
function _ref_axis(s::ScintEntry, k::Int, ref)
    if ref === :mother || ref === :world
        ref === :world && s.mother != "world" && error(
            "extent: scintillator '$(s.name)' is placed in '$(s.mother)', " *
            "not world — a world-frame extent needs the full geometry; pass " *
            "the mother tile as `ref` instead, or query in :mother")
        return zero(SVector{3,Float64}),
               SVector(ntuple(i -> i == k ? 1.0 : 0.0, 3)...)
    elseif ref isa ScintEntry
        ref.mother == s.mother || error(
            "extent: reference tile '$(ref.name)' (placed in '$(ref.mother)') " *
            "does not share the frame of '$(s.name)' (placed in '$(s.mother)')")
        return SVector{3,Float64}(ref.pos), _rotmatrix(ref.rot)[:, k]
    else
        error("extent: ref must be :world, :mother or a ScintEntry, got $ref")
    end
end

"""
    extent(s::ScintEntry, axis::Symbol, ref=:mother) -> (; min, max)

The span of scintillator tile `s` along one coordinate axis, in millimetres.
`axis` is `:x`, `:y` or `:z`. `ref` selects the frame the axis and the span are
measured in:

  * `:mother` (default) — the tile's own mother frame (the frame its `pos`/`rot`
    are stored in);
  * `:world` — the world frame (requires `s` to be world-placed, so the mother
    frame *is* the world frame);
  * another `ScintEntry` — that tile's local frame (it must share `s`'s mother
    frame, i.e. both placed in the same volume).

For an axis-aligned tile this is just the two face positions; for a rotated tile
it is the projected span (the support of the oriented box along the axis).
"""
function extent(s::ScintEntry, axis::Symbol, ref=:mother)
    k = _axis_index(axis)
    origin, d = _ref_axis(s, k, ref)
    As = _rotmatrix(s.rot)
    h  = SVector{3,Float64}(s.dims) ./ 2
    center = dot(d, SVector{3,Float64}(s.pos) - origin)
    radius = h[1] * abs(dot(As[:, 1], d)) +
             h[2] * abs(dot(As[:, 2], d)) +
             h[3] * abs(dot(As[:, 3], d))
    return (; min = center - radius, max = center + radius)
end

"""
    center_pos(s::ScintEntry) -> G4Coordinate

The tile's centre as a [`G4Coordinate`](@ref) in its mother frame — bundles
`s.pos` with `s.mother`. Pairs with [`face_center`](@ref).
"""
center_pos(s::ScintEntry) = G4Coordinate(s.pos, s.mother)

"""
    face_center(s::ScintEntry, axis::Symbol, side::Symbol) -> G4Coordinate

The centre of one face of scintillator tile `s`, as a [`G4Coordinate`](@ref) in
the tile's mother frame. `axis` is `:x`, `:y` or `:z` and `side` is `:max` (the
+`axis` face of the tile's *own* frame) or `:min` (the −`axis` face). For an
axis-aligned tile this is the face at the larger / smaller coordinate; for a
rotated tile it is still the centre of that local face.
"""
function face_center(s::ScintEntry, axis::Symbol, side::Symbol)
    k = _axis_index(axis)
    sgn = side === :max ?  1.0 :
          side === :min ? -1.0 :
          error("face_center: side must be :min or :max, got :$side")
    offset = SVector(ntuple(i -> i == k ? sgn * s.dims[k] / 2 : 0.0, 3)...)
    point  = SVector{3,Float64}(s.pos) + _rotmatrix(s.rot) * offset
    return G4Coordinate(Tuple(point), s.mother)
end
