# dubins.jl
#
# Dubins shortest-path router — one of the two fibre routers (see also
# `add_fiber_path!` in fiber_routing.jl).
#
# "Shortest smooth path between two directed endpoints with bounded curvature"
# is exactly a Dubins path — provably one of six words (LSL/LSR/RSL/RSR/RLR/LRL),
# closed-form. The solver itself is the `Dubins.jl` package; this file is the
# 2D<->3D embedding and the conversion of the Dubins solution into the manifest
# straight/bent fibre primitives ([`RouteSegment`](@ref)s defined in
# fiber_routing.jl).
#
# The router is planar: both endpoints must lie in a common working plane
# (normal `plane_normal`). Routing happens in that plane's 2D coordinates; arcs
# come back as `bent` segments. The turn -> bend_axis sign is chosen so the arc
# round-trips through `geometry_check.jl`'s `_fiber_polyline`: a left turn maps
# to `bend_axis = +plane_normal`, a right turn to `-plane_normal` (verified
# against the manifest arc convention — centre = midpoint + R cos(θ/2)(n̂×ĉ),
# swept right-hand about n̂).

# Dubins path word -> its three segment types, in order.
const _SEGMENT_TYPES = Dict{DubinsPathType,NTuple{3,Symbol}}(
    LSL => (:L, :S, :L), LSR => (:L, :S, :R), RSL => (:R, :S, :L),
    RSR => (:R, :S, :R), RLR => (:R, :L, :R), LRL => (:L, :R, :L))

# Right-handed in-plane basis (e1, e2, n̂) with e1 × e2 = n̂.
function _plane_basis(n::SVector{3,Float64})
    n̂ = normalize(n)
    # reference axis least aligned with n̂, so the cross product is well-conditioned
    ax, ay, az = abs(n̂[1]), abs(n̂[2]), abs(n̂[3])
    a = (ax <= ay && ax <= az) ? SVector(1.0, 0.0, 0.0) :
        (ay <= az)             ? SVector(0.0, 1.0, 0.0) :
                                 SVector(0.0, 0.0, 1.0)
    e1 = normalize(cross(a, n̂))
    e2 = cross(n̂, e1)
    return e1, e2, n̂
end

# Resolve the Dubins working plane for one (start, stop) pose pair: take a
# user-supplied `plane_normal` as-is, or derive one that contains both points
# and both tangents. The derived normal is the first non-degenerate of
#   d0 × d1                          # turn's natural normal
#   (p1 - p0) × d0                   # fallback if tangents are parallel
#   (p1 - p0) × d1                   # fallback if chord is parallel to d0
# Falling back to any axis perpendicular to d0 when all three are degenerate
# (everything collinear → a pure straight, the plane choice is irrelevant).
function _resolve_plane_normal(start::G4DirectedPoint, stop::G4DirectedPoint,
                               provided, planar_tol::Real)
    provided === nothing || return SVector{3,Float64}(provided)
    p0 = _vec(start.point.pos)
    p1 = _vec(stop.point.pos)
    d0 = start.direction.vec
    d1 = stop.direction.vec
    chord = p1 - p0
    for c in (cross(d0, d1), cross(chord, d0), cross(chord, d1))
        n = norm(c)
        n > planar_tol && return c / n
    end
    return _plane_basis(d0)[1]                       # any axis ⟂ d0
end

"""
    route_fiber(start::G4DirectedPoint, stop::G4DirectedPoint; min_radius,
                plane_normal=nothing, planar_tol=1e-6) -> Vector{RouteSegment}

Compute the shortest bounded-curvature (Dubins) fibre path from `start` to
`stop`, turning no tighter than `min_radius`. The result is an ordered list of
≤3 `RouteSegment`s (zero-length segments dropped).

The router is planar (2-D embedded in 3-D). `plane_normal` is the working
plane's normal; when left as `nothing` (the default) it is derived per call
from the two poses — the plane containing both points and both tangents. Pass
an explicit value to pin a working plane (e.g. for poses that happen to be
collinear and so admit many).

Pair with [`fiber_entries`](@ref) to turn the segments into named `FiberEntry`s
for a geometry manifest.
"""
function route_fiber(start::G4DirectedPoint, stop::G4DirectedPoint;
                     min_radius, plane_normal=nothing,
                     planar_tol::Real=1e-6)
    min_radius = _to_mm(min_radius)                  # accept unitful lengths
    min_radius > 0 || error("route_fiber: min_radius must be positive")
    start.point.ref == stop.point.ref || error(
        "route_fiber: 'start' is in frame '$(start.point.ref)' but 'stop' " *
        "is in frame '$(stop.point.ref)' — both endpoints must share a frame")
    e1, e2, n̂ = _plane_basis(_resolve_plane_normal(start, stop, plane_normal,
                                                   planar_tol))
    o  = _vec(start.point.pos)
    sp = _vec(stop.point.pos)

    # Planarity check: with an auto-derived normal this passes by construction,
    # but it still validates a user-supplied `plane_normal` that does not
    # actually contain the poses.
    dout = abs(dot(sp - o, n̂))
    dout <= planar_tol || error("route_fiber: end point is $(dout) mm out of " *
        "the working plane (normal $(Tuple(n̂))); the router is planar")
    abs(dot(start.direction.vec, n̂)) <= planar_tol ||
        error("route_fiber: start direction is not in the working plane")
    abs(dot(stop.direction.vec, n̂)) <= planar_tol ||
        error("route_fiber: end direction is not in the working plane")

    # endpoints in the plane's 2D coordinates (start sits at the 2D origin)
    θ0 = atan(dot(start.direction.vec, e2), dot(start.direction.vec, e1))
    θ1 = atan(dot(stop.direction.vec, e2), dot(stop.direction.vec, e1))
    q0 = [0.0, 0.0, θ0]
    q1 = [dot(sp - o, e1), dot(sp - o, e2), θ1]

    errcode, path = dubins_shortest_path(q0, q1, min_radius)
    errcode == EDUBOK || error("route_fiber: Dubins solver returned error " *
        "code $(errcode) — no path for these poses at min_radius=$(min_radius)")

    types  = _SEGMENT_TYPES[path.path_type]
    seglen = ntuple(i -> dubins_segment_length(path, i), 3)
    bound  = (0.0, seglen[1], seglen[1] + seglen[2],
              seglen[1] + seglen[2] + seglen[3])
    lift(c) = o + c[1] * e1 + c[2] * e2

    segs = RouteSegment[]
    for i in 1:3
        path.params[i] <= planar_tol && continue           # drop zero-length
        _, c0 = dubins_path_sample(path, bound[i])
        _, c1 = dubins_path_sample(path, bound[i + 1])
        p0, p1 = lift(c0), lift(c1)
        if types[i] === :S
            push!(segs, RouteSegment(STRAIGHT, p0, p1, 0.0,
                                     zero(SVector{3,Float64})))
        else
            axis = types[i] === :L ? n̂ : -n̂
            push!(segs, RouteSegment(BENT, p0, p1, path.params[i], axis))
        end
    end
    return segs
end
