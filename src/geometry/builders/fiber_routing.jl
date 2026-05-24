# fiber_routing.jl
#
# Fibre routing — the shared vocabulary plus the waypoint-path routing logic.
#
# Two routers turn a route description into manifest fibre primitives:
#   * `route_fiber`  — shortest bounded-curvature path between two directed
#                      endpoints (Dubins). Lives in dubins.jl.
#   * `add_fiber_path!` — a single continuous fibre threaded through a list of
#                      waypoints; lives in manifest_builder.jl with the other
#                      `add_*!` builders. The waypoint-fillet routing logic it
#                      calls lives here.
#
# Both speak the two types defined here: `G4DirectedPoint` (a point plus a
# direction) and `RouteSegment` (one straight or circular-arc piece of a route).
# `fiber_entries` lowers a `RouteSegment` list into manifest `FiberEntry`s.

"""
    G4DirectedPoint(point::G4Coordinate, direction)
    G4DirectedPoint(point, direction)            # `point` taken in the world frame
    G4DirectedPoint(point, ref, direction)       # `point` in frame `ref`

A point with a direction: a [`G4Coordinate`](@ref) `point` and a
[`G4Direction`](@ref) `direction`, both in the same reference frame (enforced
at construction). Unlike a full pose it pins only one tangent, not a complete
orientation. For [`route_fiber`](@ref) the two endpoints must share a frame,
lie in the working plane, and both directions must be in that plane.

`direction` may be a `G4Direction` or any bare 3-element value (tuple, vector,
`SVector`); a bare value is coerced to a `G4Direction` in `point`'s frame.
"""
struct G4DirectedPoint
    point::G4Coordinate
    direction::G4Direction
    function G4DirectedPoint(point::G4Coordinate, direction::G4Direction)
        point.ref == direction.ref || error(
            "G4DirectedPoint: point and direction must share a reference frame " *
            "(point.ref='$(point.ref)', direction.ref='$(direction.ref)')")
        return new(point, direction)
    end
end

# convenience: a bare direction (tuple/vector/SVector) inherits the point's
# frame — the common case in detector specs
G4DirectedPoint(point::G4Coordinate, direction) =
    G4DirectedPoint(point, G4Direction(direction, point.ref))

# convenience: a bare `point` (tuple/vector) is taken in the world frame
G4DirectedPoint(point, direction) =
    G4DirectedPoint(G4Coordinate(point, "world"), direction)
# convenience: an explicit (point, ref, direction) triple
G4DirectedPoint(point, ref::Union{Symbol,AbstractString}, direction) =
    G4DirectedPoint(G4Coordinate(point, ref), direction)

# Translate a `G4DirectedPoint` by a displacement; the direction is a tangent
# and is left unchanged. Adding/subtracting a `G4Vector` requires the two to
# share a reference frame (delegated to `G4Coordinate`'s arithmetic); a bare
# tuple/vector is treated as a relative offset in the point's own frame.
Base.:+(p::G4DirectedPoint, v::G4Vector) = G4DirectedPoint(p.point + v, p.direction)
Base.:+(v::G4Vector, p::G4DirectedPoint) = p + v
Base.:-(p::G4DirectedPoint, v::G4Vector) = G4DirectedPoint(p.point - v, p.direction)

Base.:+(p::G4DirectedPoint, v::Union{Tuple,AbstractVector}) =
    G4DirectedPoint(p.point + v, p.direction)
Base.:+(v::Union{Tuple,AbstractVector}, p::G4DirectedPoint) = p + v
Base.:-(p::G4DirectedPoint, v::Union{Tuple,AbstractVector}) =
    G4DirectedPoint(p.point - v, p.direction)

"""
    Waypoint

Alias for `Union{G4Coordinate, G4DirectedPoint}` — the element type accepted
in [`add_fiber_path!`](@ref)'s `waypoints` list. Use `Waypoint[…]` to
construct a vector that can hold either form so successive `push!`es with
mixed types don't trip the element-type inference of a bare `[…]` literal.
"""
const Waypoint = Union{G4Coordinate, G4DirectedPoint}

"""
    FiberKind

The two fibre-segment shapes: `STRAIGHT` and `BENT`. A [`RouteSegment`](@ref)
carries one; the manifest [`FiberEntry`](@ref) stores the same distinction as
its `kind` string (`"straight"` / `"bent"`) at the C++ serialization boundary.
"""
@enum FiberKind STRAIGHT BENT

# manifest `kind` string (the C++ text format) for a fibre-segment shape
_kindstr(k::FiberKind) = k === BENT ? "bent" : "straight"

"""
One segment of a routed fibre path: `kind` is a [`FiberKind`](@ref) (`STRAIGHT`
or `BENT`); `start` and `stop` are 3D endpoints; `bend_angle` (rad) and
`bend_axis` (unit vector) are meaningful for `BENT` segments only.
"""
struct RouteSegment
    kind::FiberKind
    start::SVector{3,Float64}
    stop::SVector{3,Float64}
    bend_angle::Float64
    bend_axis::SVector{3,Float64}
end

"""
    fiber_entries(segs; name, material_file, mother="world", reference="",
                  glued=false, glue_file="", glue_profile="", loop_id=-1)
        -> Vector{FiberEntry}

Convert routed `RouteSegment`s into manifest `FiberEntry`s. Segments are named
`"\$(name)_1"`, `"\$(name)_2"`, … and share the metadata keyword arguments; pass
a common `loop_id` so [`check_geometry`](@ref) treats them as one physical
fibre (its segments chain end-to-end and must not be flagged against each
other).
"""
function fiber_entries(segs::Vector{RouteSegment}; name::AbstractString,
                       material_file::AbstractString,
                       mother::AbstractString="world",
                       reference::AbstractString="", glued::Bool=false,
                       glue_file::AbstractString="",
                       glue_profile::AbstractString="",
                       start_reflectivity::Real=NaN,
                       end_reflectivity::Real=NaN,
                       loop_id::Int=-1)
    # The fibre's physical endpoints are the first segment's `start` and the
    # last segment's `stop`; internal joins are not physical fibre ends, so
    # reflectivities are applied only at the boundaries.
    n = length(segs)
    return FiberEntry[
        FiberEntry(
            name = "$(name)_$(i)",
            kind = _kindstr(s.kind),
            mother = mother,
            start = Tuple(s.start), stop = Tuple(s.stop),
            bend_angle = s.bend_angle, bend_axis = Tuple(s.bend_axis),
            material_file = material_file, reference = reference,
            glued = glued, glue_file = glue_file, glue_profile = glue_profile,
            start_reflectivity = i == 1 ? Float64(start_reflectivity) : NaN,
            end_reflectivity   = i == n ? Float64(end_reflectivity)   : NaN,
            loop_id = loop_id)
        for (i, s) in enumerate(segs)]
end

# ---------------------------------------------------------------------------
#  Waypoint path router
# ---------------------------------------------------------------------------
#
# `add_fiber_path!` (manifest_builder.jl) threads one continuous fibre through a
# waypoint list. The geometry is here: a gap between two pinned waypoints
# (`G4DirectedPoint`) is handed to the Dubins `route_fiber`; a run of free gaps
# becomes straight segments with a `min_radius` fillet at each corner. The
# routed straights are then split at scintillator boundaries and each piece is
# assigned its mother volume by containment.

const _ZERO3 = zero(SVector{3,Float64})

# world position of a waypoint already resolved into the world frame
_wp_pos(w::G4DirectedPoint) = _vec(w.point.pos)
_wp_pos(w::G4Coordinate)    = _vec(w.pos)

# fail if two unit directions disagree by more than `tol`
_check_dir(a, b, msg::AbstractString, tol::Real) =
    norm(a - b) <= tol || error(msg)

"""
    route_path(waypoints; min_radius, plane_normal=nothing, continuity_tol=1e-3)
        -> Vector{RouteSegment}

Route one continuous path through `waypoints` — each a [`G4Coordinate`](@ref)
(free point) or [`G4DirectedPoint`](@ref) (pinned tangent), **already resolved
into the world frame**. A gap between two pinned waypoints is a Dubins curve
([`route_fiber`](@ref)); a maximal run of straight gaps is straight segments
filleted with `min_radius` arcs at each free corner. Errors on a tangent kink
at a pinned waypoint. `plane_normal` is propagated to each Dubins section —
`nothing` (default) lets every section derive its own working plane from its
two poses. See [`add_fiber_path!`](@ref) for the builder front end.
"""
function route_path(waypoints::AbstractVector; min_radius,
                    plane_normal=nothing, continuity_tol::Real=1e-3)
    min_radius = _to_mm(min_radius)                  # accept unitful lengths
    n = length(waypoints)
    n >= 2 || error("route_path: need at least 2 waypoints, got $n")
    pinned(w) = w isa G4DirectedPoint

    segs = RouteSegment[]
    i = 1
    while i <= n - 1
        if pinned(waypoints[i]) && pinned(waypoints[i + 1])
            append!(segs, route_fiber(waypoints[i], waypoints[i + 1];
                          min_radius = min_radius, plane_normal = plane_normal))
            i += 1
        else
            j = i                                  # extend a maximal straight run
            while j <= n - 1 && !(pinned(waypoints[j]) && pinned(waypoints[j + 1]))
                j += 1
            end
            append!(segs, _fillet_run(view(waypoints, i:j),
                                      min_radius, continuity_tol))
            i = j
        end
    end
    return segs
end

# Route a run of waypoints joined by straight gaps: straight segments with a
# fillet arc at each free interior corner. Pinned interior/endpoint waypoints
# are not filleted — the path must pass through them without a kink.
function _fillet_run(run, r::Real, tol::Real)
    m = length(run)
    pts = SVector{3,Float64}[_wp_pos(w) for w in run]

    if run[1] isa G4DirectedPoint
        _check_dir(normalize(pts[2] - pts[1]), run[1].direction.vec,
            "add_fiber_path!: fibre leaves a pinned waypoint at a kink", tol)
    end
    if run[m] isa G4DirectedPoint
        _check_dir(normalize(pts[m] - pts[m - 1]), run[m].direction.vec,
            "add_fiber_path!: fibre arrives at a pinned waypoint at a kink", tol)
    end

    m == 2 && return RouteSegment[RouteSegment(STRAIGHT, pts[1], pts[2], 0.0, _ZERO3)]

    segs = RouteSegment[]
    cursor = pts[1]
    for k in 2:(m - 1)
        if run[k] isa G4DirectedPoint
            din  = normalize(pts[k] - pts[k - 1])
            dout = normalize(pts[k + 1] - pts[k])
            _check_dir(din, dout,
                "add_fiber_path!: kink at a pinned interior waypoint", tol)
            _check_dir(din, run[k].direction.vec,
                "add_fiber_path!: a pinned waypoint's direction disagrees " *
                "with the straight run through it", tol)
        else
            f = fillet_corner(pts[k - 1], pts[k], pts[k + 1], r)
            f === nothing && continue              # collinear free point
            edir = normalize(pts[k] - pts[k - 1])
            v = f.tangent_in - cursor
            dot(v, edir) < -tol && error("add_fiber_path!: corner fillets " *
                "overlap — min_radius too large for the waypoint spacing")
            norm(v) > tol &&
                push!(segs, RouteSegment(STRAIGHT, cursor, f.tangent_in, 0.0, _ZERO3))
            push!(segs, RouteSegment(BENT, f.tangent_in, f.tangent_out,
                                     f.angle, f.axis))
            cursor = f.tangent_out
        end
    end
    norm(pts[m] - cursor) > tol &&
        push!(segs, RouteSegment(STRAIGHT, cursor, pts[m], 0.0, _ZERO3))
    return segs
end

# Reconstruct a RouteSegment's centreline as a polyline (shared `arc_polyline`).
function _route_polyline(s::RouteSegment; arc_segments::Int=16)
    s.kind === STRAIGHT && return SVector{3,Float64}[s.start, s.stop]
    return arc_polyline(s.start, s.stop, s.bend_angle, s.bend_axis;
                        arc_segments = arc_segments)
end

# Resolve a list of (possibly mother-relative) waypoints into the world frame.
function _resolve_path_waypoints(placements, waypoints)
    fc = _FrameCache(GeometryManifest(placements = placements))
    out = Union{G4Coordinate,G4DirectedPoint}[]
    for w in waypoints
        if w isa G4DirectedPoint
            tf = _world_transform(fc, w.point.ref)
            push!(out, G4DirectedPoint(
                G4Coordinate(Tuple(tf(w.point.pos)), "world"),
                tf.rotation * w.direction.vec))
        elseif w isa G4Coordinate
            tf = _world_transform(fc, w.ref)
            push!(out, G4Coordinate(Tuple(tf(w.pos)), "world"))
        else
            error("add_fiber_path!: a waypoint must be a G4Coordinate or a " *
                  "G4DirectedPoint, got $(typeof(w))")
        end
    end
    return out
end

# World-frame oriented bounding box of every scintillator tile.
function _scint_obbs(placements)
    fc = _FrameCache(GeometryManifest(placements = placements))
    obbs = Tuple{String,OBB}[]
    for p in placements
        p isa ScintEntry || continue
        tf = _world_transform(fc, p.name)
        push!(obbs, (p.name, OBB(tf.translation,
                                 SMatrix{3,3,Float64}(tf.rotation),
                                 SVector{3,Float64}(p.dims) ./ 2)))
    end
    return obbs
end

# Name of the first scintillator OBB containing `p`, else "world".
function _mother_of(p::SVector{3,Float64}, obbs)
    for (nm, box) in obbs
        point_in_obb(p, box) && return nm
    end
    return "world"
end

# Assign each routed segment its mother volume by midpoint containment. A
# straight is *not* split at scintillator boundaries even if it crosses one:
# every G4Fibre carries lossy roughened endcaps at its physical ends, so
# splitting a fibre at the scint wall puts those endcaps inside the wrap
# shell and effectively kills light propagation across the wall. Instead,
# the straight stays whole and is placed in the scintillator it passes
# through — GODDESS's G4Fibre then carves the cylinder against that scint
# (`GetOutermostVolumeOutsideMother_physicalVolume`), matching B1/B2's
# "inner fibre extends past the scint wall" pattern. The wrapping cut
# auto-derive in DetectorConstruction.cc picks up such fibres via
# `mother == w.scint` and gets a clean through-hole. Assumes a single
# straight passes through at most one scintillator — true for every design
# to date; revisit if a route ever pierces two.
#
# Bent segments are likewise not split; they are assigned by midpoint and
# validated to stay within one volume (a fillet straddling a scint wall
# would have the same endcap-loss pathology and is rejected up front).
function _split_path_mothers(segs::Vector{RouteSegment}, obbs)
    pieces = Tuple{RouteSegment,String}[]
    for s in segs
        if s.kind === STRAIGHT
            mother = _mother_of((s.start + s.stop) / 2, obbs)
            push!(pieces, (s, mother))
        else
            pl = _route_polyline(s; arc_segments = 16)
            mother = _mother_of(pl[length(pl) ÷ 2 + 1], obbs)
            for q in pl
                _mother_of(q, obbs) == mother || error(
                    "add_fiber_path!: a bent segment crosses a scintillator " *
                    "boundary — add a waypoint at the crossing so the fillet " *
                    "arc stays within one volume")
            end
            push!(pieces, (s, mother))
        end
    end
    return pieces
end
