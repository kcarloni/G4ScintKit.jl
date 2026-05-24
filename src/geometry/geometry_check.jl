# geometry_check.jl
#
# Pre-flight geometry validation: detect physics-affecting volume overlaps in a
# geometry manifest *before* any Geant4 job is launched.
#
# Why this exists. Geant4 navigation requires that daughter volumes within a
# mother do not overlap; an overlap makes the navigator ambiguous and silently
# produces wrong physics (wrong material at boundaries, mis-attributed optical
# photons, energy in the wrong volume — no exception is thrown). Two cases
# matter here:
#   * fibre-fibre  — wrong light propagation. For fibres GODDESS hard-disables
#     Geant4's SearchOverlaps and its construction-time intersection test only
#     carves a fibre against its *aunt* volumes, never against other fibres, so
#     nothing in the C++/GODDESS path catches two fibres routed through the same
#     region. (Geant4's own CheckOverlaps is also unusable here — every fibre is
#     a boolean solid, which it cannot sample.)
#   * scintillator-scintillator — overlapping active volumes double-count or
#     mis-attribute energy deposits.
# This module is the Julia-side gate for both.
#
# Bent fibres are reconstructed with the shared `arc_polyline` (geometry_helpers.jl)
# — the manifest's bent-fibre arc convention, verified against the B2 reference
# manifest (both the 90° front bends and the 180° loops).

# The vector / rotation aliases (`_Vec`, `_Rot`), the row-major rotation lift
# (`_vec`, `_rot`) and the world-frame resolver (`_FrameCache`,
# `_world_transform`) live in frames.jl, included earlier.

# ---------------------------------------------------------------------------
#  Fibre centreline reconstruction
# ---------------------------------------------------------------------------

"""
    _fiber_polyline(f::FiberEntry; arc_segments=24) -> Vector{SVector{3,Float64}}

Reconstruct a fibre segment's centreline as a polyline: the two endpoints for a
straight segment, the sampled circular arc for a bent one (see
[`arc_polyline`](@ref) for the arc convention).
"""
function _fiber_polyline(f::FiberEntry; arc_segments::Int=24)
    p0, p1 = _vec(f.start), _vec(f.stop)
    f.kind == "bent" || return _Vec[p0, p1]
    return arc_polyline(p0, p1, f.bend_angle, _vec(f.bend_axis);
                        arc_segments = arc_segments)
end

# ---------------------------------------------------------------------------
#  Segment / polyline distance
# ---------------------------------------------------------------------------

# Shortest distance between segment p1->p2 and segment q1->q2 (Ericson,
# Real-Time Collision Detection — closest points of two segments).
function _segment_distance(p1::_Vec, p2::_Vec, q1::_Vec, q2::_Vec)
    d1 = p2 - p1                    # direction of segment P
    d2 = q2 - q1                    # direction of segment Q
    r  = p1 - q1
    a  = dot(d1, d1)
    e  = dot(d2, d2)
    f  = dot(d2, r)
    local s::Float64, t::Float64
    if a <= eps() && e <= eps()
        return norm(r)              # both segments are points
    end
    if a <= eps()
        s = 0.0
        t = clamp(f / e, 0.0, 1.0)
    else
        c = dot(d1, r)
        if e <= eps()
            t = 0.0
            s = clamp(-c / a, 0.0, 1.0)
        else
            b = dot(d1, d2)
            denom = a * e - b * b
            s = denom != 0 ? clamp((b * f - c * e) / denom, 0.0, 1.0) : 0.0
            t = (b * s + f) / e
            if t < 0.0
                t = 0.0
                s = clamp(-c / a, 0.0, 1.0)
            elseif t > 1.0
                t = 1.0
                s = clamp((b - c) / a, 0.0, 1.0)
            end
        end
    end
    return norm((p1 + s * d1) - (q1 + t * d2))
end

# Minimum distance between two polylines.
function _polyline_distance(a::Vector{_Vec}, b::Vector{_Vec})
    best = Inf
    for i in 1:(length(a) - 1), j in 1:(length(b) - 1)
        d = _segment_distance(a[i], a[i + 1], b[j], b[j + 1])
        d < best && (best = d)
    end
    return best
end

# ---------------------------------------------------------------------------
#  Public API
# ---------------------------------------------------------------------------

"""A detected fibre-fibre overlap: two fibre segments in the same mother volume
whose centrelines pass closer than the clearance threshold."""
struct FiberClash
    a::String           # name of the first fibre segment
    b::String           # name of the second fibre segment
    mother::String      # the shared mother volume
    distance::Float64   # closest centreline approach (mm)
end

# Logical-fibre identity. Segments sharing a loop_id >= 0 are pieces of one
# continuous physical fibre (chained end-to-end — they *must* touch); they are
# never compared against each other. loop_id < 0 means a stand-alone fibre, so
# each such entry is its own logical fibre.
_logical_fiber(f::FiberEntry, idx::Int) =
    f.loop_id >= 0 ? (true, f.loop_id) : (false, idx)

"""
    fiber_clashes(m::GeometryManifest; clearance=1.0, arc_segments=24) -> Vector{FiberClash}

Find every pair of fibre segments that overlap. Two segments overlap when they
share a mother volume, belong to different logical fibres, and their
centrelines — lifted into world coordinates — pass closer than `clearance` (mm).
The default 1.0 mm is the nominal fibre diameter (centre-to-centre < diameter ⇒
the fibres intersect).

Only segments in the *same* mother are compared: Geant4 navigation only forbids
overlap among sibling daughters, and fibres in different scintillator bars
cannot reach each other. In practice this is the `world`-mother routing
segments — exactly where B2's u-turn fibres are placed. Each fibre's start/stop
are first resolved out of its `reference`/`mother` frame into world coordinates
(see `_world_transform`), so different bars are correctly offset.

Limitation: segments of one logical fibre (shared `loop_id`) are never compared,
so a fibre that self-intersects is not flagged — only distinct fibres are.
"""
function fiber_clashes(m::GeometryManifest; clearance=1.0, arc_segments::Int=24)
    clearance = _to_mm(clearance)                            # accept Real (mm) or a unitful length
    fc = _FrameCache(m)

    # collect each fibre with its world-coordinate centreline polyline
    entries = Tuple{FiberEntry,Int,Vector{_Vec}}[]
    for (idx, p) in enumerate(m.placements)
        p isa FiberEntry || continue
        frame = isempty(p.reference) ? p.mother : p.reference
        tf = _world_transform(fc, frame)
        world_pl = _Vec[tf(q)
                        for q in _fiber_polyline(p; arc_segments=arc_segments)]
        push!(entries, (p, idx, world_pl))
    end

    # group by Geant4 mother volume
    by_mother = Dict{String,Vector{Int}}()
    for (k, e) in enumerate(entries)
        push!(get!(by_mother, e[1].mother, Int[]), k)
    end

    clashes = FiberClash[]
    for (mother, ks) in by_mother
        for ii in eachindex(ks), jj in (ii + 1):lastindex(ks)
            fi, idxi, pli = entries[ks[ii]]
            fj, idxj, plj = entries[ks[jj]]
            _logical_fiber(fi, idxi) == _logical_fiber(fj, idxj) && continue
            d = _polyline_distance(pli, plj)
            if d < clearance
                push!(clashes, FiberClash(fi.name, fj.name, mother, d))
            end
        end
    end
    sort!(clashes, by = c -> c.distance)
    return clashes
end

"""A detected scintillator-scintillator overlap: two scintillator tiles in the
same mother volume whose (oriented) boxes intersect."""
struct ScintOverlap
    a::String           # name of the first scintillator tile
    b::String           # name of the second scintillator tile
    mother::String      # the shared mother volume
end

"""
    scint_overlaps(m::GeometryManifest) -> Vector{ScintOverlap}

Find every pair of scintillator tiles whose volumes intersect. Two tiles overlap
when they share a mother volume and their oriented boxes (resolved into world
coordinates) intersect — an overlap in the *active* detector volume, which makes
Geant4 navigation ambiguous and double-counts or mis-attributes energy deposits.

Like [`fiber_clashes`](@ref), only tiles in the same mother are compared.
"""
function scint_overlaps(m::GeometryManifest)
    fc = _FrameCache(m)

    # collect each tile with its world-frame oriented box
    scints = Tuple{ScintEntry,OBB}[]
    for p in m.placements
        p isa ScintEntry || continue
        push!(scints, (p, OBB(p, _world_transform(fc, p.name))))
    end

    # group by mother volume
    by_mother = Dict{String,Vector{Int}}()
    for (k, s) in enumerate(scints)
        push!(get!(by_mother, s[1].mother, Int[]), k)
    end

    overlaps = ScintOverlap[]
    for (mother, ks) in by_mother
        for ii in eachindex(ks), jj in (ii + 1):lastindex(ks)
            sa, oa = scints[ks[ii]]
            sb, ob = scints[ks[jj]]
            if obb_overlap(oa, ob)
                push!(overlaps, ScintOverlap(sa.name, sb.name, mother))
            end
        end
    end
    return overlaps
end

"""
    check_geometry(m::GeometryManifest; clearance=1.0, arc_segments=24, on_clash=:error) -> m

Run the pre-flight geometry check on a manifest and return it (so the call
chains: `write_manifest(path, check_geometry(build_manifest(spec)))`).

Two physics-affecting overlaps are checked: fibre-fibre ([`fiber_clashes`](@ref)
— wrong light propagation) and scintillator-scintillator ([`scint_overlaps`](@ref)
— mis-attributed energy deposits). `on_clash` selects what happens when either
is found: `:error` throws (default — an overlap is a silently-wrong-physics bug),
`:warn` logs and continues, `:silent` just returns the manifest.
"""
function check_geometry(m::GeometryManifest; clearance=1.0,
                        arc_segments::Int=24, on_clash::Symbol=:error)
    on_clash in (:error, :warn, :silent) ||
        error("check_geometry: on_clash must be :error, :warn or :silent")
    clearance = _to_mm(clearance)                            # accept Real (mm) or a unitful length
    clashes  = fiber_clashes(m; clearance=clearance, arc_segments=arc_segments)
    overlaps = scint_overlaps(m)
    if isempty(clashes) && isempty(overlaps)
        return m
    end
    lines = String["geometry check failed:"]
    if !isempty(clashes)
        push!(lines, "  $(length(clashes)) fibre overlap(s) " *
                     "(clearance = $(clearance) mm):")
        for c in clashes
            push!(lines, "    $(c.a) <-> $(c.b) in '$(c.mother)'  " *
                         "min centreline distance = $(round(c.distance, digits=4)) mm")
        end
    end
    if !isempty(overlaps)
        push!(lines, "  $(length(overlaps)) scintillator overlap(s):")
        for o in overlaps
            push!(lines, "    $(o.a) <-> $(o.b) in '$(o.mother)'")
        end
    end
    msg = join(lines, "\n")
    on_clash === :error && error(msg)
    on_clash === :warn  && @warn msg
    return m
end
