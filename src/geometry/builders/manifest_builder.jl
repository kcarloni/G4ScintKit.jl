# manifest_builder.jl
#
# `ManifestBuilder` — a mutable accumulator that makes writing a new
# `build_manifest(::DetectorSpec)` declarative. It removes the three jobs every
# spec otherwise re-derives by hand:
#
#   * material resolution   — `ResolvedMaterials`, done once at construction;
#   * placement ordering    — the order you call `add_*!` *is* the manifest
#                             placement order (geometry-critical: GODDESS
#                             wrapping cuts are construction-order dependent);
#   * cross-reference wiring — names are registered on add, and every
#                             mother / reference / fiber / scint reference is
#                             validated in `to_manifest`, so a typo fails loudly
#                             instead of silently producing wrong geometry.
#
# The layout trig — the actual design — still lives in `build_manifest`. B1/B2
# predate this layer and keep their hand-written form; the builder targets new
# designs.
#
# Assumes the manifest structs (manifest.jl), the material resolution +
# `ResolvedMaterials` (materials.jl), the routers (dubins.jl, fiber_routing.jl)
# and `check_geometry` (geometry_check.jl) are already in module scope.

# ---------------------------------------------------------------------------
#  The builder
# ---------------------------------------------------------------------------

"""
    ManifestBuilder(spec)
    ManifestBuilder(materials::ResolvedMaterials)

A mutable accumulator for a [`GeometryManifest`](@ref). Add placements with the
`add_*!` functions (called in geometry-construction order), then call
[`to_manifest`](@ref) to validate and emit the manifest.
"""
mutable struct ManifestBuilder
    placements::Vector{PlacementEntry}
    materials::ResolvedMaterials
    names::Set{String}
    next_loop_id::Int
    casing::CasingSpec       # default CasingSpec() = disabled; set via add_casing!
end

ManifestBuilder(m::ResolvedMaterials) =
    ManifestBuilder(PlacementEntry[], m, Set{String}(), 0, CasingSpec())
ManifestBuilder(spec) = ManifestBuilder(ResolvedMaterials(spec))

# --- internal helpers ------------------------------------------------------

const _IDENTITY_ROT = (1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)

_ntup3(v)::NTuple{3,Float64} = (Float64(v[1]), Float64(v[2]), Float64(v[3]))
_ntup9(v)::NTuple{9,Float64} = ntuple(i -> Float64(v[i]), 9)

# Length-tuple form: each element may be a number (mm) or a unitful length;
# `_ntup3` stays for unitless 3-vectors (directions / bend axes).
_dims3(v) = (_to_mm(v[1]), _to_mm(v[2]), _to_mm(v[3]))

# Resolve a material argument: a `Symbol` keys into the builder's bundle, a
# string is taken as a literal `.properties` path.
function _matfile(b::ManifestBuilder, m::Symbol)
    m === :scint  && return b.materials.scint
    m === :wrap   && return b.materials.wrap
    m === :wls    && return b.materials.wls
    m === :cement && return b.materials.cement
    error("ManifestBuilder: unknown material key :$m " *
          "(expected :scint, :wrap, :wls or :cement)")
end
_matfile(::ManifestBuilder, m::AbstractString) = String(m)

# Register a placement name; reject empties and duplicates up front.
function _register!(b::ManifestBuilder, name::AbstractString)
    isempty(name) && error("ManifestBuilder: placement name must be non-empty")
    name in b.names &&
        error("ManifestBuilder: duplicate placement name '$name'")
    push!(b.names, String(name))
    return String(name)
end

# Lower a fibre's two endpoint `G4Coordinate`s to the entry's `reference`
# field. Both endpoints must share one frame; that frame becomes `reference`,
# or "" when it coincides with the placement `mother` (GODDESS treats an unset
# reference as the mother volume — see docs/placement_rules.md).
function _fiber_frame(start::G4Coordinate, stop::G4Coordinate,
                      mother::AbstractString, name)
    start.ref == stop.ref || error(
        "add_fiber_*!: fibre '$name' start is in frame '$(start.ref)' but " *
        "stop is in frame '$(stop.ref)' — both endpoints must share one frame")
    return start.ref == mother ? "" : start.ref
end

# Routed fibres carry their endpoints as G4DirectedPoints — same frame logic.
_fiber_frame(start::G4DirectedPoint, stop::G4DirectedPoint,
             mother::AbstractString, name) =
    _fiber_frame(start.point, stop.point, mother, name)

# Derive `fiber_is_base` from the frame the SiPM's `coupling_pos` is given in:
# the optical coupling is referenced to its base volume, which is either the
# fibre or the SiPM itself.
function _sipm_base(coupling_ref::AbstractString, fiber::AbstractString,
                    name::AbstractString)
    coupling_ref == fiber && return true
    coupling_ref == name  && return false
    error("add_sipm!: SiPM '$name' coupling_pos is in frame '$coupling_ref', " *
          "but the optical coupling must be referenced to either the fibre " *
          "('$fiber') or the SiPM itself ('$name')")
end

"""
    new_loop_id!(b::ManifestBuilder) -> Int

Allocate a fresh `loop_id`. Fibre segments sharing a `loop_id` are treated by
[`check_geometry`](@ref) as one physical fibre and never flagged against each
other — use one id per multi-segment fibre. [`add_fiber_routed!`](@ref) calls this
automatically.
"""
function new_loop_id!(b::ManifestBuilder)
    id = b.next_loop_id
    b.next_loop_id += 1
    return id
end

# ---------------------------------------------------------------------------
#  add_*! — append placements in construction order
# ---------------------------------------------------------------------------

"""
    add_scint!(b; name, dims, pos=G4Coordinate((0,0,0),"world"), rot=I,
               material=:scint, sensitive=true, g4name="") -> ScintEntry

Append a scintillator tile. `dims` is the full (x,y,z) size. `pos` is a
[`G4Coordinate`](@ref); its reference frame is the tile's mother volume
(`"world"` or an earlier scintillator). `material` is a key into the builder's
resolved materials (`:scint`/`:wrap`/`:wls`/`:cement`) or a literal
`.properties` path.
"""
function add_scint!(b::ManifestBuilder; name, dims,
                    pos::G4Coordinate=G4Coordinate((0.0, 0.0, 0.0), "world"),
                    rot=_IDENTITY_ROT, material=:scint,
                    sensitive::Bool=true, g4name::AbstractString="")
    _register!(b, name)
    entry = ScintEntry(
        name = String(name), g4name = String(g4name),
        dims = _dims3(dims), pos = pos.pos, rot = _ntup9(rot),
        mother = pos.ref, material_file = _matfile(b, material),
        sensitive = sensitive)
    push!(b.placements, entry)
    return entry
end

"""
    add_fiber_straight!(b; name, start::G4Coordinate, stop::G4Coordinate,
                        mother="world", glued=false, glue_profile="round",
                        end_reflectivity=NaN, loop_id=-1, material=:wls,
                        glue_material=:cement) -> FiberEntry

Append a straight fibre segment. `start`/`stop` are [`G4Coordinate`](@ref)s and
must share one reference frame; that frame becomes the entry's `reference`
(GODDESS frame for the endpoints). `mother` is the volume the fibre is placed
into and clipped against.
"""
function add_fiber_straight!(b::ManifestBuilder; name, start::G4Coordinate,
                             stop::G4Coordinate, mother="world",
                             glued::Bool=false, glue_profile="round",
                             start_reflectivity::Real=NaN,
                             end_reflectivity::Real=NaN, loop_id::Integer=-1,
                             material=:wls, glue_material=:cement)
    _register!(b, name)
    mother = _volname(mother)                        # accept a name or an entry
    entry = FiberEntry(
        name = String(name), kind = "straight", mother = mother,
        start = start.pos, stop = stop.pos,
        material_file = _matfile(b, material),
        reference = _fiber_frame(start, stop, mother, name),
        glued = glued,
        glue_file    = glued ? _matfile(b, glue_material) : "",
        glue_profile = glued ? String(glue_profile) : "",
        start_reflectivity = Float64(start_reflectivity),
        end_reflectivity = Float64(end_reflectivity),
        loop_id = Int(loop_id))
    push!(b.placements, entry)
    return entry
end

"""
    add_fiber_bent!(b; name, start::G4Coordinate, stop::G4Coordinate,
                    bend_angle, bend_axis, mother="world", glued=false,
                    glue_profile="round", end_reflectivity=NaN, loop_id=-1,
                    material=:wls, glue_material=:cement) -> FiberEntry

Append a bent (circular-arc) fibre segment. `start`/`stop` are
[`G4Coordinate`](@ref)s sharing one reference frame (see
[`add_fiber_straight!`](@ref)). `bend_angle` is in radian, `bend_axis` a
3-vector.
"""
function add_fiber_bent!(b::ManifestBuilder; name, start::G4Coordinate,
                         stop::G4Coordinate, bend_angle, bend_axis,
                         mother="world", glued::Bool=false,
                         glue_profile="round",
                         start_reflectivity::Real=NaN,
                         end_reflectivity::Real=NaN,
                         loop_id::Integer=-1, material=:wls,
                         glue_material=:cement)
    _register!(b, name)
    mother = _volname(mother)                        # accept a name or an entry
    entry = FiberEntry(
        name = String(name), kind = "bent", mother = mother,
        start = start.pos, stop = stop.pos,
        bend_angle = Float64(bend_angle), bend_axis = _ntup3(bend_axis),
        material_file = _matfile(b, material),
        reference = _fiber_frame(start, stop, mother, name),
        glued = glued,
        glue_file    = glued ? _matfile(b, glue_material) : "",
        glue_profile = glued ? String(glue_profile) : "",
        start_reflectivity = Float64(start_reflectivity),
        end_reflectivity = Float64(end_reflectivity),
        loop_id = Int(loop_id))
    push!(b.placements, entry)
    return entry
end

"""
    add_fiber_routed!(b; name, start::G4DirectedPoint, stop::G4DirectedPoint,
               min_radius, plane_normal=nothing, mother="world", glued=false,
               glue_profile="", material=:wls, glue_material=:cement,
               loop_id=new_loop_id!(b)) -> Vector{FiberEntry}

Route a smooth bounded-curvature fibre between two endpoints (see
[`route_fiber`](@ref)) and append every segment in order. Both endpoints must
share a reference frame ([`G4DirectedPoint`](@ref) carries it); that frame
becomes the entries' `reference`. The whole run shares one `loop_id`
(auto-allocated by default). Returns the appended segment entries.
"""
function add_fiber_routed!(b::ManifestBuilder; name, start::G4DirectedPoint,
                    stop::G4DirectedPoint, min_radius,
                    plane_normal=nothing, mother="world",
                    glued::Bool=false, glue_profile="",
                    start_reflectivity::Real=NaN,
                    end_reflectivity::Real=NaN,
                    material=:wls, glue_material=:cement,
                    loop_id::Integer=new_loop_id!(b))
    mother = _volname(mother)                        # accept a name or an entry
    reference = _fiber_frame(start, stop, mother, name)
    segs = route_fiber(start, stop; min_radius = min_radius,
                       plane_normal = plane_normal)
    entries = fiber_entries(segs; name = String(name),
        material_file = _matfile(b, material), mother = mother,
        reference = reference, glued = glued,
        glue_file = glued ? _matfile(b, glue_material) : "",
        glue_profile = glued ? String(glue_profile) : "",
        start_reflectivity = start_reflectivity,
        end_reflectivity = end_reflectivity,
        loop_id = Int(loop_id))
    for e in entries
        _register!(b, e.name)
        push!(b.placements, e)
    end
    return entries
end

"""
    add_fiber_path!(b; name, waypoints, min_radius, plane_normal=nothing,
               material=:wls, glued=false, glue_profile="round",
               glue_material=:cement, loop_id=new_loop_id!(b),
               continuity_tol=1e-3) -> Vector{FiberEntry}

Thread one continuous fibre through `waypoints` and append it as a `loop_id`-
tagged run of `FiberEntry`s. Each waypoint is a [`G4Coordinate`](@ref) (a free
point) or a [`G4DirectedPoint`](@ref) (a point with a pinned tangent), in any
reference frame. The route (see [`route_path`](@ref)) is: a gap between two
pinned waypoints is a Dubins curve; a run of free gaps is straight segments
with a `min_radius` fillet at each corner. Each segment is assigned its mother
volume by midpoint containment — straights that cross a scintillator wall are
left whole and placed in that scintillator, matching B1/B2's pattern of an
inner fibre extending past the wall (GODDESS's G4Fibre handles the boolean
carve, and the lossy fibre endcaps land outside the scint instead of at the
wall where the wrapping would otherwise trap light). Assumes a single straight
passes through at most one scintillator.

Each piece's endpoints are emitted in its **mother's** local frame (with
`reference=""`, so GODDESS treats the mother as the reference). This keeps
endpoints frame-correct for non-origin scintillators and avoids GODDESS's
"reference is neither mother nor a child of mother" warning. Returns the
appended segment entries. Errors on a tangent kink at a pinned waypoint.
"""
function add_fiber_path!(b::ManifestBuilder; name, waypoints, min_radius,
                         plane_normal=nothing, material=:wls,
                         glued::Bool=false, glue_profile="round",
                         glue_material=:cement,
                         start_reflectivity::Real=NaN,
                         end_reflectivity::Real=NaN,
                         loop_id::Integer=new_loop_id!(b),
                         continuity_tol::Real=1e-3)
    resolved = _resolve_path_waypoints(b.placements, waypoints)
    segs = route_path(resolved; min_radius = min_radius,
                      plane_normal = plane_normal, continuity_tol = continuity_tol)
    pieces = _split_path_mothers(segs, _scint_obbs(b.placements))
    isempty(pieces) && error("add_fiber_path!: routed path has no segments")

    mat   = _matfile(b, material)
    gfile = glued ? _matfile(b, glue_material) : ""
    gprof = glued ? String(glue_profile) : ""
    npieces = length(pieces)
    # World-frame endpoints from the router are lowered into each piece's
    # mother frame so the FiberEntry's `reference` can stay empty (i.e. the
    # mother *is* the reference). Necessary for non-origin scintillator
    # mothers, and silences GODDESS's "reference is neither mother nor a
    # child of mother" warning even for at-origin ones.
    fc = _FrameCache(GeometryManifest(placements = b.placements))
    entries = FiberEntry[]
    for (k, (seg, mother)) in enumerate(pieces)
        ename = "$(String(name))_$(k)"
        _register!(b, ename)
        if mother == "world"
            sstart, sstop, sbaxis = Tuple(seg.start), Tuple(seg.stop), Tuple(seg.bend_axis)
        else
            inv_tf = inv(_world_transform(fc, mother))
            sstart = Tuple(inv_tf(seg.start))
            sstop  = Tuple(inv_tf(seg.stop))
            sbaxis = Tuple(inv_tf.rotation * seg.bend_axis)   # direction: rotate only
        end
        # Reflectivities apply only at the routed fibre's physical endpoints:
        # the first segment's `start` and the last segment's `stop`. Interior
        # joins are roughened scatter junctions, not physical fibre ends.
        entry = FiberEntry(
            name = ename,
            kind = _kindstr(seg.kind),
            mother = mother,
            start = sstart, stop = sstop,
            bend_angle = seg.bend_angle, bend_axis = sbaxis,
            material_file = mat, reference = "",        # mother is the reference
            glued = glued, glue_file = gfile, glue_profile = gprof,
            start_reflectivity = k == 1       ? Float64(start_reflectivity) : NaN,
            end_reflectivity   = k == npieces ? Float64(end_reflectivity)   : NaN,
            loop_id = Int(loop_id))
        push!(b.placements, entry)
        push!(entries, entry)
    end
    return entries
end

"""
    add_wrapping!(b; scint, g4name="", material=:wrap) -> WrapEntry

Append a reflective wrapping around an already-added scintillator tile. The
`cut` list is left empty so `PlaceManifest` auto-derives the wrapping cuts.
"""
function add_wrapping!(b::ManifestBuilder; scint, g4name::AbstractString="",
                   material=:wrap)
    entry = WrapEntry(
        scint = _volname(scint), g4name = String(g4name),
        material_file = _matfile(b, material))
    push!(b.placements, entry)
    return entry
end

"""
    add_sipm!(b; name, fiber, face_dir, rel_pos::G4Coordinate, edge_length,
              coupling_normal, coupling_pos::G4Coordinate, coupling_width) -> SipmEntry

Append a SiPM photon detector and its optical coupling to `fiber`. `rel_pos` is
a [`G4Coordinate`](@ref) whose frame is the SiPM's reference volume (`"world"`
or an earlier scintillator). `coupling_pos` is a `G4Coordinate` whose frame is
the optical coupling's base volume — either `fiber` or the SiPM itself
(`name`); that choice sets `fiber_is_base`. `face_dir` and `coupling_normal`
are direction vectors in those respective frames.
"""
function add_sipm!(b::ManifestBuilder; name, fiber, face_dir,
                   rel_pos::G4Coordinate, edge_length, coupling_normal,
                   coupling_pos::G4Coordinate, coupling_width)
    _register!(b, name)
    fiber = _volname(fiber)                          # accept a name or an entry
    entry = SipmEntry(
        name = String(name), ref_volume = rel_pos.ref,
        face_dir = _ntup3(face_dir), rel_pos = rel_pos.pos,
        edge_length = _to_mm(edge_length), fiber = fiber,
        coupling_normal = _ntup3(coupling_normal),
        coupling_pos = coupling_pos.pos,
        coupling_width = _to_mm(coupling_width),
        fiber_is_base = _sipm_base(coupling_pos.ref, fiber, String(name)))
    push!(b.placements, entry)
    return entry
end

# ---------------------------------------------------------------------------
#  Module extent -> casing
# ---------------------------------------------------------------------------

# Axis-aligned bounding box (mm) of every world-frame scintillator tile and
# fibre endpoint. Scintillator corners are rotated by the tile's `rot`.
# Mother-relative placements are skipped — their coordinates are not in the
# world frame, so they cannot be bounded without resolving the mother chain.
function _module_extent(placements)
    lo = MVector(Inf, Inf, Inf)
    hi = MVector(-Inf, -Inf, -Inf)
    grow!(p) = @inbounds for k in 1:3
        lo[k] = min(lo[k], p[k])
        hi[k] = max(hi[k], p[k])
    end
    for e in placements
        if e isa ScintEntry && e.mother == "world"
            R = _rotmatrix(e.rot)                  # row-major NTuple{9} -> matrix
            c = SVector{3,Float64}(e.pos...)
            h = SVector{3,Float64}(e.dims...) ./ 2
            for sx in (-1.0, 1.0), sy in (-1.0, 1.0), sz in (-1.0, 1.0)
                grow!(c + R * SVector(sx * h[1], sy * h[2], sz * h[3]))
            end
        elseif e isa FiberEntry && e.mother == "world"
            grow!(SVector{3,Float64}(e.start...))
            grow!(SVector{3,Float64}(e.stop...))
        end
    end
    isfinite(lo[1]) ||
        error("ManifestBuilder: no scintillators or fibres to bound")
    return SVector{3,Float64}(lo), SVector{3,Float64}(hi)
end

"""
    casing_from_extent(b; aluminum_thickness=0, lead_thickness=0, num_bars=0,
                       bar_width=0, scinti_z=0, margin=0) -> CasingSpec

Build a [`CasingSpec`](@ref) whose module bounding box is the axis-aligned
extent of every scintillator and fibre added so far (plus `margin`).
`module_half_x`/`module_half_z` assume the module is centred on the origin in x
and z (as B1/B2 are); `module_min_y`/`module_max_y` are the true y bounds. The
remaining fields are semantic and passed straight through.
"""
function casing_from_extent(b::ManifestBuilder; aluminum_thickness::Real=0.0,
                            lead_thickness::Real=0.0, num_bars::Integer=0,
                            bar_width::Real=0.0, scinti_z::Real=0.0,
                            margin::Real=0.0)
    lo, hi = _module_extent(b.placements)
    return CasingSpec(
        module_half_x = max(abs(lo[1]), abs(hi[1])) + margin,
        module_min_y  = lo[2] - margin,
        module_max_y  = hi[2] + margin,
        module_half_z = max(abs(lo[3]), abs(hi[3])) + margin,
        aluminum_thickness = Float64(aluminum_thickness),
        lead_thickness     = Float64(lead_thickness),
        num_bars  = Int(num_bars),
        bar_width = Float64(bar_width),
        scinti_z  = Float64(scinti_z))
end

"""
    add_casing!(b; aluminum_thickness=0u"mm", lead_thickness=0u"mm",
                lead_overhang_x=0u"mm", lead_overhang_z=0u"mm",
                margin=0u"mm") -> CasingSpec

Compute the [`CasingSpec`](@ref) for the detector module already placed in `b`
and store it on the builder. The module bounding box is auto-derived from the
axis-aligned extent of every world-frame scintillator and fibre placement
(plus `margin` on every side).

Set `aluminum_thickness > 0` to enable the aluminum box (it sits one fixed mm
of air gap outside the module bbox, with this thickness as its wall).

Set `lead_thickness > 0` to enable the lead roof sheet sitting on top of the
aluminum box. By default the sheet matches the module footprint in x and z;
`lead_overhang_x` / `lead_overhang_z` add the same amount of extra coverage
on each side of the respective axis (B1/B2 use a non-zero overhang so the
lead extends beyond the bars).

The casing is picked up automatically by `to_manifest(b)` (which now defaults
its `casing` kwarg to `b.casing`). Pass an explicit `casing = …` to
`to_manifest` to override.
"""
function add_casing!(b::ManifestBuilder; aluminum_thickness=0.0u"mm",
                     lead_thickness=0.0u"mm",
                     lead_overhang_x=0.0u"mm", lead_overhang_z=0.0u"mm",
                     margin=0.0u"mm")
    aluminum_thickness = _to_mm(aluminum_thickness)
    lead_thickness     = _to_mm(lead_thickness)
    lead_overhang_x    = _to_mm(lead_overhang_x)
    lead_overhang_z    = _to_mm(lead_overhang_z)
    margin             = _to_mm(margin)

    lo, hi = _module_extent(b.placements)
    module_half_x = max(abs(lo[1]), abs(hi[1])) + margin
    module_half_z = max(abs(lo[3]), abs(hi[3])) + margin

    # Encode the lead sheet's half-extents into the C++ format. Only meaningful
    # when lead_thickness > 0 — the C++ code skips the sheet otherwise, so the
    # values are inert when lead is disabled.
    lead_half_x = module_half_x + lead_overhang_x
    lead_half_z = module_half_z + lead_overhang_z
    num_bars  = 1
    bar_width = lead_half_x                   # sheet_half_x = num_bars * bar_width = lead_half_x
    scinti_z  = 2 * lead_half_z               # sheet_half_z = scinti_z / 2 = lead_half_z

    b.casing = CasingSpec(
        module_half_x      = module_half_x,
        module_min_y       = lo[2] - margin,
        module_max_y       = hi[2] + margin,
        module_half_z      = module_half_z,
        aluminum_thickness = aluminum_thickness,
        lead_thickness     = lead_thickness,
        num_bars  = num_bars,
        bar_width = bar_width,
        scinti_z  = scinti_z)
    return b.casing
end

# ---------------------------------------------------------------------------
#  Reference validation + finalisation
# ---------------------------------------------------------------------------

# Check that `value` names a volume placed *before* the referencing entry:
# `seen` = names of that kind placed so far, `known` = all names of that kind.
# Distinguishes an unknown name (a typo) from one merely placed later.
function _check_ref(kind, owner, field, value, seen, known)
    value in seen && return nothing
    value in known && error(
        "ManifestBuilder: $kind '$owner' references $field='$value', which is " *
        "placed after it — a referenced volume must be added before the entry " *
        "that references it")
    error("ManifestBuilder: $kind '$owner' references $field='$value', " *
          "which is not a placed volume")
end

# Validate every mother / reference / fiber / scint reference: each must name a
# volume placed *earlier* in the list. Catches both typos and out-of-order adds
# — the manifest order is the C++ construction order, so a forward reference
# would fail in `PlaceManifest`.
function _validate_references(placements)
    all_scints = Set{String}(p.name for p in placements if p isa ScintEntry)
    all_fibers = Set{String}(p.name for p in placements if p isa FiberEntry)
    scints = Set{String}()                          # scints placed so far
    fibers = Set{String}()                          # fibres placed so far
    for p in placements
        if p isa ScintEntry
            p.mother == "world" || _check_ref(
                "scint", p.name, "mother", p.mother, scints, all_scints)
            push!(scints, p.name)
        elseif p isa FiberEntry
            p.mother == "world" || _check_ref(
                "fiber", p.name, "mother", p.mother, scints, all_scints)
            isempty(p.reference) || p.reference == "world" || _check_ref(
                "fiber", p.name, "reference", p.reference, scints, all_scints)
            push!(fibers, p.name)
        elseif p isa WrapEntry
            _check_ref("wrap", p.scint, "scint", p.scint, scints, all_scints)
        elseif p isa SipmEntry
            p.ref_volume == "world" || _check_ref(
                "sipm", p.name, "ref_volume", p.ref_volume, scints, all_scints)
            _check_ref("sipm", p.name, "fiber", p.fiber, fibers, all_fibers)
        end
    end
    return nothing
end

"""
    to_manifest(b::ManifestBuilder; setup_label="", casing=b.casing,
                check=:error, clearance=<fibre diameter>) -> GeometryManifest

Finalise the builder into a [`GeometryManifest`](@ref). Cross-references are
validated first — every `mother`/`reference`/`fiber`/`scint` must name a volume
placed earlier in the list — then [`check_geometry`](@ref) runs with
`on_clash=check` — `:error` (default) throws on a geometry clash, `:warn` logs
and continues, `:silent` skips the report.

`casing` defaults to the builder's stored casing (set via [`add_casing!`](@ref);
disabled when not called). Pass an explicit `CasingSpec` to override.
`clearance` is the fibre-fibre minimum-separation threshold; it defaults to the
actual fibre diameter read from the builder's WLS material (rather than a
hard-coded 1 mm), so a design with a thicker fibre is checked against its true
size.
"""
function to_manifest(b::ManifestBuilder; setup_label::AbstractString="",
                     casing::CasingSpec=b.casing, check::Symbol=:error,
                     clearance=fiber_cross_section(b.materials.wls).width)
    _validate_references(b.placements)
    m = GeometryManifest(setup_label = String(setup_label),
                         placements = b.placements, casing = casing)
    check_geometry(m; on_clash = check, clearance = clearance)
    return m
end
