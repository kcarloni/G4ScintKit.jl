# detector_assembly.jl
#
# Layer-2 composite helpers: small, reusable pieces of detector-geometry
# vocabulary built on the `ManifestBuilder` primitives (manifest_builder.jl).
# Where the Layer-1 `add_*!` calls place one manifest entry each, these encode
# recurring physical patterns — a lattice of bars, a SiPM coupled to the end of
# a straight readout fibre — so a `build_manifest` reads closer to the design.
#
# This file grows as new designs reveal shared patterns; it is deliberately
# small to start.

"""
    bar_lattice_spacing(num_bars, pitch) -> Vector{Float64}

The x-coordinates of `num_bars` parallel bars, centred on the origin with
centre-to-centre spacing `pitch` (the full bar width including inter-bar gaps).
Mirrors the B2/B3 bar layout — for an even `num_bars` the bars straddle the
origin; for an odd count the layout matches B2's `(num_bars ÷ 2)` convention.
"""
function bar_lattice_spacing(num_bars::Integer, pitch::Real)
    num_bars >= 1 ||
        error("bar_lattice_spacing: num_bars must be >= 1, got $num_bars")
    x0 = -(num_bars ÷ 2) * pitch + pitch / 2
    return [Float64(x0 + i * pitch) for i in 0:(num_bars - 1)]
end

"""
    add_scint_row!(b; name_prefix="scint", num_bars, inter_bar_spacing, dims,
                   center=G4Coordinate((0,0,0), "world"), axis=:x,
                   rot=I, material=:scint, sensitive=true,
                   g4name_prefix="scintillator") -> Vector{ScintEntry}

Place `num_bars` identical scintillator bars in a 1D row along `axis`
(`:x`, `:y`, or `:z`), centred on `center`. Adjacent bars are separated by
`inter_bar_spacing` — the gap between facing surfaces along `axis`; the
centre-to-centre pitch is `dims[axis] + inter_bar_spacing` (so for an even
`num_bars` the row straddles `center`, matching [`bar_lattice_spacing`](@ref)).
Bars are named `"<name_prefix>_0"`, `"<name_prefix>_1"`, … with the same
0-indexed convention applied to `g4name_prefix`.

Returns the placed `ScintEntry`s in row order so the caller can index into
them for per-bar fibre / wrapping passes.
"""
function add_scint_row!(b::ManifestBuilder;
                        name_prefix::AbstractString = "scint",
                        num_bars::Integer,
                        inter_bar_spacing,
                        dims,
                        center::G4Coordinate = G4Coordinate((0.0, 0.0, 0.0), "world"),
                        axis::Symbol = :x,
                        rot = _IDENTITY_ROT,
                        material = :scint,
                        sensitive::Bool = true,
                        g4name_prefix::AbstractString = "scintillator")
    num_bars >= 1 ||
        error("add_scint_row!: num_bars must be >= 1, got $num_bars")
    axis_idx = axis === :x ? 1 :
               axis === :y ? 2 :
               axis === :z ? 3 :
               error("add_scint_row!: axis must be :x, :y, or :z, got :$axis")

    dims_mm = _dims3(dims)
    pitch   = dims_mm[axis_idx] + _to_mm(inter_bar_spacing)
    offsets = bar_lattice_spacing(num_bars, pitch)

    e_axis = ntuple(i -> i == axis_idx ? 1.0 : 0.0, 3)

    entries = ScintEntry[]
    for (i, off) in enumerate(offsets)
        pos = center + G4Vector(off .* e_axis, center.ref)
        push!(entries, add_scint!(b;
            name      = string(name_prefix, "_", i ),
            dims      = dims_mm,
            pos       = pos,
            rot       = rot,
            material  = material,
            sensitive = sensitive,
            g4name    = string(g4name_prefix, "_", i )))
    end
    return entries
end

# Accept either a `G4Direction` (with a frame check against the host's frame)
# or a bare 3-vector for an in-plane / direction argument; return a Float64
# SVector. Lets `bundle_fiber_endpoints` (and similar helpers) take the same
# kwargs in either the frame-aware or the legacy bare-vector form.
function _coerce_inplane_vec(v::G4Direction, host_ref::AbstractString, name::AbstractString)
    v.ref == host_ref || error(
        "bundle_fiber_endpoints: $name is in frame '$(v.ref)' but plane_center " *
        "is in frame '$host_ref' — both must share a reference frame")
    return v.vec
end
_coerce_inplane_vec(v, ::AbstractString, ::AbstractString) =
    SVector{3,Float64}(Float64(v[1]), Float64(v[2]), Float64(v[3]))

# Map -0.0 to +0.0 so a direction/coupling vector serialises cleanly
# (e.g. "0,0,-1" rather than "-0,-0,-1").
_clean_zeros(t::NTuple{3,Float64}) =
    (t[1] == 0.0 ? 0.0 : t[1],
     t[2] == 0.0 ? 0.0 : t[2],
     t[3] == 0.0 ? 0.0 : t[3])

"""
    add_inline_sipm!(b; name, fiber::FiberEntry, edge_length, coupling_width,
                     sipm_end=:stop) -> SipmEntry

Place a SiPM flush against one end of a *straight* fibre and couple it to that
fibre. `fiber` is the [`FiberEntry`](@ref) returned by `add_fiber_straight!` /
`add_fiber_bent!`; `sipm_end` picks the readout end — `:stop` (default) puts
the SiPM at `fiber.stop`, `:start` at `fiber.start`. The SiPM sits
`coupling_width` beyond that endpoint along the fibre axis, facing back down
it; the coupling is centred on the fibre's local frame. The SiPM's reference
volume is the fibre's endpoint frame (`fiber.reference` if set, otherwise
`fiber.mother`).

This packages the coupling geometry — `face_dir`, `rel_pos`, `coupling_normal`,
`coupling_pos` — that a spec would otherwise hand-derive. The optical coupling
is always referenced to the fibre (the "inline" case); if you need the SiPM as
the coupling base instead, build it with [`add_sipm!`](@ref) directly (B2 does
this).
"""
function add_inline_sipm!(b::ManifestBuilder; name, fiber::FiberEntry,
                          edge_length, coupling_width,
                          sipm_end::Symbol = :stop)
    edge_length    = _to_mm(edge_length)             # accept unitful lengths
    coupling_width = _to_mm(coupling_width)
    sipm_end in (:start, :stop) || error(
        "add_inline_sipm!: sipm_end must be :start or :stop, got :$sipm_end")
    far_pos, sipm_pos = sipm_end === :stop ? (fiber.start, fiber.stop) :
                                             (fiber.stop,  fiber.start)
    frame = isempty(fiber.reference) ? fiber.mother : fiber.reference
    s = SVector{3,Float64}(far_pos)
    e = SVector{3,Float64}(sipm_pos)
    seg = e - s
    L = norm(seg)
    L > 0 || error("add_inline_sipm!: readout fibre has zero length")
    axis = seg / L

    # Optical coupling is referenced to the fibre (fiber_is_base=true). GODDESS
    # builds the fibre with local +z along (StartPoint -> EndPoint), origin at
    # the fibre midpoint (FibreConstructor::GenerateTransformation in
    # G4Fibre.cc). The SiPM sits past `sipm_end` along the fibre axis, so the
    # coupling centre lies on that side of the fibre's local origin:
    #   sipm_end = :stop  -> SiPM on the EndPoint side  -> coupling at local +z
    #   sipm_end = :start -> SiPM on the StartPoint side -> coupling at local -z
    # The coupling's outward normal (toward the SiPM) points the same way.
    sign = sipm_end === :stop ? 1.0 : -1.0
    coupling_pos_local    = (0.0, 0.0, sign * (L / 2 + coupling_width / 2))
    coupling_normal_local = (0.0, 0.0, sign)

    return add_sipm!(b; name = name, fiber = fiber.name,
        face_dir        = _clean_zeros(Tuple(-axis)),
        rel_pos         = G4Coordinate(_clean_zeros(Tuple(e + coupling_width * axis)), frame),
        edge_length     = edge_length,
        coupling_normal = coupling_normal_local,
        coupling_pos    = G4Coordinate(coupling_pos_local, fiber.name),
        coupling_width  = coupling_width)
end

"""
    add_bundle_sipm!(b; name, anchor_fiber::FiberEntry,
                    endpoints::AbstractVector{G4Coordinate}, face_dir,
                    edge_length, coupling_width,
                    sipm_thickness=0.5u"mm",
                    planarity_tol=1e-6) -> SipmEntry

Place a SiPM flush against the terminal plane of a fibre bundle. `endpoints`
is the list of bundle fibre tips (typically the value returned by
[`bundle_fiber_endpoints`](@ref)); the bundle's plane centre is taken as
their centroid. `face_dir` is a [`G4Direction`](@ref), `G4Vector`, or bare
3-vector pointing from the SiPM toward the bundle (the direction the fibre
tips are heading when they reach the plane). The SiPM sits one
`coupling_width` back from the centroid along `-face_dir`; the optical
coupling slab fills the gap and is parented to the SiPM
(`fiber_is_base=false`) because it physically spans several fibres.

`anchor_fiber` is the one fibre wired into the `SipmEntry` — the other bundle
fibres end at the bundle plane and rely on the coupling slab for optical
contact. `edge_length` is the SiPM's square edge size; the helper validates
that every endpoint sits within `edge_length/2` of the centroid in the
bundle plane (the inscribed-circle bound — orientation-independent), and
errors otherwise. `sipm_thickness` is the SiPM block's thickness along
`face_dir` and defaults to the GODDESS convention (0.5 mm).
`planarity_tol` (mm) bounds how far any endpoint may deviate from the bundle
plane along `face_dir`.

For a single straight fibre, prefer [`add_inline_sipm!`](@ref).
"""
function add_bundle_sipm!(b::ManifestBuilder; name, anchor_fiber::FiberEntry,
                          endpoints::AbstractVector{G4Coordinate}, face_dir,
                          edge_length, coupling_width,
                          sipm_thickness = 0.5u"mm",
                          planarity_tol::Real = 1e-6)
    isempty(endpoints) &&
        error("add_bundle_sipm!: endpoints must be non-empty")

    # All endpoints must share a reference frame.
    ref = first(endpoints).ref
    all(e -> e.ref == ref, endpoints) || error(
        "add_bundle_sipm!: every endpoint must share a reference frame " *
        "(got mixed frames among the given endpoints)")

    edge_length    = _to_mm(edge_length)
    coupling_width = _to_mm(coupling_width)
    sipm_half_t    = _to_mm(sipm_thickness) / 2

    # Coerce face_dir to G4Direction in the endpoints' frame (frame-checked).
    face = G4Direction(face_dir, ref)

    # Bundle plane centre is the centroid of the endpoints. Validate planarity
    # and the in-plane footprint against the SiPM face in one pass.
    n  = length(endpoints)
    cx = sum(e.pos[1] for e in endpoints) / n
    cy = sum(e.pos[2] for e in endpoints) / n
    cz = sum(e.pos[3] for e in endpoints) / n
    plane_center = G4Coordinate((cx, cy, cz), ref)

    half_edge = edge_length / 2
    centroid = SVector{3,Float64}(cx, cy, cz)
    for e in endpoints
        d        = SVector{3,Float64}(e.pos) - centroid
        along    = dot(d, face.vec)              # component along face_dir
        in_plane = d - along * face.vec
        r        = norm(in_plane)
        abs(along) > planarity_tol && error(
            "add_bundle_sipm!: endpoint $(e.pos) is $(round(abs(along), digits=6)) " *
            "mm out of the bundle plane (tol $planarity_tol mm) — the bundle " *
            "endpoints must be coplanar")
        r > half_edge && error(
            "add_bundle_sipm!: endpoint $(e.pos) is $r mm from the bundle " *
            "centroid in-plane, but edge_length/2 = $half_edge mm — " *
            "increase edge_length or check the bundle layout")
    end

    # The SiPM block sits one coupling_width back from the bundle plane,
    # opposite face_dir, so the coupling slab fills the gap.
    sipm_pos = plane_center - coupling_width * face

    # Slab parented to the SiPM: in SiPM-local frame, coupling sits at +z (the
    # face side, by GODDESS convention) with outward normal along -z. Mirrors
    # the B2 hand-rolled placement.
    coupling_pos_local    = (0.0, 0.0, 0.5 * coupling_width + sipm_half_t)
    coupling_normal_local = (0.0, 0.0, -1.0)

    return add_sipm!(b; name = name, fiber = anchor_fiber.name,
        face_dir        = _clean_zeros(Tuple(face.vec)),
        rel_pos         = sipm_pos,
        edge_length     = edge_length,
        coupling_normal = coupling_normal_local,
        coupling_pos    = G4Coordinate(coupling_pos_local, String(name)),
        coupling_width  = coupling_width)
end

"""
    bundle_fiber_endpoints(; num_fibers, pitch, plane_center::G4Coordinate,
                           plane_normal, axis=nothing, arrangement=:hex)
        -> Vector{G4Coordinate}

Lay out `num_fibers` co-planar points spaced by `pitch` (centre-to-centre) on
the plane through `plane_center` with normal `plane_normal`. Use this to place
the terminal endpoints of a fibre bundle against a shared optical-coupling
slab — the "bundle ferrule" pattern: N fibres feed one SiPM by ending flush
against one coupling slab; the SiPM is wired (via [`add_inline_sipm!`](@ref) /
[`add_sipm!`](@ref)) to a single "anchor" fibre, and the others ride along
geometrically.

`axis` is an in-plane reference direction (the line direction for `:line`, the
first hex-row direction for `:hex`). A component along `plane_normal` is
projected out. `nothing` (default) auto-picks a standard-basis direction
orthogonal to `plane_normal` — fine for `:hex` (orientation rarely matters),
but `:line` callers should usually pass `axis` explicitly to pin the line
direction. Returned `G4Coordinate`s carry `plane_center`'s reference frame.

Arrangements:
- `:line` — N points spaced along `axis`, centred on `plane_center`.
- `:hex`  — N points filling a hexagonal close-pack lattice inward-out
  (centre, then ring of 6 nearest neighbours, then ring of 12, …). `pitch`
  is centre-to-centre between nearest neighbours.
- `:grid` — explicit `rows × cols` rectangular lattice (square pitch in both
  axes). Columns run along `axis`; each column is a vertical stack of `rows`
  points along the perpendicular in-plane direction. Requires
  `num_fibers == rows * cols`; pass `rows` and `cols` as kwargs. Emitted
  column-major so per-column points are contiguous in the result — useful
  when you want adjacent return-vector entries to share a column (e.g. pairs
  to back-loop together).

The helper does not validate the bundle against the SiPM footprint or the
fibre-clash check — pair with [`check_geometry`](@ref) at finalisation
(`to_manifest` runs it by default).
"""
function bundle_fiber_endpoints(; num_fibers::Integer, pitch,
                                plane_center::G4Coordinate,
                                plane_normal, axis=nothing,
                                arrangement::Symbol = :hex,
                                rows::Integer = 0, cols::Integer = 0)
    num_fibers >= 1 ||
        error("bundle_fiber_endpoints: num_fibers must be >= 1, got $num_fibers")
    pitch_mm = _to_mm(pitch)

    n_hat = _coerce_inplane_vec(plane_normal, plane_center.ref, "plane_normal")

    # Default `axis`: pick the standard basis direction least-aligned with the
    # normal (so the projection always has a healthy in-plane magnitude).
    a = if axis === nothing
        k = argmin((abs(n_hat[1]), abs(n_hat[2]), abs(n_hat[3])))
        SVector(ntuple(i -> i == k ? 1.0 : 0.0, 3)...)
    else
        _coerce_inplane_vec(axis, plane_center.ref, "axis")
    end
    # Project onto the bundle plane; a parallel `axis` is ambiguous.
    u_raw = a - dot(a, n_hat) * n_hat
    un = norm(u_raw)
    un > 0 || error("bundle_fiber_endpoints: axis is parallel to plane_normal " *
                    "— pick an axis with a nonzero in-plane component")
    u_hat = u_raw / un
    v_hat = cross(n_hat, u_hat)               # already unit-length (n̂⊥û)

    base   = SVector{3,Float64}(plane_center.pos)
    ref    = plane_center.ref
    coords = G4Coordinate[]

    if arrangement === :line
        for k in 0:(num_fibers - 1)
            offset = (k - (num_fibers - 1) / 2) * pitch_mm
            p = base + offset * u_hat
            push!(coords, G4Coordinate(Tuple(p), ref))
        end
    elseif arrangement === :hex
        # Triangular lattice with axial basis a1=û, a2=½û+(√3/2)v̂; cell (i,j)
        # sits at distance² = i²+ij+j² (in units of pitch²) from the centre.
        # Collect cells inside a generous ring and sort by distance to fill
        # the bundle inward-out — exactly the close-pack ring expansion.
        a1 = u_hat
        a2 = 0.5 * u_hat + (sqrt(3.0) / 2.0) * v_hat
        rmax = ceil(Int, sqrt(Float64(num_fibers)))
        cells = Tuple{Int,Int,Int}[]
        for i in -rmax:rmax, j in -rmax:rmax
            push!(cells, (i*i + i*j + j*j, i, j))   # (dist², i, j) for stable sort
        end
        sort!(cells)
        for k in 1:num_fibers
            (_, i, j) = cells[k]
            p = base + (i * a1 + j * a2) * pitch_mm
            push!(coords, G4Coordinate(Tuple(p), ref))
        end
    elseif arrangement === :grid
        # Rectangular `rows × cols` lattice, square pitch in both axes. Columns
        # are arrayed along `axis` (so the column index runs along the long
        # axis); each column is a vertical stack of `rows` points along v_hat.
        # Emitted column-major — column 0's `rows` points first, then column 1,
        # … — so per-column pairs land contiguous in the returned vector, which
        # is what callers doing back-loop pairing want.
        (rows >= 1 && cols >= 1) || error(
            "bundle_fiber_endpoints: arrangement=:grid requires rows >= 1 and " *
            "cols >= 1 (got rows=$rows, cols=$cols)")
        rows * cols == num_fibers || error(
            "bundle_fiber_endpoints: arrangement=:grid requires " *
            "num_fibers == rows * cols (got num_fibers=$num_fibers, " *
            "rows=$rows, cols=$cols, rows*cols=$(rows*cols))")
        for c in 0:(cols - 1), r in 0:(rows - 1)
            du = (c - (cols - 1) / 2) * pitch_mm
            dv = (r - (rows - 1) / 2) * pitch_mm
            p  = base + du * u_hat + dv * v_hat
            push!(coords, G4Coordinate(Tuple(p), ref))
        end
    else
        error("bundle_fiber_endpoints: arrangement must be :line, :hex, or " *
              ":grid, got :$arrangement")
    end
    return coords
end
