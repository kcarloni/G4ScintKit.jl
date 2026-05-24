# display.jl
#
# Custom pretty-printers (`Base.show(io, ::MIME"text/plain", x)`) for the
# user-facing types defined in this package. The default Julia struct print
# dumps every field by name, which is unreadable for a manifest with dozens of
# placements; the methods here produce a compact, scannable REPL summary while
# leaving the one-arg `show(io, x)` mostly intact (so arrays, `@show`, and
# logging stay terse and unambiguous).
#
# All types referenced here are defined in earlier `include`d files (manifest,
# detector_spec, materials, frames, geometry_check, builders/*), so this file
# must be included last.

# ---------------------------------------------------------------------------
#  Number formatting helpers (local to display)
# ---------------------------------------------------------------------------

# Compact number rendering for tuples in summaries: drops trailing-zero noise
# but stays readable. NaN/Inf pass through the default formatter.
_fmt(x::Float64) = isfinite(x) ? @sprintf("%g", x) : string(x)
_fmt(x::Real)    = _fmt(Float64(x))

_fmt_tup(t::NTuple{3,<:Real}) = string("(", _fmt(t[1]), ", ", _fmt(t[2]), ", ", _fmt(t[3]), ")")

# ---------------------------------------------------------------------------
#  PlacementEntry — one-line forms (used by lists, arrays, @show)
# ---------------------------------------------------------------------------

function Base.show(io::IO, s::ScintEntry)
    print(io, "ScintEntry(\"", s.name, "\", dims=", _fmt_tup(s.dims),
          ", mother=\"", s.mother, "\"")
    s.sensitive && print(io, ", sensitive=true")
    print(io, ")")
end

function Base.show(io::IO, f::FiberEntry)
    print(io, "FiberEntry(\"", f.name, "\", ", f.kind,
          ", mother=\"", f.mother, "\"")
    isempty(f.reference) || print(io, ", reference=\"", f.reference, "\"")
    if f.kind == "bent"
        print(io, ", bend_angle=", _fmt(f.bend_angle))
    end
    f.loop_id >= 0 && print(io, ", loop_id=", f.loop_id)
    print(io, ")")
end

function Base.show(io::IO, w::WrapEntry)
    print(io, "WrapEntry(scint=\"", w.scint, "\", cut=",
          isempty(w.cut) ? "auto" : string("[", join(("\"$c\"" for c in w.cut), ", "), "]"),
          ")")
end

function Base.show(io::IO, sp::SipmEntry)
    print(io, "SipmEntry(\"", sp.name, "\", ref_volume=\"", sp.ref_volume,
          "\", fiber=\"", sp.fiber, "\")")
end

# ---------------------------------------------------------------------------
#  PlacementEntry — multi-line (REPL) forms
# ---------------------------------------------------------------------------

function Base.show(io::IO, ::MIME"text/plain", s::ScintEntry)
    println(io, "ScintEntry \"", s.name, "\"")
    isempty(s.g4name)     || println(io, "  g4name        = \"", s.g4name, "\"")
    println(io, "  dims          = ", _fmt_tup(s.dims), " mm")
    println(io, "  pos           = ", _fmt_tup(s.pos), " mm")
    s.rot == (1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0) ||
        println(io, "  rot           = ", s.rot)
    println(io, "  mother        = \"", s.mother, "\"")
    println(io, "  material_file = \"", s.material_file, "\"")
      print(io, "  sensitive     = ", s.sensitive)
end

function Base.show(io::IO, ::MIME"text/plain", f::FiberEntry)
    println(io, "FiberEntry \"", f.name, "\" (", f.kind, ")")
    println(io, "  mother         = \"", f.mother, "\"")
    isempty(f.reference) || println(io, "  reference      = \"", f.reference, "\"")
    println(io, "  start          = ", _fmt_tup(f.start), " mm")
    println(io, "  stop           = ", _fmt_tup(f.stop), " mm")
    if f.kind == "bent"
        println(io, "  bend_angle     = ", _fmt(f.bend_angle), " rad")
        println(io, "  bend_axis      = ", _fmt_tup(f.bend_axis))
    end
    println(io, "  material_file  = \"", f.material_file, "\"")
    if f.glued
        println(io, "  glued          = true")
        isempty(f.glue_file)    || println(io, "  glue_file      = \"", f.glue_file, "\"")
        isempty(f.glue_profile) || println(io, "  glue_profile   = \"", f.glue_profile, "\"")
    end
    isnan(f.start_reflectivity) || println(io, "  start_refl     = ", _fmt(f.start_reflectivity))
    isnan(f.end_reflectivity)   || println(io, "  end_refl       = ", _fmt(f.end_reflectivity))
    print(io, "  loop_id        = ", f.loop_id == -1 ? "(none)" : string(f.loop_id))
end

function Base.show(io::IO, ::MIME"text/plain", w::WrapEntry)
    println(io, "WrapEntry")
    println(io, "  scint         = \"", w.scint, "\"")
    isempty(w.g4name) || println(io, "  g4name        = \"", w.g4name, "\"")
    println(io, "  material_file = \"", w.material_file, "\"")
      print(io, "  cut           = ",
            isempty(w.cut) ? "auto (derived by PlaceManifest)" :
            string("[", join(("\"$c\"" for c in w.cut), ", "), "]"))
end

function Base.show(io::IO, ::MIME"text/plain", sp::SipmEntry)
    println(io, "SipmEntry \"", sp.name, "\"")
    println(io, "  ref_volume      = \"", sp.ref_volume, "\"")
    println(io, "  face_dir        = ", _fmt_tup(sp.face_dir))
    println(io, "  rel_pos         = ", _fmt_tup(sp.rel_pos), " mm")
    println(io, "  edge_length     = ", _fmt(sp.edge_length), " mm")
    println(io, "  fiber           = \"", sp.fiber, "\" (", sp.fiber_is_base ? "base" : "tip", ")")
    println(io, "  coupling_normal = ", _fmt_tup(sp.coupling_normal))
    println(io, "  coupling_pos    = ", _fmt_tup(sp.coupling_pos), " mm")
      print(io, "  coupling_width  = ", _fmt(sp.coupling_width), " mm")
end

# ---------------------------------------------------------------------------
#  CasingSpec
# ---------------------------------------------------------------------------

_casing_disabled(c::CasingSpec) =
    c.aluminum_thickness <= 0 && c.lead_thickness <= 0 && c.num_bars == 0

function Base.show(io::IO, c::CasingSpec)
    if _casing_disabled(c)
        print(io, "CasingSpec(disabled)")
    else
        print(io, "CasingSpec(num_bars=", c.num_bars,
              ", bar_width=", _fmt(c.bar_width),
              ", aluminum=", _fmt(c.aluminum_thickness),
              ", lead=", _fmt(c.lead_thickness), ")")
    end
end

function Base.show(io::IO, ::MIME"text/plain", c::CasingSpec)
    if _casing_disabled(c)
        print(io, "CasingSpec (disabled)")
        return
    end
    println(io, "CasingSpec")
    println(io, "  module_half_x      = ", _fmt(c.module_half_x), " mm")
    println(io, "  module_min_y       = ", _fmt(c.module_min_y), " mm")
    println(io, "  module_max_y       = ", _fmt(c.module_max_y), " mm")
    println(io, "  module_half_z      = ", _fmt(c.module_half_z), " mm")
    println(io, "  aluminum_thickness = ", _fmt(c.aluminum_thickness), " mm",
            c.aluminum_thickness <= 0 ? "  (disabled)" : "")
    println(io, "  lead_thickness     = ", _fmt(c.lead_thickness), " mm",
            c.lead_thickness     <= 0 ? "  (disabled)" : "")
    println(io, "  num_bars           = ", c.num_bars)
    println(io, "  bar_width          = ", _fmt(c.bar_width), " mm")
      print(io, "  scinti_z           = ", _fmt(c.scinti_z), " mm")
end

# ---------------------------------------------------------------------------
#  GeometryManifest
# ---------------------------------------------------------------------------

# Short tag used in the placement list of the REPL summary
_kind_tag(::ScintEntry) = "SCINT"
_kind_tag(::FiberEntry) = "FIBER"
_kind_tag(::WrapEntry)  = "WRAP "
_kind_tag(::SipmEntry)  = "SIPM "

# Identifier shown in the placement list
_entry_id(p::ScintEntry) = p.name
_entry_id(p::FiberEntry) = p.name
_entry_id(p::WrapEntry)  = "scint=" * p.scint
_entry_id(p::SipmEntry)  = p.name

# Compact one-line extra info for the placement list
_entry_extra(p::ScintEntry) = string("dims=", _fmt_tup(p.dims), p.sensitive ? "  sensitive" : "")
function _entry_extra(p::FiberEntry)
    s = "(" * p.kind * ")  mother=" * p.mother
    isempty(p.reference) || (s *= "  ref=" * p.reference)
    return s
end
_entry_extra(p::WrapEntry)  = isempty(p.cut) ? "cut=auto" : "cut=[" * join(p.cut, ", ") * "]"
_entry_extra(p::SipmEntry)  = "fiber=" * p.fiber

function Base.show(io::IO, m::GeometryManifest)
    print(io, "GeometryManifest(",
          isempty(m.setup_label) ? "" : "\"$(m.setup_label)\", ",
          length(scintillators(m)), " scints, ",
          length(fibers(m)),        " fibers, ",
          length(wraps(m)),         " wraps, ",
          length(sipms(m)),         " sipms)")
end

function Base.show(io::IO, ::MIME"text/plain", m::GeometryManifest)
    n_sc = length(scintillators(m))
    n_fb = length(fibers(m))
    n_wr = length(wraps(m))
    n_sp = length(sipms(m))
    println(io, "GeometryManifest",
            isempty(m.setup_label) ? "" : " \"$(m.setup_label)\"")
    println(io, "  ", n_sc, " scintillators, ", n_fb, " fibers, ",
            n_wr, " wraps, ", n_sp, " sipms")
    if isempty(m.placements)
        println(io, "  (no placements)")
    else
        # Cap the placement list so a 200-fibre manifest doesn't drown the REPL.
        n = length(m.placements)
        cap = 20
        println(io, "  placements (", n, "):")
        width = ndigits(n)
        show_n = min(n, cap)
        for i in 1:show_n
            p = m.placements[i]
            @printf(io, "    %*d. %s  %-20s %s\n", width, i, _kind_tag(p),
                    _entry_id(p), _entry_extra(p))
        end
        n > cap && println(io, "    ", " "^width, "  … and ", n - cap, " more")
    end
    print(io, "  casing: ", _casing_disabled(m.casing) ? "disabled" :
          "num_bars=$(m.casing.num_bars), bar_width=$(_fmt(m.casing.bar_width)) mm")
end

# ---------------------------------------------------------------------------
#  DetectorSpec (covers B1Spec/B2Spec/B3Spec via the abstract type)
# ---------------------------------------------------------------------------

function Base.show(io::IO, ::MIME"text/plain", spec::DetectorSpec)
    T = typeof(spec)
    println(io, nameof(T))
    names = fieldnames(T)
    isempty(names) && return
    pad = maximum(length ∘ string, names)
    n = length(names)
    for (i, name) in enumerate(names)
        v = getfield(spec, name)
        @printf(io, "  %-*s = ", pad, string(name))
        # Strings get quoted for clarity; everything else (numbers, unitful
        # quantities) prints with its own `show`.
        if v isa AbstractString
            print(io, "\"", v, "\"")
        else
            show(io, v)
        end
        i < n && println(io)
    end
end

# ---------------------------------------------------------------------------
#  ManifestBuilder
# ---------------------------------------------------------------------------

function Base.show(io::IO, b::ManifestBuilder)
    print(io, "ManifestBuilder(", length(b.placements), " placements)")
end

function Base.show(io::IO, ::MIME"text/plain", b::ManifestBuilder)
    n_sc = count(p -> p isa ScintEntry, b.placements)
    n_fb = count(p -> p isa FiberEntry, b.placements)
    n_wr = count(p -> p isa WrapEntry,  b.placements)
    n_sp = count(p -> p isa SipmEntry,  b.placements)
    println(io, "ManifestBuilder (", length(b.placements), " placements: ",
            n_sc, " scints, ", n_fb, " fibers, ", n_wr, " wraps, ", n_sp, " sipms)")
    println(io, "  next loop_id  = ", b.next_loop_id)
    println(io, "  casing        = ", b.casing)
    println(io, "  materials:")
    println(io, "    scint  = \"", b.materials.scint,  "\"")
    println(io, "    wrap   = \"", b.materials.wrap,   "\"")
    println(io, "    wls    = \"", b.materials.wls,    "\"")
      print(io, "    cement = \"", b.materials.cement, "\"")
end

# ---------------------------------------------------------------------------
#  ResolvedMaterials
# ---------------------------------------------------------------------------

function Base.show(io::IO, ::MIME"text/plain", r::ResolvedMaterials)
    println(io, "ResolvedMaterials")
    println(io, "  scint  = \"", r.scint,  "\"")
    println(io, "  wrap   = \"", r.wrap,   "\"")
    println(io, "  wls    = \"", r.wls,    "\"")
      print(io, "  cement = \"", r.cement, "\"")
end

# ---------------------------------------------------------------------------
#  Transform
# ---------------------------------------------------------------------------

function Base.show(io::IO, tf::Transform)
    print(io, "Transform(t=", _fmt_tup(Tuple(tf.translation)), ")")
end

function Base.show(io::IO, ::MIME"text/plain", tf::Transform)
    println(io, "Transform")
    println(io, "  translation = ", _fmt_tup(Tuple(tf.translation)), " mm")
      print(io, "  rotation    = ", tf.rotation)
end

# ---------------------------------------------------------------------------
#  Routing primitives
# ---------------------------------------------------------------------------

function Base.show(io::IO, p::G4DirectedPoint)
    print(io, "G4DirectedPoint(", p.point, ", direction=", _fmt_tup(Tuple(p.direction.vec)), ")")
end

function Base.show(io::IO, s::RouteSegment)
    if s.kind === STRAIGHT
        print(io, "RouteSegment(straight, ", _fmt_tup(Tuple(s.start)),
              " → ", _fmt_tup(Tuple(s.stop)), ")")
    else
        print(io, "RouteSegment(bent, ", _fmt_tup(Tuple(s.start)),
              " → ", _fmt_tup(Tuple(s.stop)),
              ", angle=", _fmt(s.bend_angle), " rad)")
    end
end

# ---------------------------------------------------------------------------
#  Geometry-check results
# ---------------------------------------------------------------------------

function Base.show(io::IO, c::FiberClash)
    print(io, "FiberClash(", c.a, " ↔ ", c.b, " in \"", c.mother,
          "\", d=", _fmt(c.distance), " mm)")
end

function Base.show(io::IO, o::ScintOverlap)
    print(io, "ScintOverlap(", o.a, " ↔ ", o.b, " in \"", o.mother, "\")")
end
