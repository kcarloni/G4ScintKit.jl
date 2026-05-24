# G4ScintKit.jl smoke tests. Detector-spec-specific tests live in the
# downstream G4Tambo.jl package (manifest reproducibility, B1/B2/... knob
# coverage). Here we only verify the kit's generic primitives.

using G4ScintKit
using Test
using LinearAlgebra: norm, normalize
using StaticArrays: SVector, SMatrix
using AstroParticleUnits

using G4ScintKit: _fiber_polyline
using G4ScintKit: OBB, point_in_obb, segment_obb_interval, fillet_corner

# A minimal DetectorSpec just for exercising ManifestBuilder / ResolvedMaterials
# without depending on any concrete detector design.
Base.@kwdef struct _SmokeSpec <: DetectorSpec
    goddess_root::String   = G4ScintKit._default_goddess_root()
    scint_material::String = "fermilab"
    wrap_material::String  = "tio2"
    wls_material::String   = "y11-300-r1"
    cement_material::String = "air1mm"
end

@testset "G4ScintKit" begin

@testset "module loads / key exports present" begin
    for name in (:GeometryManifest, :ManifestBuilder, :DetectorSpec,
                 :build_manifest, :write_manifest, :route_fiber, :route_path,
                 :G4Coordinate, :G4Vector, :G4Direction, :Transform,
                 :fiber_clashes, :scint_overlaps, :check_geometry)
        @test isdefined(G4ScintKit, name)
    end
end

@testset "ManifestBuilder — minimal hand-built detector" begin
    b = ManifestBuilder(_SmokeSpec())
    s  = add_scint!(b; name="s", dims=(50.0, 10.0, 1000.0))
    f  = add_fiber_straight!(b; name="f", mother="s",
            start=G4Coordinate((0.0, 0.0, 500.0),  "s"),
            stop =G4Coordinate((0.0, 0.0,-500.0),  "s"),
            glued=true, end_reflectivity=1.0)
    w  = add_wrapping!(b; scint="s")
    sp = add_sipm!(b; name="sp", fiber="f", face_dir=(0.0, 0.0,-1.0),
            rel_pos=G4Coordinate((0.0,0.0,520.0), "world"), edge_length=2.0,
            coupling_normal=(0.0,0.0,-1.0),
            coupling_pos=G4Coordinate((0.0,0.0,0.0), "f"),
            coupling_width=0.25)

    m = to_manifest(b; setup_label="smoke")
    @test m isa GeometryManifest
    @test m.setup_label == "smoke"
    @test [typeof(p) for p in m.placements] ==
          [ScintEntry, FiberEntry, WrapEntry, SipmEntry]
end

@testset "fiber_router — collinear, pure left turn" begin
    zn = (0.0, 0.0, 1.0)

    segs = route_fiber(G4DirectedPoint((0.,0,0),(1.,0,0)),
                       G4DirectedPoint((100.,0,0),(1.,0,0));
                       min_radius=30.0, plane_normal=zn)
    @test length(segs) == 1
    @test segs[1].kind === STRAIGHT
    @test norm(segs[1].stop - segs[1].start) ≈ 100.0

    # (0,0)+x -> (50,50)+y : 90° left arc of radius 50
    segs = route_fiber(G4DirectedPoint((0.,0,0),(1.,0,0)),
                       G4DirectedPoint((50.,50,0),(0.,1,0));
                       min_radius=50.0, plane_normal=zn)
    @test length(segs) == 1 && segs[1].kind === BENT
    @test segs[1].bend_angle ≈ π/2
    @test segs[1].bend_axis ≈ SVector(0.0, 0.0, 1.0)
end

@testset "fiber_clashes / scint_overlaps detect bad geometry" begin
    crossing = GeometryManifest(placements=PlacementEntry[
        FiberEntry(name="A", kind="straight", mother="world",
                   start=(0.0,0.0,-100.0), stop=(0.0,0.0,100.0),
                   material_file="x", loop_id=-1),
        FiberEntry(name="B", kind="straight", mother="world",
                   start=(-100.0,0.0,0.0), stop=(100.0,0.0,0.0),
                   material_file="x", loop_id=-1)])
    @test length(fiber_clashes(crossing)) == 1

    overlap = GeometryManifest(placements=PlacementEntry[
        ScintEntry(name="s0", dims=(50.0,10.0,1000.0), pos=(0.0,0.0,0.0),  material_file="x"),
        ScintEntry(name="s1", dims=(50.0,10.0,1000.0), pos=(30.0,0.0,0.0), material_file="x")])
    @test length(scint_overlaps(overlap)) == 1

    @test_throws ErrorException check_geometry(crossing)
    @test check_geometry(crossing; on_clash=:silent) === crossing
end

end  # @testset "G4ScintKit"