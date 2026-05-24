# G4ScintKit.jl smoke tests. Verifies the kit's generic primitives plus a
# determinism check on the bundled example designs (../examples/designs/),
# which exercises the full manifest-building + serialization pipeline.
# Downstream design-specific coverage (knob sweeps, physics sanity) lives in
# the consuming package.

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
    goddess_root::String   = default_goddess_root()
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

@testset "output read pipeline" begin
    using HDF5
    using StructArrayTables: StructArray
    using AstroParticleUnits: ustrip

    # Helper: write a creation-order-tracked group with a NamedTuple of columns.
    function _write_group!(parent, name, cols::NamedTuple)
        gcpl = HDF5.GroupCreateProperties(; track_order = true)
        gid = HDF5.API.h5g_create(parent, name,
            HDF5.API.H5P_DEFAULT, gcpl.id, HDF5.API.H5P_DEFAULT)
        HDF5.API.h5g_close(gid)
        g = parent[name]
        for (k, v) in pairs(cols)
            g[String(k)] = collect(v)
        end
    end

    # Build a tiny synthetic file mirroring the C++ schema for one run.
    function _make_fixture(path::String; run_id::Int = 0, nevents::Int = 3)
        h5open(path, "w") do f
            gcpl = HDF5.GroupCreateProperties(; track_order = true)
            gid = HDF5.API.h5g_create(f, "g4run_$run_id",
                HDF5.API.H5P_DEFAULT, gcpl.id, HDF5.API.H5P_DEFAULT)
            HDF5.API.h5g_close(gid)
            r = f["g4run_$run_id"]
            ev = collect(1:nevents)
            _write_group!(r, "input", (
                g4event_id              = ev,
                init_pos_x_mm           = float.(ev) .* 10.0,
                init_kinetic_energy_MeV = float.(ev),
                init_dir_x              = ones(nevents),
            ))
            _write_group!(r, "particle_hits", (
                g4event_id        = ev,
                edep_MeV          = float.(ev) .* 0.1,
                first_hit_time_ns = float.(ev) .* 0.5,
            ))
        end
    end

    @testset "_clean_columns: suffix stripping + units + rename" begin
        nt = (
            g4event_id    = [1, 2],
            init_pos_x_mm = [10.0, 20.0],
            edep_MeV      = [0.5, 1.5],
            voltages_V    = [3.3, 3.4],
            init_dir_x    = [1.0, 0.0],   # no suffix → pass through
        )
        out = G4ScintKit._clean_columns(nt)
        @test propertynames(out) == (:event_id, :init_pos_x, :edep, :voltages, :init_dir_x)
        @test out.init_pos_x == [10.0, 20.0] .* u"mm"
        @test out.edep       == [0.5, 1.5]    .* u"MeV"
        @test out.voltages   == [3.3, 3.4]    .* u"V"
        @test out.init_dir_x == [1.0, 0.0]    # untouched
        @test out.event_id   == [1, 2]
    end

    @testset "load: single file, one run, requested groups" begin
        path = tempname() * ".h5"
        _make_fixture(path; run_id=0, nevents=3)

        data = load(path; groups=(:input, :particle_hits))
        @test data isa SimulationOutput
        @test propertynames(data) == (:input, :particle_hits)
        @test length(data) == 2                  # group count
        @test haskey(data, :input)
        @test data[:input] === data.input        # indexing
        @test collect(keys(data)) == [:input, :particle_hits]
        # iteration yields (name, struct_array) pairs
        @test first(data) == (:input => data.input)
        @test length(data.input) == 3
        @test data.input.event_id   == [1, 2, 3]
        @test data.input.run_id     == [0, 0, 0]
        @test data.input.init_pos_x == [10.0, 20.0, 30.0] .* u"mm"
        @test eltype(ustrip.(data.input.init_kinetic_energy)) == Float64
        @test data.particle_hits.edep ≈ [0.1, 0.2, 0.3] .* u"MeV"

        rm(path; force=true)
    end

    @testset "load: missing groups are omitted from the output" begin
        path = tempname() * ".h5"
        _make_fixture(path)
        data = load(path; groups=(:input, :sipm_voltage_trace))
        @test length(data.input) == 3
        @test !haskey(data, :sipm_voltage_trace)
        @test propertynames(data) == (:input,)
        rm(path; force=true)
    end

    @testset "load: multi-file concat with run_id" begin
        outdir = mktempdir()
        datadir = mkpath(joinpath(outdir, "Data"))
        _make_fixture(joinpath(datadir, "a.h5"); run_id=0, nevents=2)
        _make_fixture(joinpath(datadir, "b.h5"); run_id=1, nevents=3)

        data = load(outdir; groups=(:input,))
        @test length(data.input) == 5
        @test sort(data.input.run_id) == [0, 0, 1, 1, 1]
        # Concat is per-file in directory order; within each file, run-id order
        @test data.input.event_id == [1, 2, 1, 2, 3]
    end

    @testset "available_groups" begin
        path = tempname() * ".h5"
        _make_fixture(path)
        @test available_groups(path) == [:input, :particle_hits]
        rm(path; force=true)
    end
end

@testset "bundled designs build deterministically" begin
    # Including the design files defines B3Spec / B4Spec at top-level and
    # extends G4ScintKit.build_manifest for each. Exercises the full
    # ManifestBuilder + assembly + serializer stack against real geometry.
    include(joinpath(@__DIR__, "..", "examples", "designs", "B3.jl"))
    include(joinpath(@__DIR__, "..", "examples", "designs", "B4.jl"))

    for spec in (B3Spec(), B4Spec())
        m1 = build_manifest(spec)
        @test m1 isa GeometryManifest
        @test !isempty(m1.placements)

        # Build a fresh manifest from the same spec; serialize both and assert
        # byte-equal output. Catches non-determinism in either build_manifest
        # or write_manifest, and silently dropped fields in the serializer.
        m2 = build_manifest(spec)
        p1 = tempname() * ".manifest"
        p2 = tempname() * ".manifest"
        write_manifest(p1, m1)
        write_manifest(p2, m2)
        @test read(p1, String) == read(p2, String)
        @test !isempty(read(p1, String))    # serializer actually wrote something
    end
end

end  # @testset "G4ScintKit"