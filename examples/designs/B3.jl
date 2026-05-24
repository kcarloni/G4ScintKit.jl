
"""
    B3Spec(; kwargs...)

The "B3" design: a single scintillator tile with a reflective wrapping,
read out by two WLS fibers that meet at a shared SiPM towards the front.

Bundled with G4ScintKit as a canonical example exercising:
multi-fiber routing via `add_fiber_path!`, `bundle_fiber_endpoints`,
and `add_bundle_sipm!`. Lengths are unitful.
"""
Base.@kwdef struct B3Spec <: DetectorSpec

    # scintillator tile (width, thickness, length)
    scint_width::LengthQ     = 4.95u"cm"
    scint_thickness::LengthQ = 0.95u"cm"
    scint_length::LengthQ    = 1.875u"m"

    # outer casing (<= 0 disables the box / sheet)
    aluminum_thickness::LengthQ = 0.0u"mm"
    lead_thickness::LengthQ     = 0.0u"mm"

    # material abbreviations (resolved against goddess_root)
    goddess_root::String    = default_goddess_root()
    scint_material::String  = "fermilab"
    wrap_material::String   = "tio2"
    wls_material::String    = "y11-300-r1"
    cement_material::String = "air1mm"
end

function G4ScintKit.build_manifest(spec::B3Spec)

    b = ManifestBuilder(spec)

    # ---------------------------------------
    # set up placement + property constants

    fiber_diameter = fiber_cross_section(b.materials.wls).width
    fiber_min_radius = 5.0u"cm"
    fiber_length_buffer = 1.0u"mm"
    fiber_offcenter_shift = spec.scint_width / 4

    sipm_coupling_width = 0.25u"mm"
    sipm_edge_length    = 2 * fiber_diameter + 1u"mm"

    # ---------------------------------------

    origin = G4Coordinate(0.0, 0.0, 0.0, "world")

    # scintillator first...
    scint_dims = ( spec.scint_width, spec.scint_thickness, spec.scint_length )
    scint_1 = add_scint!(b;
        name = "scint_1",
        dims = scint_dims,
        pos = origin,
        sensitive = true
    )

    # next we decide where we will put the sipm...
    scint_far_end = face_center(scint_1, :z, :min)
    scint_near_end = face_center(scint_1, :z, :max)
    sipm_pos = scint_near_end + G4Vector( 0u"cm", 0u"cm", 10u"cm", "world" )
    sipm_face_dir = G4Direction(0, 0, -1, "world")
    coupling_face_pos = sipm_pos + sipm_coupling_width * sipm_face_dir

    # extend along z in both directions to poke out of the scintillator bar
    fiber_z_offset = G4Vector( 0u"cm", 0u"cm", fiber_length_buffer, "world" )

    # shift along x to fit two fibers inside the bar
    fiber_x_offset = G4Vector( fiber_offcenter_shift, 0u"cm", 0u"cm", "world" )

    # arrange our fibers at the sipm face
    endpoints = bundle_fiber_endpoints(;
        num_fibers = 2,
        pitch = 1.1 * fiber_diameter,
        plane_center = coupling_face_pos,
        plane_normal = sipm_face_dir,
        axis = G4Direction( 1, 0, 0, "world" )
    )
    fiber_end_backstep = 2u"cm" * sipm_face_dir

    # head towards the sipm
    fiber_heading_dir = G4Direction(0, 0, 1, "world")
    fiber_towards_step = 2u"cm" * fiber_heading_dir

    fiber_1_waypoints = [
        scint_far_end - fiber_z_offset + fiber_x_offset,
        scint_near_end + fiber_z_offset + fiber_x_offset,
        G4DirectedPoint(
            scint_near_end + fiber_z_offset + fiber_x_offset + fiber_towards_step,
            fiber_heading_dir,
        ),
        G4DirectedPoint(
            endpoints[1] + fiber_end_backstep,
            G4Direction( -1 * sipm_face_dir )
        ),
        endpoints[1]
    ]
    fiber_2_waypoints = [
        scint_far_end - fiber_z_offset - fiber_x_offset,
        scint_near_end + fiber_z_offset - fiber_x_offset,
        G4DirectedPoint(
            scint_near_end + fiber_z_offset - fiber_x_offset + fiber_towards_step,
            fiber_heading_dir,
        ),
        G4DirectedPoint(
            endpoints[2] + fiber_end_backstep,
            -1 * sipm_face_dir
        ),
        endpoints[2]
    ]

    fiber_1 = add_fiber_path!( b;
        name = "fiber_1",
        waypoints = fiber_1_waypoints,
        min_radius = fiber_min_radius,
        start_reflectivity = 1.0,
    )
    add_fiber_path!( b;
        name = "fiber_2",
        waypoints = fiber_2_waypoints,
        min_radius = fiber_min_radius,
        start_reflectivity = 1.0,
    )

    # wrapping has to come *after* the fiber, fiber pokes out
    add_wrapping!(b; scint=scint_1 )

    add_bundle_sipm!( b;
        name = "sipm",
        anchor_fiber = last( fiber_1 ),
        endpoints = endpoints,
        face_dir = sipm_face_dir,
        edge_length = sipm_edge_length,
        coupling_width = sipm_coupling_width
    )

    add_casing!(b;
        aluminum_thickness = spec.aluminum_thickness,
        lead_thickness     = spec.lead_thickness
    )

    return to_manifest(b; setup_label = "B3")
end
