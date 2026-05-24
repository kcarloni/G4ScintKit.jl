
"""
    B4Spec(; kwargs...)

The "B4" design: four scintillator bars with a reflective wrapping. Each bar is
read out by two WLS fibers; fibers loop at the back in an AABB AABB pattern
and meet at the front at a shared SiPM.

Bundled with G4ScintKit as a canonical example exercising:
multi-bar `add_scint_row!`, looped routing, and a grid-arrangement SiPM bundle.
Lengths are unitful.
"""
Base.@kwdef struct B4Spec <: DetectorSpec

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

function G4ScintKit.build_manifest(spec::B4Spec)

    b = ManifestBuilder(spec)

    # ---------------------------------------
    # set up placement + property constants

    fiber_diameter = fiber_cross_section(b.materials.wls).width
    fiber_min_radius = 5.0u"cm"
    fiber_length_buffer = 1.0u"mm"
    fiber_offcenter_shift = spec.scint_width / 4

    num_bars = 4
    num_fibers = 4

    sipm_coupling_width = 0.25u"mm"
    sipm_edge_length    = 4 * fiber_diameter + 1u"mm"

    # ---------------------------------------

    origin = G4Coordinate(0.0, 0.0, 0.0, "world")

    # add scintillator bars first...
    scint_dims = ( spec.scint_width, spec.scint_thickness, spec.scint_length )
    inter_bar_spacing = 1u"mm"
    scints = add_scint_row!(b;
        num_bars, inter_bar_spacing,
        dims = scint_dims,
        center = origin, axis = :x
    )

    # next we decide where we will put the sipm...
    scint_near_end = midpoint(face_center(scints[1],  :z, :max), face_center(scints[end], :z, :max))
    sipm_pos = scint_near_end + G4Vector( 0u"cm", 0u"cm", 20u"cm", "world" )
    sipm_face_dir = G4Direction(0, 0, -1, "world")
    coupling_face_pos = sipm_pos + sipm_coupling_width * sipm_face_dir

    # extend along z in both directions to poke out of the scintillator bar
    fiber_z_offset = fiber_length_buffer * G4Direction( 0, 0, 1, "world" )

    # shift along x to fit two fibers inside the bar
    fiber_x_offset = fiber_offcenter_shift * G4Direction( 1, 0, 0, "world" )

    # arrange our fibers at the sipm face
    endpoints = bundle_fiber_endpoints(;
        num_fibers = num_fibers * 2,            # both ends
        pitch = 1.1 * fiber_diameter,
        plane_center = coupling_face_pos,
        plane_normal = sipm_face_dir,
        axis = G4Direction( 0, 1, 0, "world" ),
        arrangement = :grid,
        rows = 2, cols = 4
    )
    sort!(endpoints; by = e -> (e.x, e.y))

    # route from grid bundle to row
    fiber_end_backstep = 1u"cm" * sipm_face_dir
    y_dir = G4Direction( 0, 1, 0, "world" )
    fiber_heading_dir = G4Direction(0, 0, 1, "world")
    fiber_towards_step = 1u"cm" * fiber_heading_dir

    # place all the fibers in the scintillator
    scint_near_positions = face_center.( scints, :z, :max )
    scint_far_positions = face_center.( scints, :z, :min )

    fiber_near_positions = G4Coordinate[]
    fiber_far_positions = G4Coordinate[]
    for i in 1:8
        j = (i+1)÷2                              # scintillator index
        sx = 1 - 2 * (i % 2)                     # x-shift sign: -1 or +1
        x_offset = sx * fiber_x_offset
        y_offset = endpoints[i].y * y_dir

        push!( fiber_near_positions, scint_near_positions[j] + fiber_z_offset + x_offset + y_offset )
        push!( fiber_far_positions,  scint_far_positions[j]  - fiber_z_offset + x_offset + y_offset )
    end

    fibers = []
    for i in 1:4
        v = Waypoint[]

        # head back from the sipm
        ep = endpoints[i]
        push!(v, ep, G4DirectedPoint(
            ep + fiber_end_backstep,
            sipm_face_dir
        ))

        # head towards + through the scintillator bar
        for pt in (fiber_near_positions[i] + fiber_towards_step, fiber_near_positions[i], fiber_far_positions[i])
            push!( v, G4DirectedPoint( pt, -1 * fiber_towards_step ))
        end

        # u-turn and head back
        j = i + num_fibers
        for pt in ( fiber_far_positions[j], fiber_near_positions[j], fiber_near_positions[j] + fiber_towards_step )
            push!( v, G4DirectedPoint( pt,  1 * fiber_towards_step ))
        end

        # back towards the sipm
        ep = endpoints[j]
        push!(v,
            G4DirectedPoint(
                ep + fiber_end_backstep,
                -1 * sipm_face_dir
            ), ep
        )

        fiber = add_fiber_path!( b;
            name = "fiber_$i",
            waypoints = v,
            min_radius = fiber_min_radius,
        )
        push!( fibers, fiber )
    end

    # wrapping has to come *after* the fiber, fiber pokes out
    for scint in scints
        add_wrapping!(b; scint=scint )
    end

    add_bundle_sipm!( b;
        name = "sipm",
        anchor_fiber = last( fibers[1] ),
        endpoints = endpoints,
        face_dir = sipm_face_dir,
        edge_length = sipm_edge_length,
        coupling_width = sipm_coupling_width
    )

    add_casing!(b;
        aluminum_thickness = spec.aluminum_thickness,
        lead_thickness     = spec.lead_thickness
    )

    return to_manifest(b; setup_label = "B4" )
end
