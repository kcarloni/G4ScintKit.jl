module G4ScintKit

# All external dependencies are imported here; the included files below assume
# these names are in module scope.
using AstroParticleUnits
using Dubins
using HDF5
using LinearAlgebra
using Printf
using Rotations
using StaticArrays
using StructArrayTables

# --- output readers ---
# Public surface: `load(outdir)` reads the C++ HDF5 output into a NamedTuple
# of StructArrays, attaching units and concatenating across runs/files.
# `h5read(StructArray, file, group)` is the lower-level ad-hoc reader for
# arbitrary groups; `read_sipm_voltage_trace` is a legacy typed reader kept
# for compatibility.
include("output/hdf5.jl")
export h5display, h5read
include("output/columns.jl")
include("output/load.jl")
export load, available_groups, ALL_GROUPS, SimulationOutput
include("output/read.jl")
export read_sipm_voltage_trace

# --- geometry: low-level manifest representation ---
include("geometry/manifest.jl")
export GeometryManifest, PlacementEntry
export ScintEntry, FiberEntry, WrapEntry, SipmEntry, CasingSpec
export scintillators, fibers, wraps, sipms
export fiber_length, fiber_lengths
export write_manifest

# --- geometry: position / displacement / direction bound to a reference frame ---
include("geometry/g4coordinate.jl")
export G4Coordinate, G4Vector, G4Direction, midpoint

# --- geometry: reference-frame algebra (Transform, world-frame resolution) ---
include("geometry/frames.jl")
export Transform

# --- geometry: reusable geometric primitives (OBBs, fillets, tile extent) ---
include("geometry/geometry_helpers.jl")
export extent, face_center, center_pos

# --- geometry: material-file resolution ---
include("geometry/materials.jl")
export ResolvedMaterials, fiber_cross_section, default_goddess_root

# --- geometry: high-level DetectorSpec interface (abstract base + build_manifest) ---
include("geometry/specs/detector_spec.jl")
export DetectorSpec, LengthQ, EnergyQ, TimeQ, VoltageQ, build_manifest, strip_units

# --- geometry: g4sipm model registry (depends on LengthQ) ---
include("geometry/sipms.jl")
export SipmModelInfo, SIPM_MODELS, sipm_model_info, sipm_edge_length

# --- particle-list CSV writer (mirrors C++ ParticleListSource) ---
# Included before g4run.jl so run_simulation can dispatch on a Vector of entries.
include("run/particle_list.jl")
export ParticleListEntry, write_particle_list

# --- bash-pipeline run wrappers (depend on GeometryManifest, DetectorSpec) ---
include("run/g4run.jl")
export run_simulation, run_visu

# --- geometry: pre-flight validation ---
include("geometry/geometry_check.jl")
export FiberClash, fiber_clashes, ScintOverlap, scint_overlaps, check_geometry

# --- geometry builders: fibre routing vocabulary + waypoint router ---
include("geometry/builders/fiber_routing.jl")
export G4DirectedPoint, Waypoint, RouteSegment, FiberKind, STRAIGHT, BENT
export fiber_entries, route_path

# --- geometry builders: Dubins shortest-path router ---
include("geometry/builders/dubins.jl")
export route_fiber

# --- geometry builders: ManifestBuilder ---
include("geometry/builders/manifest_builder.jl")
export ManifestBuilder
export add_scint!, add_fiber_straight!, add_fiber_bent!, add_fiber_routed!
export add_fiber_path!, add_wrapping!, add_sipm!
export new_loop_id!, casing_from_extent, add_casing!, to_manifest

# --- geometry builders: composite assembly helpers ---
include("geometry/builders/detector_assembly.jl")
export bar_lattice_spacing, add_scint_row!
export add_inline_sipm!, add_bundle_sipm!, bundle_fiber_endpoints

# --- custom REPL pretty-printers for the user-facing structs ---
# (must be last: every type it specialises on is defined above.)
include("display.jl")


end # module G4ScintKit