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

# --- output readers (inactive) ---
# include("output/hdf5.jl")
# export h5display, h5read
# include("output/read.jl")
# export read_sipm_voltage_trace

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
export ResolvedMaterials

# --- bash-pipeline run wrappers (depend on GeometryManifest) ---
include("run/g4run.jl")
export run_simulation, run_visu

# --- geometry: high-level DetectorSpec interface (abstract base + build_manifest) ---
include("geometry/specs/detector_spec.jl")
export DetectorSpec, build_manifest

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

# --- geometry builders: Layer-1 ManifestBuilder ---
include("geometry/builders/manifest_builder.jl")
export ManifestBuilder
export add_scint!, add_fiber_straight!, add_fiber_bent!, add_fiber_routed!
export add_fiber_path!, add_wrapping!, add_sipm!
export new_loop_id!, casing_from_extent, add_casing!, to_manifest

# --- geometry builders: Layer-2 composite assembly helpers ---
include("geometry/builders/detector_assembly.jl")
export bar_lattice_spacing, add_scint_row!
export add_inline_sipm!, add_bundle_sipm!, bundle_fiber_endpoints

# --- custom REPL pretty-printers for the user-facing structs ---
# (must be last: every type it specialises on is defined above.)
include("display.jl")


end # module G4ScintKit