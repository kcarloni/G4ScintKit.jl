# Open the B3 design in the Geant4 visualizer.
#
# Usage (from anywhere):
#   julia G4ScintKit.jl/examples/run_b3_visu.jl
#
# Requires the C++ binary to be built (bash bash_scripts/2_compile.sh) and
# bash_scripts/setup_paths.sh.example to find your Geant4 install (either
# edit it or export GEANT4_INSTALL_DIR).

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using AstroParticleUnits
using G4ScintKit

# include the B3 design spec = a single bar with two fibers, leading to one sipm
include(joinpath(@__DIR__, "designs", "B3.jl"))

spec = B3Spec()

# from a spec, we build a manifest = detailed description of the detector geometry, for Geant
manifest = build_manifest( spec )

# we can compute properties of our design, e.g.
fiber_lengths( manifest )

# and we can run Geant4 in visual mode to inspect it: 
run_visu( manifest )

# we can also run Geant4 in batch mode: 
outdir = joinpath(@__DIR__, "out", "B3")
run_simulation(B3Spec();
    outdir = outdir,
    nevents = 10,
    injparticle = "mu-",
    injenergy = "3_GeV",         # underscore = space; run.sh substitutes
    trackphotons = false,
)

# and inspect the output!
out = load( outdir; )

out.input
out.particle_hits
out.optical_photons
out.sipm_hits