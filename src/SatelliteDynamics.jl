__precompile__(true)
module SatelliteDynamics

# Usings
using Dates
using Printf
using LinearAlgebra
using SOFA

using StaticArrays: SVector, SMatrix

# Define constants for array indexing
const idx1t3 = SVector(1, 2, 3)
const idx4t6 = SVector(4, 5, 6)
const idx1t6 = SVector(1, 2, 3, 4, 5, 6)

# Includes
include("constants.jl")
include("universe.jl")
include("time.jl")
include("reference_systems.jl")
include("attitude.jl")
include("coordinates.jl")
include("astrodynamics.jl")
include("orbit_dynamics.jl")
include("sgp_models.jl")

# Export EarthEnvironment submodule
include(joinpath(".", "earth_environment", "earth_environment.jl"))
include(joinpath(".", "simulation", "simulation.jl"))

end # module
