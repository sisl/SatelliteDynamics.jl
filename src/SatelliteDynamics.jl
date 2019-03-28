__precompile__(true)
module SatelliteDynamics

# Usings
using Reexport

# Includes
include("constants.jl")
include("universe.jl")
include("time.jl")
include("reference_systems.jl")
include("attitude.jl")
include("coordinates.jl")
include("astrodynamics.jl")
include("orbit_dynamics.jl")

# Export Values
@reexport using SatelliteDynamics.Constants
@reexport using SatelliteDynamics.Universe
@reexport using SatelliteDynamics.Time
@reexport using SatelliteDynamics.ReferenceSystems
@reexport using SatelliteDynamics.Attitude
@reexport using SatelliteDynamics.Coordinates
@reexport using SatelliteDynamics.Astrodynamics
@reexport using SatelliteDynamics.OrbitDynamics

# Export EarthEnvironment submodule
include(joinpath(".", "earth_environment", "earth_environment.jl"))
include(joinpath(".", "simulation", "simulation.jl"))

end # module
