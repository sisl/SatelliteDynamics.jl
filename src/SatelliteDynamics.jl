module SatelliteDynamics

# Usings
using Reexport
using LinearAlgebra

# Includes
include("constants.jl")
include("universe.jl")
# include("time.jl")
# include("refsys.jl")
# include("coordinates.jl")
# include("astrodynamics.jl")
# include("orbit_dynamics.jl")
# include("simulation.jl")
# include("data_structures.jl")

# Export Values
@reexport using SatelliteDynamics.Constants
@reexport using SatelliteDynamics.Universe
# @reexport using SatelliteDynamics.Time
# @reexport using SatelliteDynamics.RefSys
# @reexport using SatelliteDynamics.Coordinates
# @reexport using SatelliteDynamics.Astrodynamics
# @reexport using SatelliteDynamics.OrbitDynamics
# @reexport using SatelliteDynamics.Simulation
# @reexport using SatelliteDynamics.DataStructures

# Update package data files at build time
update_eop(:C04_14)
update_eop(:C04_80)
update_eop(:FINALS_2000)

end # module
