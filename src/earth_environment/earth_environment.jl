__precompile__(true)
module EarthEnvironment

using Reexport

# EarthEnvironment submodules
include("nrlmsise00.jl")

@reexport using SatelliteDynamics.EarthEnvironment.NRLMSISE00

end # EarthEnvironment Module