__precompile__(true)
module EarthEnvironment

using Reexport

# EarthEnvironment submodules
include("space_weather.jl")
include("nrlmsise00.jl")

@reexport using SatelliteDynamics.EarthEnvironment.SpaceWeather
@reexport using SatelliteDynamics.EarthEnvironment.NRLMSISE00

end # EarthEnvironment Module