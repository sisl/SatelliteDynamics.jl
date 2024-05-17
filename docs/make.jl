using Documenter
using SatelliteDynamics

include("src/makeplots.jl")

DocMeta.setdocmeta!(SatelliteDynamics, :DocTestSetup, :(using SatelliteDynamics); recursive=true)

# Generate plots
makeplots()

# Generate documents
makedocs(
    modules   = [SatelliteDynamics],
    doctest   = false,
    clean     = true,
    linkcheck = false,
    checkdocs = :none,
    format    = Documenter.HTML(
        canonical="https://sisl.github.io/SatelliteDynamics.jl",
        edit_link="main",
        assets=String[],
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
    sitename  = "SatelliteDynamics.jl",
    authors   = "Duncan Eddy",
    pages     = Any[
        "Home" => "index.md",
        "Modules" => Any[
            "modules/constants.md",
            "modules/universe.md",
            "modules/time.md",
            "Reference Systems" => "modules/reference_systems.md",
            "modules/attitude.md",
            "modules/coordinates.md",
            "modules/astrodynamics.md",
            "modules/sgp_models.md",
            "Orbit Dynamics" => "modules/orbit_dynamics.md",
            "Earth Environment" => Any[
                "Space Weather" => "modules/earth_environment/space_weather.md",
                "NRLMSISE00" => "modules/earth_environment/nrlmsise00.md"
            ],
            "Simulation" => Any[
                "modules/simulation/integrators.md",
                "modules/simulation/propagators.md"
            ]
        ],
        "Tutorials" => Any[
            "tutorials/orbit_propagation_example.md",
        ],
        "Function Index" => "function_index.md",
    ]
)

deploydocs(
    repo = "github.com/sisl/SatelliteDynamics.jl",
    devbranch = "main",
)