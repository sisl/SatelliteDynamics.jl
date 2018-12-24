using Documenter, SatelliteDynamics

makedocs(
    modules   = [SatelliteDynamics],
    doctest   = false,
    clean     = true,
    linkcheck = true,
    format    = Documenter.HTML(),
    sitename  = "SatelliteDynamics.jl",
    authors   = "Duncan Eddy",
    pages     = Any[
        "Home" => "index.md",
        "Modules" => Any[
            "modules/constants.md",
            "modules/universe.md",
            "modules/time.md",
            "modules/refsys.md",
            "modules/attitude.md",
            "modules/coordinates.md",
            "modules/astrodynamics.md",
            "Function Index" => "modules/function_index.md",
        ],
        # "Tutorials" => Any[
        #     "tutorials/universe_example.md",
        #     "tutorials/epoch_example.md",
        # ],
    ]
)

deploydocs(
    repo = "github.com/sisl/SatelliteDynamics.jl",
    devbranch = "master",
    devurl = "latest"
)