# SatelliteDynamics.jl Documentation

Welcome to the SatelliteDynamics.jl documentation! The package is meant to address
the needs of the satellite operator, academic researcher, and public enthusiast 
communities. In particular, it is designed and built to fill the following needs:

1. Open-source, MIT-lencesed software to remove barriers to entry of doing high-fidelity satellite modeling.
2. High-quality, validated, tested library for orbit and attitude dynamics modeling.
3. Easily acceible API design to make implementation of simulation and analysis intuitive.

## Getting Started: Installation And First Steps

To install the package, use the following command inside the Julia REPL:
```julia
Pkg.add("SatelliteDynamics")
```

To load the package, use the command:

```julia
using SatelliteDynamics
```

## Package Structure

The package is divided into a number of submodules each designed to provide a single, 
well-defined set of functions. The details of each specific module can be found 
on the relevant documentation page.

```@contents
Pages = [
"modules/constants.md",
"modules/universe.md",
"modules/time.md",
"modules/reference_systems.md",
"modules/attitude.md",
"modules/coordinates.md",
"modules/astrodynamics.md",
"modules/orbit_dynamics.md",
]
Depth = 2
```

There are also two subpackages which provide additional functionality in a 
self-contained manner. These are `EarthEnvrionment` and `Simulation`.

**EarthEnvrionment:**
```@contents
Pages = [
"modules/earth_environment/space_weather.md",
"modules/earth_environment/nrlmsise00.md",
]
Depth = 3
```

**Simulation:**
```@contents
Pages = [
"modules/simulation/integrators.md",
"modules/simulation/propagators.md",
]
Depth = 3
```

## Examples

The best way to learn how to use any software is to actually see it in action and
use it for yourself. A number of tutorials showing how to use the modules provided
as part of SatelliteDynamics.jl in your software are included as part of this
documentation. You can find the list of tutorials below:

```@contents
Pages = [
    "tutorials/orbit_propagation_example.md",
]
Depth = 2
```

## Citing

The software in this package was developed as part of academic research.
If you would like to help support it, please star the repository as such metrics on
usage help guage interest, secure funding, and continue development. 

If you use the software as part of your research, teaching, or other activities, we would be grateful if you could cite our work.

## Contributing

Contributions are always welcome! For feature requests, questions, or if a bug is found please create a github issue.