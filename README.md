# SatelliteDynamics.jl [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://sisl.github.io/SatelliteDynamics.jl/stable/) [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://sisl.github.io/SatelliteDynamics.jl/dev/) [![Build Status](https://github.com/sisl/SatelliteDynamics.jl/actions/workflows/CI.yaml/badge.svg?branch=main)](https://github.com/sisl/SatelliteDynamics.jl/actions/workflows/CI.yml?query=branch%3Amain) [![Coverage](https://codecov.io/gh/sisl/SatelliteDynamics.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/sisl/SatelliteDynamics.jl)

> [!TIP]
> This package is no longer under active development, but still generally works.
>
> Please consider using [brahe](https://github.com/duncaneddy/brahe) as an alternative for new projects. It is actively maintained and has a broader feature set.

SatelliteDynamics package is meant to address the needs of the satellite operator, academic researcher, and public enthusiast communities.

Current high-fidelity satellite dynamics modeling software generally falls into two camps:
1. It is commercially licensed and closed-source code, making it difficult, if not impossible, to be used by hobbyists and academic researchers.
2. There is a steep learning curve which makes correct use of the underlying libraries difficult.

These two challenges make it an unfortunately common occurance that guidance, navigation, and control engineers will frequently reimplement common astrodynamics libraries for each new project.

With these two deficienties in mind, SatelliteDynamics.jl aims to provide a open-source, MIT-licensed, high fidelity astrodynamics toolbox to help make it easy to perform high-quality simulation and analysis of satellite attitude and orbit dynamics.

## Getting Started: Installation And First Steps

To install the package, use the following command inside the Julia REPL:
```julia
Pkg.add("SatelliteDynamics")
```

To load the package, use the command:

```julia
using SatelliteDynamics
```

## Documentation

The documentation for the package can be found here: <https://sisl.github.io/SatelliteDynamics.jl/latest>

More example code and more thorough documentation will be added as time permits.

## Contributing

Contributions are welcome! 

### Local Development

**Running Unit Tests**

The package has a suite of unit tests that can be run using the following command:

```
julia --project -e 'using Pkg; Pkg.test("SatelliteDynamics")'
```

Or alternatively, you can run the tests from the package REPL:

```zsh
julia --project
```

```julia
]; test
```

**Building the Documentation**

To build the documentation locally, you can use the following command:

```
cd docs
```

Then

```
julia --project make.jl
```

The documentation will be built in the `docs/build` directory and can be viewed by opening the `index.html` file in a web browser.
