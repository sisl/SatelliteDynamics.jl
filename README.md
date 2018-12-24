| Testing | Coverage | Documentation |
| :-----: | :------: | :-----------: |
| [![Build Status](https://travis-ci.org/sisl/SatelliteDynamics.jl.svg?branch=master)](https://travis-ci.org/sisl/SatelliteDynamics.jl) | [![Coverage Status](https://coveralls.io/repos/github/sisl/SatelliteDynamics.jl/badge.svg?branch=master)](https://coveralls.io/github/sisl/SatelliteDynamics.jl?branch=master) | [![](https://img.shields.io/badge/docs-latest-blue.svg)](https://sisl.github.io/SatelliteDynamics.jl/latest) |

# SatelliteDynamics.jl
SatelliteDynamics package is meant to address the needs of the satellite operator, academic researcher, and public enthusiast communities.

Current high-fidelity satellite dynamics modeling software generally falls into two camps:
1. It is commercially licensed, closed-source code, making it difficult if not prohibitive to be used by hobbiests, academic researchers, or anyone.
2. There is a steep learning curve which makes correct use of the underlying libraries difficult.

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

More example code and use will be added as time permits.