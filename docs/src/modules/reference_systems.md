# ReferenceSystems

The ReferenceSystems submodule contains precision transformations between common standard reference systems.

Currently only Earth-based reference systems are supported.

Most of the transformations are make backend calls to the SOFA C-library functions provide the package [SOFA.jl](https://github.com/sisl/SOFA.jl "SOFA.jl")

```@docs
rRTNtoECI
rECItoRTN
sECItoRTN
sRTNtoECI
bias_precession_nutation
earth_rotation
polar_motion
rECItoECEF
rECEFtoECI
sECItoECEF
sECEFtoECI
```