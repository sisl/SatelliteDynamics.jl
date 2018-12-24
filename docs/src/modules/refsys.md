# ReferenceSystems

The ReferenceSystems submodule contains precision transformations between common standard reference systems.

Currently only Earth-based reference systems are supported.

Most of the transformations are make backend calls to the SOFA C-library functions provide the package [SOFA.jl](https://github.com/sisl/SOFA.jl "SOFA.jl")

```@docs
rotation_rtn_to_eci
rotation_eci_to_rtn
state_eci_to_rtn
state_rtn_to_eci
bias_precession_nutation
earth_rotation
polar_motion
rotation_eci_ecef
rotation_ecef_eci
state_eci_to_ecef
state_ecef_to_eci
```