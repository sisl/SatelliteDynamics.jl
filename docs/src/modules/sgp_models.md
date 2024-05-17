# SGP Models

The `sgp_models` module contains the implementation of the Simplified General Perturbations (SGP) models for orbit propagation. These models are based on the work of [Vallado et al. (2001)](https://celestrak.com/publications/AIAA/2006-6753/). The SGP models are a set of simplified models for propagating the orbits of Earth-orbiting satellites. The models are based on the two-line element (TLE) format used by the US Department of Defense to distribute satellite orbit data. The SGP models are widely used in the satellite tracking community due to their simplicity and computational efficiency.

```@docs
TLE
state
ecef
eci
```