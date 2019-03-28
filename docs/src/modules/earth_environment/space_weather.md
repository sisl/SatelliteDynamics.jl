# Space Weather

The Space Weather submodule provides classes to store and access space weather 
data files. In particular Solar Flux and Geomagnetic index data.

The module defines the global variables `GEOMAGNETIC_DATA` and `SOLAR_FLUX_DATA` which are loaded at runtime and 

`GEOMAGNETIC_DATA` stores Geomagnetic Kp and Ap indicies used in models of 
Earth's atmosphere.

`SOLAR_FLUX_DATA` stores F10.7 cm solar flux data which is a measurement of 
solar activity and another input to most atmospheric models.

```@docs
GeomagneticIndexData
GEOMAGNETIC_DATA
KpIndices
KpIndex
KpDailyIndex
ApIndices
ApIndex
ApDailyIndex
SolarFluxData
SOLAR_FLUX_DATA
f107Data
f107Observed
f107Adjusted
f107ObservedAvg
f107AdjustedAvg
```