# Space Weather

The Space Weather submodule provides classes to store and access space weather 
data files. In particular Solar Flux and Geomagnetic index data.

The module defines the global variables `GEOMAGNETIC_DATA` and `SOLAR_FLUX_DATA` which are loaded at runtime and 

`GEOMAGNETIC_DATA` stores Geomagnetic Kp and Ap indicies used in models of 
Earth's atmosphere.

`SOLAR_FLUX_DATA` stores F10.7 cm solar flux data which is a measurement of 
solar activity and another input to most atmospheric models.

```@docs
SpaceWeatherData
SPACE_WEATHER_DATA
load_space_weather_data
KpIndices
KpIndex
KpDailyIndex
ApIndices
ApIndex
ApDailyIndex
f107Data
f107Observed
f107Adjusted
f107ObservedAvg
f107AdjustedAvg
```