# Coordinates

The coordiantes module provides functionatiy for transforming between various
body-fixed reference frame coordinate systems. In particular, geodetic and geocentric transformations are included. Topocentric transformations are also included.

```@docs
geocentric_to_ecef
ecef_to_geocentric
geodetic_to_ecef
ecef_to_geodetic
rotation_ecef_to_enz
rotation_enz_to_ecef
state_ecef_to_enz
state_enz_to_ecef
rotation_ecef_to_sez
rotation_sez_to_ecef
state_ecef_to_sez
state_sez_to_ecef
enz_to_azel
sez_to_azel
```