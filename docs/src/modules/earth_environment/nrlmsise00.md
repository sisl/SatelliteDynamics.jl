# NRLMSISE00

The NRLMSISE00 module provides a pure-Julia implementation of the NRLMSISE00 
Earth atmospheric model.

The implementation is based off of Daniel Brodo's C implementation which can be 
found [here](https://www.brodo.de/space/nrlmsise/).

Most of the module code is for purple internal only. The density should only
be retrieved through the `density_nrlmsise00` function which will properly 
compute the required model inputs for the given time.

Unforunately because the Geomagnetic and solar flux data has no predicted 
component this density model is currently restricted to computing density in the
past.

```@docs
density_nrlmsise00
```