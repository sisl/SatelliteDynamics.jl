# Universe

The Universe submodule defines simulation-specific data files which are constants 
of most simulations. In particular it provides data structures for storing and
accessing Earth orientation parameters and spherical harmonic gravity field 
models.

The module defines the global variables `EOP` and `GRAVITY_MODEL` which load
default data values at runtime.

`EOP` defaults to use the rapid Earth orientation data file `finals.all (IAU 2000)`
distributed by the IERS. The module also supports IERS C04 product files.

`GRAVITY_MODEL` defaults to use the EGM2008 spherical harmonic gravity model, 
truncated to order and degree 90.

All data files in the module can be updated by running the command `download_all_data()` in the Julia REPL.

```@docs
download_kp
download_solar_flux
download_all_data
EarthOrientationData
EOP
UT1_UTC
POLE_LOCATOR
XP
YP
set_eop
load_eop
GravModel
GRAVITY_MODEL
load_gravity_model
```