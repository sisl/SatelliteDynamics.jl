# Universe

The Universe module defines simulation-specific data files which are constants 
of most simulations. In particular it provides data structures for storing and
accessing Earth orientation parameters and spherical harmonic gravity field 
models.

The module defines the global variables `EOP` and `GRAVITY_MODEL` which are
loaded at runtime and .

`EOP` defaults to use the rapid Earth orientation data file `finals.all (IAU 2000)`
distributed by the IERS. The module also supports IERS C04 product files.

`GRAVITY_MODEL` defaults to use the EGM2008 spherical harmonic gravity model, 
truncated to order and degree 90.

```@docs
EarthOrientationData
EOP
UT1_UTC
POLE_LOCATOR
XP
YP
set_eop
load_eop
update_eop
GravModel
GRAVITY_MODEL
load_gravity_model
```