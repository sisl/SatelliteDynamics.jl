# Time

The Time submodule contains common time transformations such as converting between different date representations or converting a specific instant in time between different time systems.

The module also defines the `Epoch` class which provides the standard representation of instants in time used throughout the module. 

Most of the transformations are make backend calls to the SOFA C-library functions provide the package [SOFA.jl](https://github.com/sisl/SOFA.jl "SOFA.jl")

```@docs
caldate_to_mjd
mjd_to_caldate
caldate_to_jd
jd_to_caldate
elapsed_from_epoch
days_from_elapsed
time_system_offset
Epoch
caldate
jd
mjd
day_of_year
gmst
gast
epoch_to_jdfd
```