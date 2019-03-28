# Propagators

The Propagators module aims to provide straight forward, easy to use, interfaces for simulations of orbit and attitude dynamics. The module both defines the state 
derivative dynamics functions as well as methods to propagate the state vector
in time.

It is important to note the default behavior and documentation for each function
because functionality can change dramatically (simulation of state-only for full
integration of the variational equations) depending on the inputs provided.

```@docs
fderiv_earth_orbit
EarthInertialState
step!
stepto!
sim!
reinit!
```