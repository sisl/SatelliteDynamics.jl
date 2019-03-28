# Integrators

The Integrators module is meant to provide built in numerical integration 
capabilities for the module. While the dynamics functions provided as part of
of the Simulation module are compatible with the more complete set of numerical 
integrators provided in [DifferentialEquations.jl](https://github.com/JuliaDiffEq/DifferentialEquations.jl), the integrators implemented here make other useful
functionality such as integration of the Variational Equations, and incorporating
impulsive maneuvers easier to implement.

Currently only the `RK4` integrator is supported.

Normally it is easier to use the propagators defined in the Propagators module 
for simulation, but the class and method provided here can be used for the 
implementation of custom numerical propagators.

```@docs
RK4
istep
```