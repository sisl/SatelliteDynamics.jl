# Simulation

The Simulation module provides the capability to simulate satellite orbit and 
attitude dynamics by combining the dynamics functions provided in this module 
with the numerical integration capabilities of DifferentialEquations.jl.

The goal is to provide straight forward, and easily customizable interfaces for simulations of orbit and attitude dynamics.

```@docs
propagate_orbit
```