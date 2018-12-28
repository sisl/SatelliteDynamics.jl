# Orbit Propagation

One of the primary features that SatelliteDynamics.jl aims to make easily accessible to users, is the ability to perform high-fidelity orbit propagation in an easy and customizable manner. In the example we will show how to simulate a satellite orbits using the tools provided in this module.

To start out we will perform an orbit propagation using the most basic orbit 
model possible: a point-mass approximation of Earth's gravity without any 
other perturbation models added on. 

First, we must declare the initial conditions for the simulation. This entails 
declaring an initial `Epoch` as well as the inertial Cartesean state of the 
statellite at that Epoch. 

```julia
# Declare simulation initial Epoch
epc0 = Epoch(2019, 1, 1, 12, 0, 0, 0.0) 

# Declare initial state in terms of osculating orbital elements
oe0  = [R_EARTH + 500e3, 0.0, 90.0, 0, 0, 0]

# Convert osculating elements to Cartesean state
eci0 = sOSCtoCART(oe0, use_degrees=true)
```

Next, simulate the orbit:

```julia
# Set the propagation end time to one orbit period after the start
epcf = epc0 + orbit_period(oe0[1])

# Propagate the orbit
t, epc, eci = simulate(orb, epcf, timestep=1, dtmax=1)
```

Putting it all together we have:

```julia
# Declare simulation initial Epoch
epc0 = Epoch(2019, 1, 1, 12, 0, 0, 0.0) 

# Declare initial state in terms of osculating orbital elements
oe0  = [R_EARTH + 500e3, 0.0, 90.0, 0, 0, 0]

# Convert osculating elements to Cartesean state
eci0 = sOSCtoCART(oe0, use_degrees=true)

# Set the propagation end time to one orbit period after the start
epcf = epc0 + orbit_period(oe0[1])

# Propagate the orbit
t, epc, eci = simulate(epc0, eci0, epcf, timestep=1, dtmax=1)
```

And that's it! All it took was 5 lines of code with the SatelliteDynamics 
module to propagate an orbit. 