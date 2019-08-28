###################
# Dynamics Models #
###################

export fderiv_earth_orbit
"""
Compute the state derivative.

Arguments:
- `epc::Epoch`: Current epoch
- `x::Array{<:Real, 1}`: Satellite state vector
- `mass::Real`: Satellite mass [kg]
- `area_drag`: Velocity-facing area affected by drag. [m^2]
- `coef_drag`: Coefficient of drag [dimensionless]
- `area_srp`: Velocity-facing area affected by drag. [m^2]
- `coef_srp`: Coefficient of reflectivity [dimensionless]  
- `n_grav::Integer`: Gravity model degree (Default: 20)
- `m_grav::Integer`: Gravity model order (Default: 20)
- `drag::Bool`: Include cannonball atomospheric drag in force model (Default: `true`)
- `srp::Bool`: Include flat-plate solar radiation pressure in force model (Default: `true`)
- `moon::Bool`: Include thridbody lunar gravity in force model (Default: `true`)
- `sun::Bool`: Include thirdbody solar in force model (Default: `true`)
- `relativity::Bool`: Include relativistic effects in force model (Default: `true`)

Returns:
- `dx::Array{<:Float64, 1}`: Satellite state derivative, velocity and accelerations [m; m/s]
"""
function fderiv_earth_orbit(epc::Epoch, x::Array{<:Real} ;
             mass::Real=1.0, area_drag::Real=1.0, coef_drag::Real=2.3, 
             area_srp::Real=1.0, coef_srp::Real=1.8, 
             n_grav::Integer=20, m_grav::Integer=20, 
             drag::Bool=true, srp::Bool=true, moon::Bool=true, sun::Bool=true, 
             relativity::Bool=true)
    
    # Extract position and velocity
    r = x[1:3]
    v = x[4:6]

    # Compute ECI to ECEF Transformation -> IAU2010 Theory
    PN = bias_precession_nutation(epc)
    E  = earth_rotation(epc)
    W  = polar_motion(epc)
    R  = W * E * PN

    # Compute geolocation
    geod = sECEFtoGEOD(R*r)

    # Compute sun and moon position
    r_sun  = sun_position(epc)
    r_moon = moon_position(epc)

    # Compute acceleration
    a = zeros(Float64, 3)

    a += accel_gravity(x, R, n_grav, m_grav)

    # Drag
    if drag
        # Use PN*x[1:3] to compute the satellite position in the true-of-date inertial frame
        rho = density_nrlmsise00(epc, geod)
        a  += accel_drag(x, rho, mass, area_drag, coef_drag, Array{Real, 2}(PN))
    end

    # Solar Radiation Pressure
    if srp
        nu = eclipse_conical(x, r_sun)
        a += nu*accel_srp(x, r_sun, mass, area_srp, coef_srp)
    end

    # Thirdbody
    if sun
        a += accel_thirdbody_sun(x, r_sun)
    end

    if moon
        a += accel_thirdbody_moon(x, r_moon)
    end

    # Relativity
    if relativity
        a += accel_relativity(x)
    end

    return vcat(v, a)
end

##############
# Propagator #
##############

export EarthInertialState
"""
Satellite orbit state vector using an Earth inertial state representation and
dynamics model.

If an initial state transition matrix is provided it will be used for propagation
if there is no state transition matrix then only the state will be propagated.

Attributes:
- `rk4::RK4` Internal numerical integrator used for state propagation
- `dt::Real` Default propagation time step. All steps will be this size unless 
the state vector is requested to propagate to a time smaller than this step size,
which it will do.
- `epc::Epoch` Epoch of state
- `x::Array{Float64, 1}` State vector. Earth-centered inertial Cartesian state.
- `phi::Union{Nothing, Array{Float64, 2}}` State transition matrix, or the matrix
of partial derivatives of the state at the current time with respect to the 
start of propagation.

The following force model parametters can be set as keyword arguments

Parameters:
- `mass::Real`: Satellite mass [kg]
- `area_drag`: Velocity-facing area affected by drag. [m^2]
- `coef_drag`: Coefficient of drag [dimensionless]
- `area_srp`: Velocity-facing area affected by drag. [m^2]
- `coef_srp`: Coefficient of reflectivity [dimensionless]  
- `n_grav::Integer`: Gravity model degree (Default: 20)
- `m_grav::Integer`: Gravity model order (Default: 20)
- `drag::Bool`: Include cannonball atomospheric drag in force model (Default: `true`)
- `srp::Bool`: Include flat-plate solar radiation pressure in force model (Default: `true`)
- `moon::Bool`: Include thridbody lunar gravity in force model (Default: `true`)
- `sun::Bool`: Include thirdbody solar in force model (Default: `true`)
- `relativity::Bool`: Include relativistic effects in force model (Default: `true`)

"""
mutable struct EarthInertialState
    rk4::RK4
    dt::Real
    epc::Epoch
    x::Array{Float64, 1}
    phi::Union{Nothing, Array{Float64, 2}}
end

function EarthInertialState(epc::Epoch, x::Array{<:Real, 1}, phi::Union{Nothing, Array{Float64, 2}}=nothing; dt::Real=30.0, kwargs...)
    rk4 = RK4(fderiv_earth_orbit; kwargs...)
    return EarthInertialState(rk4, dt, epc, x, phi)
end

export step!
"""
Step dynamics state the specified time.

Arguments:
- `state::EarthInertialState` State vector to propagate
- `dt::Real` Step size. If no value is provided will use default state stepsize
"""
function step!(state::EarthInertialState, dt::Real=0.0)
    # Set step size
    if dt == 0.0
        dt = state.h
    end

    if state.phi == nothing
        # If no state transition matrix is provided integrate 
        # without propagating variational equations
        state.x = istep(state.rk4, state.epc, dt, state.x)
    else
        state.x, state.phi = istep(state.rk4, state.epc, dt, state.x, state.phi)
    end

    # Advance state propagator time
    state.epc += dt
end

export stepto!
"""
Step dynamics state to the specified time. Will take as many internal steps as 
required to advance the propagator to the correct time. No internal step will 
exceed the size specified in the state propagator.

Arguments:
- `state::EarthInertialState` State vector to propagate
- `time::Union{Real, Epoch}` Time to propagate internal state too. Can be either
a real number to advance the state by or the Epoch
"""
function stepto!(state::EarthInertialState, time::Union{Real, Epoch}=0.0)
    if typeof(time) <: Real
        time = state.epc + time
    end

    if time < state.epc
        @warn "Propagation time is before current state epoch. Will integrate dyanmics in reverse."
    end
    
    # Propagate state
    while state.epc < time
        # Set propagation state size
        dt = min(state.dt, time - state.epc)

        # Propagate state
        step!(state, dt)
    end
end

export sim!
"""
Simulate the sttate dynamics to the specified time. Takes same inputs as `stepto!`
but instead of just updating the state vector to the specified time, this 
function also returns the timesteps, state values, and state transition
matrices for each time step.

Arguments:
- `state::EarthInertialState` State vector to propagate
- `time::Union{Real, Epoch}` Time to propagate internal state too. Can be either
a real number to advance the state by or the Epoch

Returns:
- `t::Array{Float64, 1}` Elapsed time as a scalar from the initial simulation epoch 
- `epc::Array{Epoch, 1}` Epoch at each timestep 
- `x::Array{Float64, 2}` State vectors at each time step. Time is along second axis
- `Phi::Array{Float64, 2}` Stacked array of state transition matrices
"""
function sim!(state::EarthInertialState, time::Union{Real, Epoch}=0.0)
    if typeof(time) <: Real
        time = state.epc + time
    end

    if time < state.epc
        @warn "Propagation time is before current state epoch. Will integrate dyanmics in reverse."
    end

    # Number of state variables
    n = length(state.x)

    # Compute number of time steps in propagation
    n_steps = ceil(Int, (time - state.epc)/state.dt) + 1

    # Initialize containers to store output
    t   = zeros(Float64, n_steps)
    epc = Array{Epoch}(undef, n_steps)
    x   = zeros(Float64, n, n_steps)
    A   = zeros(Float64, n*n_steps, n)

    # Save initial state
    idx = 1

    epc[idx]  = state.epc
    t[idx]    = epc[idx] - epc[1]
    x[:, idx] = state.x

    if state.phi != nothing
        A[(1+n*(idx-1)):(n*idx), :] = state.phi
    end

    # Propagate state
    while state.epc < time
        # Increase storage index
        idx += 1

        # Set propagation state size
        dt = min(state.dt, time - state.epc)

        # Propagate state
        step!(state, dt)

        # Save output
        epc[idx]  = state.epc
        t[idx]    = epc[idx] - epc[1]
        x[:, idx] = state.x
        if state.phi != nothing
            A[(1+n*(idx-1)):(n*idx), :] = state.phi
        end
    end

    # Output
    if state.phi != nothing
        # Return with state transition matrix if used
        return t, epc, x, A
    end

    # Default return without STM 
    return t, epc, x
end

export reinit!
"""
Reinitialize State transition matrix to identity at the current time step.

Used to reinitialize the state transition matrix value to identity
"""
function reinit!(state::EarthInertialState)
    state.phi = diagm(0 => ones(Float64, length(state.x)))
end