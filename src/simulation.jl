__precompile__(true)
module Simulation

using SatelliteDynamics.Time: Epoch
using SatelliteDynamics.OrbitDynamics: deriv_orbit_earth, accel_point_mass
using OrdinaryDiffEq # Don't use all of DifferentialEquations due to long precompile time

"""
f_orbit_dynamics is an internal helper function used to interface the OrbitDynamics module with DifferentialEquations.jl solvers.

Arguments:
- `u::Array{<:Real, 1}`: Current satellite inertial state
- `p::Array{<:Real, 1}`: Parameter vector
- `t::Real`: Time offset from initial time

Returns:
- `du::Array{<:Real, 1}`: Curre]nt state derivative

Parameters:
1. `epc0::Epoch`: Initial Epoch
2. `mass::Real`: Satellite mass [kg]
3. `Ad::Real`: Area of drag (cannon-ball model) [m^2]
4. `Cd::Real`: Coefficient of drag [dimensionless]
5. `Asrp::Real`: Area of solar radiation prerssure [m^2]
6. `Cr::Real`: Coefficient of reflectivity [dimensionless]
7. `n_grav::Integer`: Degree of gravity field [dimensionless]
8. `m_grav::Integer`: Order of gravity field [dimensionless]
9. `drag::Bool`: Apply drag perturbation
10. `srp::Bool`: Apply solar radiation pressure perturbation
11. `moon::Bool`: Apply third-body lunar gravity perturbation
12. `sun::Bool`: Apply third-body solar gravity perturbation
13. `relativity::Bool`: Apply relativistic perturbations
"""
function f_orbit_dynamics(u, p, t)
    epc0 = p[1]
    epc  = epc0 + t
    du = deriv_orbit_earth(epc, u, mass=p[2], area_drag=p[3], coef_drag=p[4], 
                            area_srp=p[5], coef_srp=p[6], n_grav=p[7], 
                            m_grav=p[8], drag=p[9], srp=p[10], moon=p[11], sun=p[12], relativity=p[13])
    
    return du
end

export propagate_orbit
"""
Simulate orbit dynamics

Arguments:
- `epc0::Epoch`: Propagation start Epoch
- `eci0::Epoch`: Initial Cartesean inertial state [m; m/s]
- `epcf::Epoch`: Final to simulate to.
- `timestep::Real`: Timestep to use for simulation (Default: 5.0)
- `solver`: Solver to use to solve ODEProblem
- `atol::Real`: Absolute tolerate limit for differential equation solution (Default: 1.0e-9)
- `atol::Real`: Absolute tolerate limit for differential equation solution (Default: 1.0e-9)
- `mass::Real`: Satellite mass [kg] (Default: 100.0)
- `area_drag::Real`: Area of drag (cannon-ball model) [m^2] (Default: 1.0)
- `coef_drag::Real`: Coefficient of drag [dimensionless] (Default: 2.3)
- `arear_srp::Real`: Area of solar radiation prerssure [m^2] (Default: 1.0)
- `coef_srp::Real`: Coefficient of reflectivity [dimensionless] (Default: 1.8)
- `n_grav::Integer`: Degree of gravity field [dimensionless] (Default: 0.0)
- `m_grav::Integer`: Order of gravity field [dimensionless] (Default: 0.0)
- `drag::Bool`: Apply drag perturbation (Default: false)
- `srp::Bool`: Apply solar radiation pressure perturbation (Default: false)
- `moon::Bool`: Apply third-body lunar gravity perturbation (Default: false)
- `sun::Bool`: Apply third-body solar gravity perturbation (Default: false)
- `relativity::Bool`: Apply relativistic perturbations (Default: false)

Returns:
- `t::Array{Float64, 1}`: Simulation output times at elapsed seconds from initial Epoch
- `epc::Array{Epoch, 1}`: Simulation output times as absolute Epochs
- `eci::Array{Float64, 2}`: Propgated inertial state. Time is aligned with column dimension.
"""
function propagate_orbit(epc0::Epoch, x0::Array{<:Real, 1}, epcf::Epoch; timestep=5.0, dtmax=60.0, solver=RK4(), rtol=1.0e-9, atol=1.0e-9, mass=100.0::Real, area_drag=1.0::Real, coef_drag=2.3::Real, arear_srp=1.0::Real, coef_srp=1.8::Real, n_grav=0::Integer, m_grav=0::Integer, drag=false::Bool, srp=false::Bool, moon=false::Bool, sun=false::Bool, relativity=false::Bool)
    tspan  = (0, epcf-epc0)
    params = [epc0, mass, area_drag, coef_drag, arear_srp, coef_srp, n_grav, m_grav, drag, srp, moon, sun, relativity]
    prob   = ODEProblem(f_orbit_dynamics, x0, tspan, params)

    sol = solve(prob, solver, rtol=rtol, atol=atol, save_start=true, save_end=true, saveat=timestep, dtmax=dtmax)

    n_step = length(sol.t)
    t      = zeros(Float64, n_step)
    epc    = Array{Epoch, 1}(undef, n_step)
    eci    = zeros(Float64, 6, n_step)

    for i in 1:n_step
        t[i]      = sol.t[i]
        epc[i]    = epc0 + t[i]
        eci[:, i] = sol[i]
    end

    return t, epc, eci
end

end # End module