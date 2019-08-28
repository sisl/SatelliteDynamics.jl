# Exports
export RK4, istep

#######
# RK4 #
#######

"""
RK4 is a numerical integrator object which stores a function of signature 
f(epc, x,...; kwargs...) along with its associated settings parameters.

Arguments:
- `f::Function` Function for integration. 
- `auxParams::Dict`
"""
struct RK4
    f::Function
    auxParams
end

function RK4(f::Function; kwargs...)
    return RK4(f, kwargs)
end

"""
Internal funciton used to compute the variational matrix.

Arguments:
- `rk4::RK4` 4-th order Runge-Kutta numerical integrator object
- `epc::Union{Real, Epoch}` Absolute time of start of integration step
- `x::Array{<:Real, 1}` State vector for linearizing perturbations around
- `apert::Union{Real, Array{<:Real,1}}` Absolute perturbations to state elements
used to numerically calculate the partial derivatives through a forward difference
method.

Returns:
- `V:Array{Float64, 2}` Variational equation drivative matrix. This is the 
matrix of partial derivatives of the dynamics function _f_ with respect to the 
input state _x_.
"""
function varmat(rk4::RK4, epc::Union{Real,Epoch}, x::Array{<:Real, 1}; apert=1.0::Union{Real, Array{<:Real,1}})
    # Get state length
    n = length(x)

    # Expand perturbation step-size as vector
    if typeof(apert) <: Real
        apert = ones(typeof(apert), n)
    end

    # Initialize variational matrix
    V = zeros(Float64, n, n)

    # Evaluate unperturbed state derivative function
    fx = rk4.f(epc, x; rk4.auxParams...)

    # Compute perturbed acceleration for each field
    for i in 1:n
        # Perturbed state vector
        px     = copy(x)
        px[i] += apert[i]

        # Perturbed state derivative
        pfx = rk4.f(epc, px; rk4.auxParams...)

        # Compute partial of state derivative with respect to state using 
        # forward difference method
        dfdx = (pfx - fx)/apert[i]

        # Set i-th column as state partial 
        V[:, i] = dfdx 
    end

    return V
end

function istep(rk4::RK4, epc::Union{Real,Epoch}, dt::Real, x::Array{<:Real, 1})
    # Compute State coefficients

    # Perform internal steps
    k1 = rk4.f(epc,                      x; rk4.auxParams...)
    k2 = rk4.f(epc + dt/2.0, x+(dt/2.0)*k1; rk4.auxParams...)
    k3 = rk4.f(epc + dt/2.0, x+(dt/2.0)*k2; rk4.auxParams...)
    k4 = rk4.f(epc + dt,           x+dt*k3; rk4.auxParams...)

    # Compute updated state
    xu = x + (dt/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4)

    return xu
end

"""
Performs single integration step of numerical integrator.

Arguments:
- `rk4::RK4` 4-th order Runge-Kutta numerical integrator object
- `epc::Union{Real, Epoch}` Absolute time of start of integration step
- `dt::Real` Integration step size
- `x::Array{<:Real, 1}` State vector
- `phi::Array{<:Real, 2}` State transition matrix at start of integration step
"""
function istep(rk4::RK4, epc::Union{Real,Epoch}, dt::Real, x::Array{<:Real, 1}, phi::Array{<:Real, 2}; apert=1.0::Union{Real, Array{<:Real,1}})
    # Compute State coefficients

    # Perform state integration
    xk1 = rk4.f(epc,                         x; rk4.auxParams...)
    xk2 = rk4.f(epc + dt/2.0, x + (dt/2.0)*xk1; rk4.auxParams...)
    xk3 = rk4.f(epc + dt/2.0, x + (dt/2.0)*xk2; rk4.auxParams...)
    xk4 = rk4.f(epc + dt,           x + dt*xk3; rk4.auxParams...)

    # Perform Integration of variational equations
    phik1 = varmat(rk4, epc, x, apert=apert)*phi
    phik2 = varmat(rk4, epc + dt/2.0, x + dt/2.0*xk1, apert=apert)*(phi + (dt/2.0)*phik1)
    phik3 = varmat(rk4, epc + dt/2.0, x + dt/2.0*xk2, apert=apert)*(phi + (dt/2.0)*phik2)
    phik4 = varmat(rk4,     epc + dt,     x + dt*xk3, apert=apert)*(phi + dt*phik3)
 
    # Compute updated state
    xu   = x   + (dt/6.0)*(xk1   + 2.0*xk2   + 2.0*xk3   + xk4)
    phiu = phi + (dt/6.0)*(phik1 + 2.0*phik2 + 2.0*phik3 + phik4)

    return xu, phiu
end