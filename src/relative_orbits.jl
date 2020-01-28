##############
# RTN | LVLH #
##############

export rRTNtoECI
export rECItoRTN


"""
Compute the radial, along-track, cross-track (RTN) to Earth-centered inertial
rotation matrix. If applied to a position vector in the RTN frame, it will
transform that vector to into the equivalent position vector in the ECI frame.

The RTN frame is also commonly refered to as the local-vertical, local-horizontal (LVLH) frame.

Arguments:
- `x::AbstractVector{<:Real}`: Inertial state (position and velocity) of primary (observing) satellite

Returns:
- `R_rtn2eci::SMatrix{3,3}`: Rotation matrix transforming _from_ the RTN frame _to_ the ECI frame
"""
function rRTNtoECI(x::AbstractVector{<:Real})
    r = x[idx1t3]
    v = x[idx4t6]
    n = cross(r, v)

    R = normalize(r)
    N = normalize(n)
    T = cross(N, R)

    R_rtn2eci = hcat(R, T, N)

    return R_rtn2eci
end


"""
Compute the Earth-centered inertial to radial, along-track, cross-track (RTN)
rotation matrix. If applied to a position vector in the ECI frame, it will
transform that vector into the equivalent position vector in the RTN frame.

The RTN frame is also commonly refered to as the local-vertical, local-horizontal (LVLH) frame.

Arguments:
- `x::AbstractVector{<:Real}`: Inertial state (position and velocity) of primary (observing) satellite

Returns:
- `R_eci2rtn::SMatrix{3,3}`: Rotation matrix transforming _from_ the ECI frame _to_ the RTN frame
"""
function rECItoRTN(x::AbstractVector{<:Real})
    r = x[idx1t3]
    v = x[idx4t6]
    n = cross(r, v)

    R = normalize(r)
    N = normalize(n)
    T = cross(N, R)

    R_eci2rtn = hcat(R, T, N)'

    return R_eci2rtn
end


export sECItoRTN
"""
Compute the radial, along-track, cross-track (RTN) coordinates of a target satellite in the primary satellites RTN frame.

The RTN frame is also commonly refered to as the local-vertical, local-horizontal (LVLH) frame.

Arguments:
- `x::AbstractVector{<:Real}`: Inertial state (position and velocity) of primary (observing) satellite
- `xt::AbstractVector{<:Real}`: Inertial state (position and velocity) of the target satellite

Returns:
- `xrtn::AbstractVector{<:Real}`: Position and velocity of the target relative of the observing satellite in the RTN.
"""
function sECItoRTN(x::AbstractVector{<:Real}, xt::AbstractVector{<:Real})
    if length(xt) >= 6
        return sECItoRTN(x, xt[idx1t6])
    else
        return sECItoRTN(x, xt[idx1t3])
    end
end

"""
Arguments:
- `x::AbstractVector{<:Real}`: Inertial state (position and velocity) of primary (observing) satellite
- `xt::SVector{3,<:Real}`: Inertial state (position) of the target satellite

Returns:
- `xrtn::SVector{3,<:Real}`: Position of the target relative of the observing satellite in the RTN.
"""
function sECItoRTN(x::AbstractVector{<:Real}, xt::SVector{3,<:Real})
    # Create RTN rotation matrix
    R_eci2rtn = rECItoRTN(x)

    # Transform Position
    r     = x[idx1t3]
    rho   = xt[idx1t3] - r
    r_rtn = R_eci2rtn*rho

    return r_rtn
end

"""
Arguments:
- `x::AbstractVector{<:Real}`: Inertial state (position and velocity) of primary (observing) satellite
- `xt::SVector{6,<:Real}`: Inertial state (position and velocity) of the target satellite

Returns:
- `xrtn::SVector{6,<:Real}`: Position and velocity of the target relative of the observing satellite in the RTN.
"""
function sECItoRTN(x::AbstractVector{<:Real}, xt::SVector{6,<:Real})
    # Create RTN rotation matrix
    R_eci2rtn = rECItoRTN(x)

    # Transform Position
    r     = x[idx1t3]
    rho   = xt[idx1t3] - r
    r_rtn = R_eci2rtn*rho

    v       = x[idx4t6]
    f_dot   = norm(cross(r, v)) / norm(r)^2
    omega   = SVector(0, 0, f_dot)
    rho_dot = xt[idx4t6] - v
    v_rtn   = R_eci2rtn*rho_dot - cross(omega, r_rtn)

    return vcat(r_rtn, v_rtn)
end


export sRTNtoECI
"""
Compute the Earth-center

The RTN frame is also commonly refered to as the local-vertical, local-horizontal (LVLH) frame.

Arguments:
- `x::AbstractVector{<:Real}`: Inertial state (position and velocity) of primary (observing) satellite
- `xrtn::AbstractVector{<:Real}`: Inertial state (position and velocity) of the target satellite

Returns:
- `xt::AbstractVector{<:Real}`: Position and velocity of the target relative of the observing satellite in the RTN.
"""
function sRTNtoECI(x::AbstractVector{<:Real}, xrtn::AbstractVector{<:Real})
    if length(xrtn) >= 6
        return sRTNtoECI(x, xrtn[idx1t6])
    else
        return sRTNtoECI(x, xrtn[idx1t3])
    end
end

"""
Arguments:
- `x::AbstractVector{<:Real}`: Inertial state (position and velocity) of primary (observing) satellite
- `xrtn::SVector{3,<:Real}`: Inertial state (position) of the target satellite

Returns:
- `xt::SVector{3,<:Real}`: Position of the target relative of the observing satellite in the RTN.
"""
function sRTNtoECI(x::AbstractVector{<:Real}, xrtn::SVector{3,<:Real})
    # Create RTN rotation matrix
    R_rtn2eci = rRTNtoECI(x)

    # Transform position
    r   = x[idx1t3]
    r_t = R_rtn2eci*xrtn + r

    return r_t
end

"""
Arguments:
- `x::AbstractVector{<:Real}`: Inertial state (position and velocity) of primary (observing) satellite
- `xrtn::SVector{6,<:Real}`: Inertial state (position and velocity) of the target satellite

Returns:
- `xt::SVector{6,<:Real}`: Position and velocity of the target relative of the observing satellite in the RTN.
"""
function sRTNtoECI(x::AbstractVector{<:Real}, xrtn::SVector{6,<:Real})
    # Create RTN rotation matrix
    R_rtn2eci = rRTNtoECI(x)

    # Transform position
    r     = x[idx1t3]
    r_rtn = xrtn[idx1t3]
    r_t   = R_rtn2eci*r_rtn + r

    # Transform velocity
    v     = x[idx4t6]
    v_rtn = xrtn[idx4t6]
    f_dot = norm(cross(r, v)) / norm(r)^2
    omega = SVector(0, 0, f_dot)
    v_t   = R_rtn2eci*(v_rtn + cross(omega, r_rtn)) + v

    return vcat(r_t, v_t)
end

#######
# ROE #
#######

"""
Convert Keplerian relative orbital element state to ROE relative state
"""
function sOSCtoROE(oe_c::SVector{6,<:Real}, oe_d::SVector{6,<:Real})
end

"""
Convert ROE relative state and chief relative orbital elements and return
osculating orbital elements
"""
function sROEtoOSC(oe_c::SVector{6,<:Real}, roe_d::SVector{6,<:Real})
end

"""
Convert inertial states to relative orbital element 
"""
function sECItoROE(x_c::SVector{6,<:Real}, x_d::SVector{6,<:Real})
end

"""
Convert ROE relative state and chief absolute inertial state and return
osculating orbital elements
"""
function sROEtoECI(x_c::SVector{6,<:Real}, roe_d::SVector{6,<:Real})
end

"""
Convert RTN state and chief absolute inertial state and return relative
orbital element state
"""
function sRTNtoROE(x_c::SVector{6,<:Real}, rtn_d::SVector{6,<:Real})
end

"""
Covert relative orbital element state and chief inertial state and return
RTN state.
"""
function sROEtoRTN(x_c::SVector{6,<:Real}, roe_d::SVector{6,<:Real})
end