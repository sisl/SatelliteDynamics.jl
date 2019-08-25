__precompile__(true)
module ReferenceSystems

using LinearAlgebra

using StaticArrays: SVector, SMatrix

using SOFA
using SatelliteDynamics.Constants
using SatelliteDynamics.Universe: UT1_UTC, POLE_LOCATOR
using SatelliteDynamics.Time: Epoch, mjd

const idx123    = SVector(1, 2, 3)
const idx456    = SVector(4, 5, 6)
const idx123456 = SVector(1, 2, 3, 4, 5, 6)


##############
# RTN | LVLH #
##############

export rRTNtoECI
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
    r = x[idx123]
    v = x[idx456]
    n = cross(r, v)

    R = normalize(r)
    N = normalize(n)
    T = cross(N, R)

    R_rtn2eci = hcat(R, T, N)

    return R_rtn2eci
end


export rECItoRTN
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
    r = x[idx123]
    v = x[idx456]
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
        return sECItoRTN(x, xt[idx123456])
    else
        return sECItoRTN(x, xt[idx123])
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
    r     = x[idx123]
    rho   = xt[idx123] - r
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
    r     = x[idx123]
    rho   = xt[idx123] - r
    r_rtn = R_eci2rtn*rho

    v       = x[idx456]
    f_dot   = norm(cross(r, v)) / norm(r)^2
    omega   = SVector(0, 0, f_dot)
    rho_dot = xt[idx456] - v
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
        return sRTNtoECI(x, xrtn[idx123456])
    else
        return sRTNtoECI(x, xrtn[idx123])
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
    r   = x[idx123]
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
    r     = x[idx123]
    r_rtn = xrtn[idx123]
    r_t   = R_rtn2eci*r_rtn + r

    # Transform velocity
    v     = x[idx456]
    v_rtn = xrtn[idx456]
    f_dot = norm(cross(r, v)) / norm(r)^2
    omega = SVector(0, 0, f_dot)
    v_t   = R_rtn2eci*(v_rtn + cross(omega, r_rtn)) + v

    return vcat(r_t, v_t)
end


#######################################
# IAU 2010 | Inertial <-> Earth-Fixed #
#######################################

export bias_precession_nutation
"""
Computes the Bias-Precession-Nutation matrix transforming the GCRS to the 
CIRS intermediate reference frame. This transformation corrects for the 
bias, precession, and nutation of Celestial Intermediate Origin (CIO) with
respect to inertial space.

Arguments:
- `epc::Epoch`: Epoch of transformation

Returns:
- `rc2i::Matrix{<:Real}`: 3x3 Rotation matrix transforming GCRS -> CIRS
"""
function bias_precession_nutation(epc::Epoch)
    # Constants of IAU 2006A transofrmation
    DMAS2R =  4.848136811095359935899141e-6 / 1.0e3
    dx06   =  0.0001750*DMAS2R
    dy06   = -0.0002259*DMAS2R

    # Compute X, Y, s terms using low-precision series terms
    x, y, s = iauXys00b(MJD_ZERO, mjd(epc, tsys="TT"))

    # Apply IAU2006 Offsets
    x += dx06
    y += dy06

    # Compute transformation and return
    rc2i = iauC2ixys(x, y, s)

    return rc2i
end


export earth_rotation
"""
Computes the Earth rotation matrix transforming the CIRS to the TIRS
intermediate reference frame. This transformation corrects for the Earth
rotation.

Arguments:
- `epc::Epoch`: Epoch of transformation

Returns:
- `r::Matrix{<:Real}`: 3x3 Rotation matrix transforming CIRS -> TIRS
"""
function earth_rotation(epc::Epoch)
    # Compute Earth rotation angle
    era = iauEra00(MJD_ZERO, mjd(epc, tsys="UT1"))

    # Rotate Matrix and return
    r = iauRz(era, Matrix{Float64}(I, 3, 3))

    return r
end


export polar_motion
"""
Computes the Earth rotation matrix transforming the TIRS to the ITRF reference 
frame.

# Arguments
- `epc::Epoch`: Epoch of transformation

# Returns
- `rpm::Matrix{<:Real}`: 3x3 Rotation matrix transforming TIRS -> ITRF
"""
function polar_motion(epc::Epoch)
    xp, yp = POLE_LOCATOR(mjd(epc, tsys="UTC"))

    # Compute transformation and return
    rpm = iauPom00(xp, yp, iauSp00(MJD_ZERO, mjd(epc, tsys="TT")))

    return rpm
end


export rECItoECEF
"""
Computes the combined rotation matrix from the inertial to the Earth-fixed
reference frame. Applies corrections for bias, precession, nutation,
Earth-rotation, and polar motion.

The transformation is accomplished using the IAU 2006/2000A, CIO-based 
theory using classical angles. The method as described in section 5.5 of 
the SOFA C transformation cookbook.

# Arguments
- `epc::Epoch`: Epoch of transformation

# Returns
- `r::Matrix{<:Real}`: 3x3 Rotation matrix transforming GCRF -> ITRF
"""
function rECItoECEF(epc::Epoch)
    # Compute intermediate transformations
    rc2i = bias_precession_nutation(epc)
    r    = earth_rotation(epc)
    rpm  = polar_motion(epc) 

    return rpm * r * rc2i
end


export rECEFtoECI
"""
Computes the combined rotation matrix from the Earth-fixed to the inertial
reference frame. Applies corrections for bias, precession, nutation,
Earth-rotation, and polar motion.

The transformation is accomplished using the IAU 2006/2000A, CIO-based 
theory using classical angles. The method as described in section 5.5 of 
the SOFA C transformation cookbook.

# Arguments
- `epc::Epoch`: Epoch of transformation

# Returns
- `r::Matrix{<:Real}`: 3x3 Rotation matrix transforming ITRF -> GCRF
"""
function rECEFtoECI(epc::Epoch)
    # Compute intermediate transformations
    rc2i = bias_precession_nutation(epc)
    r    = earth_rotation(epc)
    rpm  = polar_motion(epc) 

    return rc2i' * r' * rpm'
end


export sECItoECEF
"""
Transforms an Earth inertial state into an Earth fixed state

The transformation is accomplished using the IAU 2006/2000A, CIO-based 
theory using classical angles. The method as described in section 5.5 of 
the SOFA C transformation cookbook.

# Arguments
- `epc::Epoch`: Epoch of transformation
- `x::AbstractVector{<:Real}`: Inertial state (position, velocity) [m; m/s]

# Returns
- `x_ecef::Vector{<:Real}`: Earth-fixed state (position, velocity)
"""
function sECItoECEF(epc::Epoch, x::AbstractVector{<:Real})
    dim_x  = length(x)
    x_ecef = zeros(Float64, dim_x)

    # Extract State Components
    r_eci = x[1:3]

    if dim_x >= 6
        v_eci = x[4:6]
    end

    if dim_x == 9
        a_eci = x[7:9]
    end

    # Compute Sequential Transformation Matrices
    rc2i = bias_precession_nutation(epc)
    r    = earth_rotation(epc)
    pm   = polar_motion(epc)

    # Create Earth's Angular Rotation Vector
    omega_vec = [0, 0, OMEGA_EARTH] # Neglect LOD effect

    # Calculate ECEF State
    x_ecef[1:3] = pm *  r * rc2i * r_eci

    if dim_x == 6
        x_ecef[4:6] = pm * (r * rc2i * v_eci - cross(omega_vec, r * rc2i * r_eci))
    end
    
    if dim_x == 9
        x_ecef[7:9] = pm * (r * rc2i * a_eci - cross(omega_vec, cross(omega_vec, r * rc2i * r_eci)) 
                                         - 2 * cross(omega_vec, r * rc2i * v_eci))
    end

    return x_ecef
end


export sECEFtoECI
"""
Transforms an Earth fixed state into an Inertial state

The transformation is accomplished using the IAU 2006/2000A, CIO-based 
theory using classical angles. The method as described in section 5.5 of 
the SOFA C transformation cookbook.

# Arguments
- `epc::Epoch`: Epoch of transformation
- `x::AbstractVector{<:Real}`: Earth-fixed state (position, velocity) [m; m/s]

# Returns
- `x_eci::Vector{<:Real}`: Inertial state (position, velocity)
"""
function sECEFtoECI(epc::Epoch, x::AbstractVector{<:Real})
    # Set state variable size
    dim_x = length(x)
    x_eci = zeros(Float64, dim_x)

    # Extract State Components
    r_ecef = x[1:3]

    if dim_x >= 6
        v_ecef = x[4:6]
    end

    if dim_x == 9
        a_ecef = x[7:9]
    end

    # Compute Sequential Transformation Matrices
    bpn = bias_precession_nutation(epc)
    rot = earth_rotation(epc)
    pm  = polar_motion(epc)

    # Create Earth's Angular Rotation Vector
    omega_vec = [0, 0, OMEGA_EARTH] # Neglect LOD effect
    
    # Calculate ECEF State
    x_eci[1:3] = (pm * rot * bpn)' * r_ecef

    if dim_x >= 6
        x_eci[4:6] = (rot * bpn)' * ((pm)' * v_ecef + cross(omega_vec, (pm)' * r_ecef))
    end

    if dim_x >= 9
        x_eci[7:9] = (rot * bpn)' * ((pm)' * a_ecef + cross(omega_vec, cross(omega_vec, (pm)' * x_eci[1:3])) 
                                 + 2 * cross(omega_vec, (pm)' * x_eci[4:6])) 
    end

    return x_eci
end

end
