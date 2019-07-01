__precompile__(true)
module ReferenceSystems

using LinearAlgebra

using SOFA
using SatelliteDynamics.Constants
using SatelliteDynamics.Universe: UT1_UTC, POLE_LOCATOR
using SatelliteDynamics.Time: Epoch, mjd

##############
# RTN | LVLH #
##############

export rRTNtoECI
"""
Compute the radial, along-track, cross-track (RTN) rotation matrix. Which,
if applied to a position vector in the RTN frame, will transform that vector to
beinto the equivalent relative position in the ECI frame.

The RTN frame is also commonly refered to as the local-vertical, local-horizontal (LVLH) frame.

Arguments:
- `x::Array{<:Real, 1}`: Inertial state (position and velocity) of primary (observing) satellite
- `xt::Array{<:Real, 1}`: Inertial state (position and velocity) of the target satellite

Returns:
- `R_rtn_to_eci::Array{<:Real, 1}`: Rotation matrix transforming _from_ the RTN frame _to_ the ECI frame.
"""
function rRTNtoECI(x::Array{<:Real, 1})
    r = x[1:3]
    v = x[4:6]

    R = r/norm(r)
    N = cross(r, v)/norm(cross(r, v))
    T = cross(N, R)

    R_rtn2eci = hcat(R, T, N)

    return R_rtn2eci
end

export rECItoRTN
"""
Compute the Earth-centered inertial to radial, along-track, cross-track (RTN) 
rotation matrix. Which, if applied to a position vector in the ECI frame, will 
transform that vector into the equivalent position vector in the RTN frame.

The RTN frame is also commonly refered to as the local-vertical, local-horizontal (LVLH) frame.

Arguments:
- `x::Array{<:Real, 1}`: Inertial state (position and velocity) of primary (observing) satellite
- `xt::Array{<:Real, 1}`: Inertial state (position and velocity) of the target satellite

Returns:
- `R_eci_to_rtn::Array{<:Real, 1}`: Rotation matrix transforming _from_ the ECI frame _to_ the RTN frame.
"""
function rECItoRTN(x::Array{<:Real, 1})
    return rRTNtoECI(x)'
end

export sECItoRTN
"""
Compute the radial, along-track, cross-track (RTN) coordinates of a target satellite in the primary satellites RTN frame.

The RTN frame is also commonly refered to as the local-vertical, local-horizontal (LVLH) frame.

Arguments:
- `x::Array{<:Real, 1}`: Inertial state (position and velocity) of primary (observing) satellite
- `xt::Array{<:Real, 1}`: Inertial state (position and velocity) of the target satellite

Returns:
- `rtn::Array{<:Real, 1}`: Position and velocity of the target relative of the observing satellite in the RTN.
"""
function sECItoRTN(x::Array{<:Real, 1}, xt::Array{<:Real, 1}; is_relative::Bool=false)
    # Create RTN rotation matrix
    R_eci2rtn = rECItoRTN(x)

    # Initialize output vector
    x_rtn = zeros(Float64, length(xt) >= 6 ? 6 : 3)

    # Transform Position
    r          = x[1:3]
    rho        = xt[1:3] - r
    x_rtn[1:3] = R_eci2rtn*rho

    # Transform velocity
    if length(xt) >= 6
        v          = x[4:6]
        f_dot      = norm(cross(r, v))/norm(r)^2
        omega      = Array{Float64, 1}([0.0, 0.0, f_dot])
        rho_dot    = xt[4:6] - v
        x_rtn[4:6] = R_eci2rtn*rho_dot - cross(omega, x_rtn[1:3])
    end

    return x_rtn
end

export sRTNtoECI
"""
Compute the Earth-center

The RTN frame is also commonly refered to as the local-vertical, local-horizontal (LVLH) frame.

Arguments:
- `x::Array{<:Real, 1}`: Inertial state (position and velocity) of primary (observing) satellite
- `xt::Array{<:Real, 1}`: Inertial state (position and velocity) of the target satellite

Returns:
- `rtn::Array{<:Real, 1}`: Position and velocity of the target relative of the observing satellite in the RTN.
"""
function sRTNtoECI(x::Array{<:Real, 1}, xrtn::Array{<:Real, 1})
    # Create RTN rotation matrix
    R_rtn2eci = rRTNtoECI(x)

    # Initialize output vector
    xt = zeros(Float64, length(xrtn) >= 6 ? 6 : 3)

    # Transform position
    r       = x[1:3]
    r_rtn   = xrtn[1:3]
    xt[1:3] = R_rtn2eci*r_rtn + r

    # Transform velocity
    if length(xrtn) >= 6
        v = x[4:6]
        v_rtn   = xrtn[4:6]
        f_dot   = norm(cross(r, v))/norm(r)^2
        omega   = Array{Float64, 1}([0.0, 0.0, f_dot])
        xt[4:6] = R_rtn2eci*(v_rtn + cross(omega, r_rtn)) + v
    end

    return xt
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
- `rc2i::Array{<:Real, 2}`: 3x3 Rotation matrix transforming GCRS -> CIRS
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
- `r::Array{<:Real, 2}`: 3x3 Rotation matrix transforming CIRS -> TIRS
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
- `rpm::Array{<:Real, 2}`: 3x3 Rotation matrix transforming TIRS -> ITRF
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
- `r::Array{<:Real, 2}`: 3x3 Rotation matrix transforming GCRF -> ITRF
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
- `r::Array{<:Real, 1}`: 3x3 Rotation matrix transforming ITRF -> GCRF
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
- `x::Array{<:Real, 1}`: Inertial state (position, velocity) [m; m/s]

# Returns
- `x_ecef::Array{<:Real, 1}`: Earth-fixed state (position, velocity)
"""
function sECItoECEF(epc::Epoch, x::Array{<:Real, 1})
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
- `x::Array{<:Real, 1}`: Earth-fixed state (position, velocity) [m; m/s]

# Returns
- `x_ecef::Array{<:Real, 1}`: Inertial state (position, velocity)
"""
function sECEFtoECI(epc::Epoch, x::Array{<:Real, 1})
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