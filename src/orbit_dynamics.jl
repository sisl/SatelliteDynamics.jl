__precompile__(true)
module OrbitDynamics

using SatelliteDynamics.Constants
using SatelliteDynamics.Universe
using SatelliteDynamics.Time: Epoch, mjd
using SatelliteDynamics.Coordinates: sECEFtoGEOD
using SatelliteDynamics.Attitude: Rx, Ry, Rz
using SatelliteDynamics.ReferenceSystems
using LinearAlgebra

#####################
# Utility Functions #
#####################

"""
Kronecker Delta Function.

Returns 1 if inputs are equal returns 0 otherwise.

Arguments:
- `a::Real`: First input argument
- `b::Real`: Second input argument 

Returns:
- `delta:Integer`: Kronecker delta result.
"""
function kron(a::Real, b::Real)
    return a == b ? 1 : 0
end

"""
Internal helper function to aid in denormalization of gravity field coefficients.

Provides method co compute the factorial ratio: (n-m)!/(n+m)!

Arguments:
- `n::Real`: Gravity model degree, n.
- `m::Real`: Gravity model order, m.

Returns:
- `p::Real`: Factorial product
"""
function facprod(n::Int, m::Int)
    p = 1.0

    for i in (n-m+1):(n+m+1) # Is this correct? Should it be (n-m+1):(n+m)
        p = p/i
    end

    return p
end

###########
# Gravity #
###########

export accel_point_mass
"""
Computes the acceleration of a satellite caused by a point-mass approximation 
of the central body. Returns the acceleration vector of the satellite.

Assumes the satellite is much, much less massive than the central body.

Arguments:
- `r_sat::Array{<:Real, 1}`: satellite position in a commonn inertial frame [m]
- `r_body::Array{<:Real, 1}`: position of body in a commonn inertial frame [m]
- `GM::Array{<:Real, 1}`: gravitational coeffient of attracting body [m^3/s^2] Default: SatelliteDynamics.Constants.GM_EARTH)
(Default: SatelliteDynamics.Constants.GM_EARTH

Return:
- `a::Array{<:Real, 1}`: Acceleration in X, Y, and Z inertial directions [m/s^2]
"""
function accel_point_mass(r_sat::Array{<:Real, 1}, r_body::Array{<:Real, 1}, gm_body::Real=GM_EARTH)
    # Restrict inputs to position only
    r_sat  = r_sat[1:3]
    r_body = r_body[1:3]

    # Relative position vector of satellite w.r.t. the attraching body
    d = r_sat - r_body

    # Acceleration
    a = -gm_body * (d/norm(d)^3 + r_body/norm(r_body)^3)

    return a
end

"""
Computes the acceleration on a satellite caused by a point-mass approximation 
of a massive body. Returns the acceleration vector of the satellite.

Arguments:
- `r_sat::Array{<:Real, 1}`: satellite position in the inertial frame [m]
- `GM::Array{<:Real, 1}`: gravitational coeffient of attracting body [m^3/s^2] Default: SatelliteDynamics.Constants.GM_EARTH)
(Default: SatelliteDynamics.Constants.GM_EARTH

Return:
- `a::Array{<:Real, 1}`: Acceleration in X, Y, and Z inertial directions [m/s^2]
"""
function accel_point_mass(x::Array{<:Real, 1}, gm_body::Real=GM_EARTH)
    # Restrict inputs to position only. Considered in body frame
    r  = x[1:3]

    # Acceleration
    a = -gm_body * r/norm(r)^3

    return a
end

"""
Compute the gravitational acceleration at a model given a spherical harmonic
gravity field model.

Arguments:
- `r::Array{<:Real, 1}`: Position of the point in the body (field) fixed frame. [m]
- `coef::Array{<:Real, 2}`: Gravity coefficients stored in dense matrix form. C_nm terns are stored along rows indexed from C_00 in coef[1, 1] to C_nm in coef[n+1, m+1]. S_nm terms are stored along matrix columns with S_n0 stored in coef[0, n+1]
- `n_max::Integer`: Maximum degree coefficient to use in expansion
- `m_max::Integer`: Maximum order coefficient to use in the expansion. Must be less than the degree.
- `r_ref::Real`: Reference distance of the gravity field.
- `GM::Real`: Gravitational constant of central body
- `normralized::Bool`: Whether the input gravity field coefficients are normalized coefficients (Default: true)

Returns:
- `a::Array{<:Real, 1}`: Acceleration in X, Y, and Z inertial directions [m/s^2]
"""
function spherical_harmonic_gravity(r::Array{<:Real, 1}, coef::Array{<:Real, 2}, n_max::Integer, m_max::Integer, r_ref::Real, GM::Real; normalized::Bool=true)
    # Intermediate computations
    r_sqr = dot(r, r)
    rho   = r_ref^2/r_sqr
    x0    = r_ref * r[1] / r_sqr
    y0    = r_ref * r[2] / r_sqr
    z0    = r_ref * r[3] / r_sqr

    # Initialize Intermediary Matrices
    V = zeros(Float64, n_max+2, n_max+2)
    W = zeros(Float64, n_max+2, n_max+2)

    # Calculate zonal terms V(n, 0). Set W(n,0)=0.0
    V[0+1, 0+1] = r_ref /sqrt(r_sqr)
    W[0+1, 0+1] = 0.0

    V[1+1, 0+1] = z0 * V[0+1, 0+1]
    W[1+1, 0+1] = 0.0

    for n in 2:(n_max+1)
        V[n+1, 0+1] = ((2*n-1)*z0*V[n+1-1, 0+1] - (n-1)*rho*V[n+1-2,0+1])/n
        W[n+1, 0+1] = 0.0
    end

    # Calculate tesseral and sectoral terms
    for m in 1:m_max+1
        # Calculate V(m,m) to V(n_max+1,m)
        V[m+1, m+1] = (2*m-1)*(x0*V[m+1-1, m+1-1] - y0*W[m+1-1, m+1-1])
        W[m+1, m+1] = (2*m-1)*(x0*W[m+1-1, m+1-1] + y0*V[m+1-1, m+1-1])

        if m <= m_max
            V[m+1+1, m+1] = (2*m+1)*z0*V[m+1,m+1]
            W[m+1+1, m+1] = (2*m+1)*z0*W[m+1,m+1]
        end

        for n in (m+2):(n_max+1)
            V[n+1,m+1] = ((2*n-1)*z0*V[n+1-1,m+1]-(n+m-1)*rho*V[n+1-2,m+1])/(n-m)
            W[n+1,m+1] = ((2*n-1)*z0*W[n+1-1,m+1]-(n+m-1)*rho*W[n+1-2,m+1])/(n-m)
        end
    end

    # Calculate accelerations
    ax = 0.0
    ay = 0.0
    az = 0.0

    for m in 0:m_max
        for n in m:n_max
            C = 0.0
            S = 0.0
            if m == 0
                # Denormalize Coeeficients
                if normalized
                    N = sqrt(2*n+1)
                    C = N*coef[n+1, 0+1]
                else
                    C = coef[n+1, 0+1]
                end

                ax -= C*V[n+1+1, 1+1]
                ay -= C*W[n+1+1, 1+1]
                az -= (n+1)*C*V[n+1+1, 0+1]

            else
                if normalized
                    N = sqrt((2 - kron(0,m))*(2*n+1)*facprod(n,m))
                    C = N*coef[n+1,   m+1]
                    S = N*coef[m+1-1, n+1]
                else
                    C = coef[n+1,   m+1]
                    S = coef[m+1-1, n+1]
                end

                fac =  0.5 * (n-m+1)*(n-m+2)
                ax += +0.5 * (-C * V[n+1+1, m+1+1] - S * W[n+1+1, m+1+1])
                      +fac * (+C * V[n+1+1, m+1-1] + S * W[n+1+1, m+1-1])
                ay += +0.5 * (-C * W[n+1+1, m+1+1] + S * V[n+1+1, m+1+1])
                      +fac * (-C * W[n+1+1, m+1-1] + S * V[n+1+1, m+1-1])
                az += (n-m+1)*(-C*V[n+1+1, m+1] - S*W[n+1+1, m+1])
            end
        end
    end

    a = (GM / (r_ref^2)) * [ax, ay, az]

    return a
end

export accel_gravity
"""
Computes the accleration caused by Earth gravity as modeled by a spherical 
harmonic gravity field.

Arguments:
- `r_sat::Array{<:Real, 1}`: Satellite position in the inertial frame [m]
- `R_eci_ecef::Array{<:Real, 2}`: Rotation matrix transforming a vector from the inertial to body-fixed reference frames. 
- `n_max::Integer`: Maximum degree coefficient to use in expansion
- `m_max::Integer`: Maximum order coefficient to use in the expansion. Must be less than the degree.
    
Return:
- `a::Array{<:Real, 1}`: Gravitational acceleration in X, Y, and Z inertial directions [m/s^2]

References:
1. O. Montenbruck, and E. Gill, _Satellite Orbits: Models, Methods and Applications_, 2012, p.56-68.
"""
function accel_gravity(x::Array{<:Real, 1}, R_eci_ecef::Array{<:Real, 2}, n_max::Int=20, m_max::Int=20)
    
    # Check Limits of Gravity Field
    if n_max > GRAVITY_MODEL.n_max
        error("Requested gravity model order $n_max is larger than the maximum order of the model maximum: ", GRAVITY_MODEL.n_max)
    end

    if m_max > GRAVITY_MODEL.m_max
        error("Requested gravity model degree $m_max is larger than the maximum degree of the model maximum: ", GRAVITY_MODEL.m_max)
    end

    # Gravitational parameters of primary body
    GM    = GRAVITY_MODEL.GM
    R_ref = GRAVITY_MODEL.R

    # Body-fixed position
    r_bf = R_eci_ecef * x[1:3]

    # Compute spherical harmonic acceleration
    a_ecef = spherical_harmonic_gravity(r_bf, GRAVITY_MODEL.data, n_max, m_max, R_ref, GM, normalized=GRAVITY_MODEL.normalized)

    # Inertial acceleration
    a_eci = R_eci_ecef' * a_ecef

    # Finished
    return a_eci
end

#########################
# Planetary Ephemerides #
#########################

export sun_position
"""
Compute the Sun's position in the EME2000 inertial frame through the use
of low-precision analytical functions.

Argument:
- `epc::Epoch`: Epoch

Returns:
- `r_sun::Array{<:Real, 1}`: Position vector of the Sun in the Earth-centered inertial fame.

Notes:
1. The EME2000 inertial frame is for most purposes equivalent to the GCRF frame.

References:
1. O. Montenbruck, and E. Gill, _Satellite Orbits: Models, Methods and Applications_, 2012, p.70-73.
"""
function sun_position(epc::Epoch)
    # Constants
    mjd_tt  = mjd(epc, tsys=:TT)      # MJD of epoch in TT
    epsilon = 23.43929111*pi/180.0     # Obliquity of J2000 ecliptic
    T       = (mjd_tt-MJD2000)/36525.0 # Julian cent. since J2000

    # Variables

    # Mean anomaly, ecliptic longitude and radius
    M = 2.0*pi * modf(0.9931267 + 99.9973583*T)[1]                 # [rad]
    L = 2.0*pi * modf(0.7859444 + M/(2.0*pi) + (6892.0*sin(M)+72.0*sin(2.0*M)) / 1296.0e3)[1] # [rad]
    r = 149.619e9 - 2.499e9*cos(M) - 0.021e9*cos(2*M)           # [m]

    # Equatorial position vector
    p_sun = Rx(-epsilon) * [r*cos(L), r*sin(L), 0.0]

    return p_sun
end

export moon_position
"""
Compute the Moon's position in the EME2000 inertial frame through the use
of low-precision analytical functions.

Argument:
- `epc::Epoch`: Epoch

Returns:
- `r_moon::Array{<:Real, 1}`: Position vector of the Moon in the Earth-centered inertial fame.

Notes:
1. The EME2000 inertial frame is for most purposes equivalent to the GCRF frame.

References:
1. O. Montenbruck, and E. Gill, _Satellite Orbits: Models, Methods and Applications_, 2012, p.70-73.
"""
function moon_position(epc::Epoch)
    # Constants
    mjd_tt  = mjd(epc, tsys=:TT)      # MJD of epoch in TT
    epsilon = 23.43929111*pi/180.0     # Obliquity of J2000 ecliptic
    T       = (mjd_tt-MJD2000)/36525.0 # Julian cent. since J2000

    # Mean elements of lunar orbit
    L_0 =     modf(0.606433 + 1336.851344*T)[1] # Mean longitude [rev] w.r.t. J2000 equinox
    l   = 2.0*pi*modf(0.374897 + 1325.552410*T)[1] # Moon's mean anomaly [rad]
    lp  = 2.0*pi*modf(0.993133 +   99.997361*T)[1] # Sun's mean anomaly [rad]
    D   = 2.0*pi*modf(0.827361 + 1236.853086*T)[1] # Diff. long. Moon-Sun [rad]
    F   = 2.0*pi*modf(0.259086 + 1342.227825*T)[1] # Argument of latitude 


    # Ecliptic longitude (w.r.t. equinox of J2000)
    dL = + 22640*sin(l) - 4586*sin(l-2*D) + 2370*sin(2*D) +  769*sin(2*l)
         - 668*sin(lp) - 412*sin(2*F) - 212*sin(2*l-2*D) - 206*sin(l+lp-2*D)
         + 192*sin(l+2*D) - 165*sin(lp-2*D) - 125*sin(D) - 110*sin(l+lp)
         + 148*sin(l-lp) - 55*sin(2*F-2*D)

    L = 2.0*pi * modf(L_0 + dL/1296.0e3)[1]  # [rad]

    # Ecliptic latitude
    S  = F + (dL+412*sin(2*F)+541*sin(lp)) * AS2RAD 
    h  = F-2*D
    N  = - 526*sin(h) + 44*sin(l+h) - 31*sin(-l+h) - 23*sin(lp+h)
         + 11*sin(-lp+h) - 25*sin(-2*l+F) + 21*sin(-l+F)
    B  = (18520.0*sin(S) + N) * AS2RAD   # [rad]

    # Distance [m]
    r = + 385000e3 - 20905e3*cos(l) - 3699e3*cos(2*D-l) - 2956e3*cos(2*D)
        - 570e3*cos(2*l) + 246e3*cos(2*l-2*D) - 205e3*cos(lp-2*D)
        - 171e3*cos(l+2*D) - 152e3*cos(l+lp-2*D)   

    # Equatorial coordinates
    p_moon = Rx(-epsilon) * [r*cos(L)*cos(B), r*sin(L)*cos(B), r*sin(B)]

    return p_moon
end

######################
# Third-Body Gravity #
######################

export accel_thirdbody_sun
"""
Computes the acceleration of a satellite in the inertial frame due to the
gravitational attraction of the Sun.

Arguments:
- `x::Array{<:Real, 1}`: Satellite Cartesean state in the inertial reference frame [m; m/s]
- `r_sun::Array{<:Real, 1}`: Position of sun in inertial frame.

Return:
- `a::Array{<:Real, 1}`: Acceleration due to the Sun's gravity in the inertial frame [m/s^2]

References:
1. O. Montenbruck, and E. Gill, _Satellite Orbits: Models, Methods and Applications_, 2012, p.69-70.
"""
function accel_thirdbody_sun(x::Array{<:Real, 1}, r_sun::Array{<:Real, 1})
    # Acceleration due to sun point mass
    a_sun = accel_point_mass(x[1:3], r_sun, GM_SUN)

    return a_sun
end

function accel_thirdbody_sun(epc::Epoch, x::Array{<:Real, 1})
    # Compute solar position
    r_sun = sun_position(epc)

    # Acceleration due to sun point mass
    a_sun = accel_point_mass(x[1:3], r_sun, GM_SUN)

    return a_sun
end

export accel_thirdbody_moon
"""
Computes the acceleration of a satellite in the inertial frame due to the
gravitational attraction of the Moon.

Arguments:
- `x::Array{<:Real, 1}`: Satellite Cartesean state in the inertial reference frame [m; m/s]
- `r_moon::Array{<:Real, 1}`: Position of moon in inertial frame.

Returns:
- `a::Array{<:Real, 1}`: Acceleration due to the Moon's gravity in the inertial frame [m/s^2]

References:
1. O. Montenbruck, and E. Gill, _Satellite Orbits: Models, Methods and Applications_, 2012, p.69-70.
"""
function accel_thirdbody_moon(x::Array{<:Real, 1}, r_moon::Array{<:Real, 1})
    # Acceleration due to moon point mass
    a_moon = accel_point_mass(x[1:3], r_moon, GM_MOON)

    return a_moon
end

function accel_thirdbody_moon(epc::Epoch, x::Array{<:Real, 1})
    # Compute solar position
    r_moon = moon_position(epc)

    # Acceleration due to moon point mass
    a_moon = accel_point_mass(x[1:3], r_moon, GM_MOON)

    return a_moon
end

####################
# Atmospheric Drag #
####################

export density_harris_priester
"""
Computes the local density using the Harris-Priester density model.

Arguments:
- `x::Array{<:Real, 1}`: Satellite Cartesean state in the inertial reference frame [m; m/s]
- `r_sun::Array{<:Real, 1}`: Position of sun in inertial frame.

Returns:
- `rho:Float64`: Local atmospheric density [kg/m^3]

References:
1. O. Montenbruck, and E. Gill, _Satellite Orbits: Models, Methods and Applications_, 2012, p.89-91.
"""
function density_harris_priester(x::Array{<:Real, 1}, r_sun::Array{<:Real, 1})
    # Harris-Priester Constants
    hp_upper_limit =   1000.0          # Upper height limit [km]
    hp_lower_limit =    100.0          # Lower height limit [km]
    hp_ra_lag      = 0.523599          # Right ascension lag [rad]
    hp_n_prm       =        3          # Harris-Priester parameter 
                                        # 2(6) low(high) inclination
    hp_N           = 50                # Number of coefficients

    # Height [km]
    hp_h = [100.0, 120.0, 130.0, 140.0, 150.0, 160.0, 170.0, 180.0, 190.0, 200.0,     
            210.0, 220.0, 230.0, 240.0, 250.0, 260.0, 270.0, 280.0, 290.0, 300.0,     
            320.0, 340.0, 360.0, 380.0, 400.0, 420.0, 440.0, 460.0, 480.0, 500.0,     
            520.0, 540.0, 560.0, 580.0, 600.0, 620.0, 640.0, 660.0, 680.0, 700.0,     
            720.0, 740.0, 760.0, 780.0, 800.0, 840.0, 880.0, 920.0, 960.0,1000.0]

    # Minimum density [g/km^3]
    hp_c_min = [4.974e+05, 2.490e+04, 8.377e+03, 3.899e+03, 2.122e+03, 1.263e+03,         
                8.008e+02, 5.283e+02, 3.617e+02, 2.557e+02, 1.839e+02, 1.341e+02,         
                9.949e+01, 7.488e+01, 5.709e+01, 4.403e+01, 3.430e+01, 2.697e+01,         
                2.139e+01, 1.708e+01, 1.099e+01, 7.214e+00, 4.824e+00, 3.274e+00,         
                2.249e+00, 1.558e+00, 1.091e+00, 7.701e-01, 5.474e-01, 3.916e-01,         
                2.819e-01, 2.042e-01, 1.488e-01, 1.092e-01, 8.070e-02, 6.012e-02,         
                4.519e-02, 3.430e-02, 2.632e-02, 2.043e-02, 1.607e-02, 1.281e-02,         
                1.036e-02, 8.496e-03, 7.069e-03, 4.680e-03, 3.200e-03, 2.210e-03,         
                1.560e-03, 1.150e-03]

    # Maximum density [g/km^3]
    hp_c_max = [4.974e+05, 2.490e+04, 8.710e+03, 4.059e+03, 2.215e+03, 1.344e+03,         
                8.758e+02, 6.010e+02, 4.297e+02, 3.162e+02, 2.396e+02, 1.853e+02,         
                1.455e+02, 1.157e+02, 9.308e+01, 7.555e+01, 6.182e+01, 5.095e+01,         
                4.226e+01, 3.526e+01, 2.511e+01, 1.819e+01, 1.337e+01, 9.955e+00,         
                7.492e+00, 5.684e+00, 4.355e+00, 3.362e+00, 2.612e+00, 2.042e+00,         
                1.605e+00, 1.267e+00, 1.005e+00, 7.997e-01, 6.390e-01, 5.123e-01,         
                4.121e-01, 3.325e-01, 2.691e-01, 2.185e-01, 1.779e-01, 1.452e-01,         
                1.190e-01, 9.776e-02, 8.059e-02, 5.741e-02, 4.210e-02, 3.130e-02,         
                2.360e-02, 1.810e-02]

    # Satellite height
    geod   = sECEFtoGEOD(x[1:3], use_degrees=true)
    height = geod[3]/1.0e3 # height in [km]

    # Exit with zero density outside height model limits
    if height > hp_upper_limit || height < hp_lower_limit
        return 0.0
    end


    # Sun right ascension, declination
    ra_sun  = atan( r_sun[2], r_sun[1] )
    dec_sun = atan( r_sun[3], sqrt( r_sun[1]^2 + r_sun[2]^2 ) )


    # Unit vector u towards the apex of the diurnal bulge
    # in inertial geocentric coordinates
    c_dec = cos(dec_sun)
    u     = [c_dec * cos(ra_sun + hp_ra_lag),
             c_dec * sin(ra_sun + hp_ra_lag),
             sin(dec_sun)]


    # Cosine of half angle between satellite position vector and
    # apex of diurnal bulge
    c_psi2 = 0.5 + 0.5 * dot(x[1:3], u)/norm(x[1:3])

    # Height index search and exponential density interpolation
    ih = 0                            # section index reset
    for i in 1:hp_N                   # loop over N_Coef height regimes
        if height >= hp_h[i] && height < hp_h[i+1] 
            ih = i                    # ih identifies height section
            break
        end
    end

    h_min = ( hp_h[ih] - hp_h[ih+1] )/log( hp_c_min[ih+1]/hp_c_min[ih] )
    h_max = ( hp_h[ih] - hp_h[ih+1] )/log( hp_c_max[ih+1]/hp_c_max[ih] )

    d_min = hp_c_min[ih] * exp( (hp_h[ih]-height)/h_min )
    d_max = hp_c_max[ih] * exp( (hp_h[ih]-height)/h_max )

    # Density computation
    density = d_min + (d_max-d_min) * c_psi2^hp_n_prm

    # Convert from g/km^3 to kg/m^3
    density *= 1.0e-12

    # Finished
    return density
end

function density_harris_priester(epc::Epoch, x::Array{<:Real, 1})
    r_sun = sun_position(epc)
    return density_harris_priester(x, r_sun)
end

export accel_drag
"""
Computes the perturbing, non-conservative acceleration caused by atmospheric
drag assuming that the ballistic properties of the spacecraft are captured by
the coefficient of drag.

Arguments:
- `x::Array{<:Real, 1}`: Satellite Cartesean state in the inertial reference frame [m; m/s]
- `rho::Real`: atmospheric density [kg/m^3]
- `mass::Real`: Spacecraft mass [kg]
- `area::Real`: Wind-facing cross-sectional area [m^2]
- `Cd::Real`: coefficient of drag [dimensionless]
- `T::Array{<:Real, 2}`: Rotation matrix from the inertial to the true-of-date frame

Return:
- `a::Array{<:Real, 1}`: Acceleration due to drag in the X, Y, and Z inertial directions. [m/s^2]

References:
1. O. Montenbruck, and E. Gill, _Satellite Orbits: Models, Methods and Applications_, 2012, p.83-86.
"""
function accel_drag(x::Array{<:Real, 1}, rho::Real, mass::Real, area::Real, Cd::Real, T::Array{<:Real, 2})

    # Constants
    omega = [0, 0, OMEGA_EARTH]

    # Position and velocity in true-of-date system
    r_tod = T * x[1:3]
    v_tod = T * x[4:6]

    # Velocity relative to the Earth's atmosphere
    v_rel = v_tod - cross(omega, r_tod)
    v_abs = norm(v_rel)

    # Acceleration 
    a_tod  = -0.5*Cd*(area/mass)*rho*v_abs*v_rel
    a_drag = T' * a_tod 

    return a_drag
end

############################
# Solar Radiation Pressure #
############################

export eclipse_cylindrical
"""
Computes the illumination fraction of a satellite in Earth orbit using a
cylindrical Earth shadow model.

Arguments:
- `x::Array{<:Real, 1}`: Satellite Cartesean state in the inertial reference frame [m; m/s]
- `r_sun::Array{<:Real, 1}`: Position of sun in inertial frame.

Return:
- `nu::Float64`: Illumination fraction (0 <= nu <= 1). nu = 0 means spacecraft in complete shadow, nu = 1 mean spacecraft fully illuminated by sun.

References:
1. O. Montenbruck, and E. Gill, _Satellite Orbits: Models, Methods and Applications_, 2012, p.80-83.
"""
function eclipse_cylindrical(x::Array{<:Real, 1}, r_sun::Array{<:Real, 1})
    # Satellite inertial position
    r = x[1:3]
    
    # Sun-direction unit-vector
    e_sun = r_sun / norm(r_sun)
    
    # Projection of spacecraft position
    s = dot(r, e_sun)

    # Compute illumination
    nu = 0.0
    if s/norm(s) >= 1.0 || norm(r - s*e_sun) > R_EARTH 
        nu = 1.0
    end

    return nu
end

function eclipse_cylindrical(epc::Epoch, x::Array{<:Real, 1})
    r_sun = sun_position(epc)
    return eclipse_cylindrical(x, r_sun)
end

export eclipse_conical
"""
Computes the illumination fraction of a satellite in Earth orbit using a
conical Earth shadow model.

Arguments:
- `x::Array{<:Real, 1}`: Satellite Cartesean state in the inertial reference frame [m; m/s]
- `r_sun::Array{<:Real, 1}`: Position of sun in inertial frame.

Return:
- `nu::Float64`: Illumination fraction (0 <= nu <= 1). nu = 0 means spacecraft in complete shadow, nu = 1 mean spacecraft fully illuminated by sun.

References:
1. O. Montenbruck, and E. Gill, _Satellite Orbits: Models, Methods and Applications_, 2012, p.80-83.
"""
function eclipse_conical(x::Array{<:Real, 1}, r_sun::Array{<:Real, 1})
    # Satellite inertial position
    r = x[1:3]

    # Occultation Geometry
    a = asin(R_SUN/norm(r_sun - r))
    b = asin(R_EARTH/norm(r))
    c = acos(dot(r, r_sun-r)/(norm(r)*norm(r_sun-r)) )

    e_sun = r_sun / norm(r_sun)

    # Test Occulation Conditions
    nu = 0.0
    if abs(a - b) < c && c < (a + b)
        # Partial occultation
    
        xx = (c^2 + a^2 - b^2)/(2*c)
        yy = sqrt(a^2 - xx^2)
        A  = a^2 * acos(xx/a) + b^2 * acos((c-xx)/b) - c * yy

        nu = 1 - A/(pi*a^2)
    elseif (a + b) <= c
        # No occultation
        nu = 1.0
    else
        # Full occultation
        nu = 0.0
    end

    return nu
end


function eclipse_conical(epc::Epoch, x::Array{<:Real, 1})
    r_sun = sun_position(epc)
    return eclipse_conical(x, r_sun)
end

export accel_srp
"""Computes the perturbing acceleration due to direct solar radiation 
pressure assuming the reflecting surface is a flat plate pointed directly at
the Sun.

Arguments:
- `x::Array{<:Real, 1}`: Satellite Cartesean state in the inertial reference frame [m; m/s]

Returns:
- `a::Array{<:Real, 1}`: Satellite acceleration due to solar radiation pressure [m/s^2]

References:
1. O. Montenbruck, and E. Gill, _Satellite Orbits: Models, Methods and Applications_, 2012, p.77-79.
"""
function accel_srp(x::Array{<:Real, 1}, r_sun::Array{<:Real, 1}, mass::Real=0, area::Real=0, CR::Real=1.8, p0::Real=P_SUN, au::Real=AU)
    # Spacecraft position vector
    r = x[1:3]

    # Relative position vector of spacecraft w.r.t. Sun
    d = r - r_sun

    # Acceleration due to moon point mass
    a_srp = d * (CR*(area/mass)*p0*AU^2 / norm(d)^3)

    # Return
    return a_srp
end

##############
# Relativity #
##############

export accel_relativity
"""Computes perturbation accleration of a satellite in the Inertial frame
due to the combined effects of special and general relativity.

Arguments:
- `x::Array{<:Real, 1}`: Satellite Cartesean state in the inertial reference frame [m; m/s]

Returns:
- `a::Array{<:Real, 1}`: Satellite acceleration due to relativity. [m/s^2]

References:
1. O. Montenbruck, and E. Gill, _Satellite Orbits: Models, Methods and Applications_, 2012, p.110-112.
"""
function accel_relativity(x::Array{<:Real, 1})
    # Extract state variables
    r = x[1:3]
    v = x[4:6]

    # Intermediate computations
    norm_r = norm(r)
    r2     = norm_r^2

    norm_v = norm(v)
    v2     = norm_v^2

    c  = C_LIGHT
    c2 = c^2

    # Compute unit vectors
    er = r/norm_r
    ev = v/norm_v

    # Compute perturbation and return
    a_rel = GM_EARTH/r2 * ( (4*GM_EARTH/(c2*norm_r) - v2/c2)*er + 4*v2/c2*dot(er, ev)*ev)

    return a_rel
end

end