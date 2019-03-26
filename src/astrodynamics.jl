__precompile__(true)
module Astrodynamics


using LinearAlgebra
using SatelliteDynamics.Constants: GM_EARTH, R_EARTH, J2_EARTH


export mean_motion
"""
Compute the mean motion given a semi-major axis.

Arguments:
- `a::Real`: Semi-major axis. [m]
- `use_degrees:Bool`: If `true` returns result in units of degrees
- `GM::Real`: Gravitational constant of central body. Defaults to `SatelliteDynamics.Constants.GM_EARTH` if none is provided.

Returns:
- `n::Real`: Orbital mean motion. [rad/s] or [deg/s]
"""
function mean_motion(a::Real; use_degrees::Bool=false, GM::Real=GM_EARTH)
    n = sqrt(GM/a^3)

    if use_degrees
        n *= 180.0/pi
    end

    return n
end

export semimajor_axis
"""
Calculate semi-major axis given mean-motion.

Arguments:
- `n::Real`: Orbital mean motion. [rad/s] or [deg/s]
- `use_degrees:Bool`: If `true` interpret input as being in degrees.
- `GM::Real`: Gravitational constant of central body. Defaults to `SatelliteDynamics.Constants.GM_EARTH` if none is provided.

Returns:
- `a::Real`: Semi-major axis. [m]
"""
function semimajor_axis(n::Real; use_degrees::Bool=false, GM::Real=GM_EARTH)
    if use_degrees
        n *= pi/180.0
    end

    a = (GM/n^2)^(1.0/3.0)

    return a
end

export orbit_period
"""
Compute the satellite orbital period given the semi-major axis.

Arguments:
- `a::Real`: Semi-major axis. [m]
- `GM::Real`: Gravitational constant of central body. Defaults to `SatelliteDynamics.Constants.GM_EARTH` if none is provided.

Returns:
- `T::Real`: Orbital period. [s]
"""
function orbit_period(a::Real; GM::Real=GM_EARTH)
    return 2.0*pi*sqrt(a^3/GM)
end

export sun_sync_inclination
"""
Compute the required inclination for a Sun-synchronous Earth orbit.

Algorithm assumes the nodal precession is entirely due to the J2 perturbation, and no other perturbations are considered.

The inclination is computed using a first-order, non-iterative approximation.

Arguments:
- `a::Real`: Semi-major axis. [m]
- `e::Real`: Eccentricity. [dimensionless]
- `use_degrees:Bool`: If `true` interpret output will be returned in degrees.

Returns:
- `iss::Real`: Requierd inclination for a sun-synchronous orbit. [rad] or [deg]
"""
function sun_sync_inclination(a::Real, e::Real; use_degrees::Bool=false)
    # Compute the required RAAN precession of a sun-synchronous orbit
    OMEGA_DOT_SS = 2*pi/365.2421897/86400.0

    # Inclination required for sun-synchronous orbits
    i_ss = acos(-2*a^(7/2)*OMEGA_DOT_SS*(1-e^2)^2/(3*(R_EARTH^2)*J2_EARTH*sqrt(GM_EARTH)))

    if use_degrees == true
        i_ss *= 180.0/pi
    end

    return i_ss
end

export anomaly_eccentric_to_mean
"""
Convert eccentric anomaly into mean anomaly.

Arguments:
- `E::Real`: Eccentric anomaly. [rad] or [deg]
- `e::Real`: Eccentricity. [dimensionless]
- `use_degrees:Bool`: If `true` interpret input will be interpreted as being in degrees, and output will be returned in degrees.

Returns:
- `M::Real`: Mean anomaly. [rad] or [deg]
"""
function anomaly_eccentric_to_mean(E::Real, e::Real; use_degrees::Bool=false)
    # Convert degree input
    if use_degrees == true
        E *= pi/180.0
    end

    # Convert eccentric to mean
    M = E - e*sin(E)

    # Convert degree output
    if use_degrees == true
        M *= 180.0/pi
    end

    return M
end

export anomaly_mean_to_eccentric
"""
Convert mean anomaly into eccentric anomaly.

Arguments:
- `M::Real`: Mean anomaly. [deg] or [deg]
- `e::Real`: Eccentricity. [dimensionless]
- `use_degrees:Bool`: If `true` interpret input will be interpreted as being in degrees, and output will be returned in degrees.

Returns:
- `E::Real`: Eccentric anomaly. [rad] or [deg]
"""
function anomaly_mean_to_eccentric(M::Real, e::Real; use_degrees::Bool=false)
    # Convert degree input
    if use_degrees == true
        M *= pi/180.0
    end

    # Convert mean to eccentric
    max_iter = 15
    epsilson = eps(Float64)*100.0

    # Initialize starting values
    M = mod(M, 2.0*pi)
    if e < 0.8
        E = M
    else
        E = pi
    end

    # Initialize working variable
    f = E - e*sin(E) - M
    i = 0

    # Iterate until convergence
    while abs(f) > epsilson
        f = E - e*sin(E) - M
        E = E - f / (1.0 - e*cos(E))

        # Increase iteration counter
        i += 1
        if i == max_iter
            error("Maximum number of iterations reached before convergence.")
        end
    end

    # Convert degree output
    if use_degrees == true
        E *= 180.0/pi
    end

    return E
end

export sOSCtoCART
"""
Given an orbital state expressed in osculating orbital elements compute the equivalent Cartesean position and velocity of the inertial state.

The osculating elements are assumed to be (in order):
1. _a_, Semi-major axis [m]
2. _e_, Eccentricity [dimensionless]
3. _i_, Inclination [rad]
4. _Ω_, Right Ascension of the Ascending Node (RAAN) [rad]
5. _ω_, Argument of Perigee [ramd]
6. _M_, Mean anomaly [rad]

Arguments:
- x_oe `x::Array{<:Real, 1}`: Osculating orbital elements. See above for desription of the elements and their required order.
- `use_degrees:Bool`: If `true` interpret input will be interpreted as being in degrees, and output will be returned in degrees.
- `GM::Real`: Gravitational constant of central body. Defaults to `SatelliteDynamics.Constants.GM_EARTH` if none is provided.

# Returns
- x `x::Array{<:Real, 1}`: Cartesean inertial state. Returns position and velocity. [m; m/s]
"""
function sOSCtoCART(x_oe::Array{<:Real, 1}; use_degrees::Bool=false, GM::Real=GM_EARTH)

    if use_degrees == true
        # Copy and convert input from degrees to radians if necessary
        oe = deepcopy(x_oe)
        oe[3:6] = oe[3:6]*pi/180.0
    else
        oe = x_oe
    end
    
    # Unpack input
    a, e, i, RAAN, omega, M = oe

    E = anomaly_mean_to_eccentric(M, e)

    # Create perifocal coordinate vectors
    P    = zeros(Float64, 3)
    P[1] = cos(omega)*cos(RAAN) - sin(omega)*cos(i)*sin(RAAN)
    P[2] = cos(omega)*sin(RAAN) + sin(omega)*cos(i)*cos(RAAN)
    P[3] = sin(omega)*sin(i)

    Q    = zeros(Float64, 3)
    Q[1] = -sin(omega)*cos(RAAN) - cos(omega)*cos(i)*sin(RAAN)
    Q[2] = -sin(omega)*sin(RAAN) + cos(omega)*cos(i)*cos(RAAN)
    Q[3] =  cos(omega)*sin(i)

    # Find 3-Dimensional Position
    x = zeros(Float64, 6)
    x[1:3] = a*(cos(E)-e)*P + a*sqrt(1-e*e)*sin(E)*Q
    x[4:6] = sqrt(GM*a)/norm(x[1:3])*(-sin(E)*P + sqrt(1-e*e)*cos(E)*Q)

    return x
end

export sCARTtoOSC
"""
Given a Cartesean position and velocity in the inertial frame, return the 
state expressed in terms of  osculating orbital elements.

The osculating elements are assumed to be (in order):
1. _a_, Semi-major axis [m]
2. _e_, Eccentricity [dimensionless]
3. _i_, Inclination [rad]
4. _Ω_, Right Ascension of the Ascending Node (RAAN) [rad]
5. _ω_, Argument of Perigee [ramd]
6. _M_, Mean anomaly [rad]

Arguments:
- x `x::Array{<:Real, 1}`: Cartesean inertial state. Returns position and velocity. [m; m/s]
- `use_degrees:Bool`: If `true` interpret input will be interpreted as being in degrees, and output will be returned in degrees.
- `GM::Real`: Gravitational constant of central body. Defaults to `SatelliteDynamics.Constants.GM_EARTH` if none is provided.

# Returns
- x_oe `x::Array{<:Real, 1}`: Osculating orbital elements. See above for desription of the elements and their required order.
"""
function sCARTtoOSC(x::Array{<:Real, 1}; use_degrees::Bool=false, GM::Real=GM_EARTH)

    # Initialize Cartesian Polistion and Velocity
    r = x[1:3]
    v = x[4:6]

    h = cross(r, v)  # Angular momentum vector
    W = h/norm(h)    # Unit vector along angular momentum vector

    i     = atan(sqrt(W[1]*W[1] + W[2]*W[2]), W[3])    # Compute inclination
    OMEGA = atan(W[1], -W[2])                          # Right ascension of ascending node                     # Compute RAAN
    p     = norm(h)*norm(h)/GM                         # Semi-latus rectum
    a     = 1.0/(2.0/norm(r) - norm(v)*norm(v)/GM)     # Semi-major axis
    n     = sqrt(GM/(a^3))                             # Mean motion

    # Numerical stability hack for circular and near-circular orbits
    # Ensures that (1-p/a) is always positive
    if isapprox(a, p, atol=1e-9, rtol=1e-8)
        p = a
    end

    e     = sqrt(1 - p/a)                              # Eccentricity
    E     = atan(dot(r, v)/(n*a^2), (1-norm(r)/a))     # Eccentric Anomaly
    M     = anomaly_eccentric_to_mean(E, e)            # Mean Anomaly
    u     = atan(r[3], -r[1]*W[2] + r[2]*W[1])         # Mean longiude
    nu    = atan(sqrt(1-e*e)*sin(E), cos(E)-e)         # True Anomaly
    omega = u - nu                                     # Argument of perigee

    # Correct angles to run from 0 to 2PI
    OMEGA = OMEGA + 2.0*pi
    omega = omega + 2.0*pi
    M     = M     + 2.0*pi

    OMEGA = mod(OMEGA, 2.0*pi)
    omega = mod(omega, 2.0*pi)
    M     = mod(M, 2.0*pi)

    # Create Orbital Element Vector
    x_oe = [a, e, i, OMEGA, omega, M]

    # Convert output to degrees if necessary
    if use_degrees == true
        x_oe[3:6] = x_oe[3:6]*180.0/pi
    end

    return x_oe
end

end