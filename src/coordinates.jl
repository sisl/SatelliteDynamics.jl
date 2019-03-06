__precompile__(true)
module Coordinates

using LinearAlgebra

using SatelliteDynamics.Constants: WGS84_a, WGS84_f

####################
# Helper Constants #
####################

# Intermidiate calculations calculations
const ECC2 = WGS84_f * (2.0 - WGS84_f) # Square of eccentricisty

##############
# Geocentric #
##############

export sGEOCtoECEF
"""
Convert geocentric position to equivalent Earth-fixed position.

Arguments:
- `geoc::Array{<:Real, 1}`: Geocentric coordinates (lon, lat, altitude) [rad] / [deg]
- `use_degrees:Bool`: If `true` interpret input as being in degrees.

Returns:
- `ecef::Array{<:Real, 1}`: Earth-fixed coordinates [m]
"""
function sGEOCtoECEF(geoc::Array{<:Real, 1} ; use_degrees::Bool=false)
    # Extract lat and lon
    lon = geoc[1]
    lat = geoc[2]
    
    # Handle non-explict use-degrees
    if length(geoc) == 3
        alt = geoc[3]
    else
        alt = 0.0
    end

    # Convert input to radians
    if use_degrees
        lat = lat*pi/180.0
        lon = lon*pi/180.0
    end

    # Check validity of input
    if lat < -pi/2 || lat > pi/2
        throw(ArgumentError("Lattiude, $lat, out of range. Must be between -90 and 90 degrees."))
    end

    # Compute Earth fixed coordinates
    r       = WGS84_a + alt
    x = r*cos(lat)*cos(lon)
    y = r*cos(lat)*sin(lon)
    z = r*sin(lat)
    
    return [x, y, z]
end

export sECEFtoGEOC
"""
Convert Earth-fixed position to geocentric location.

Arguments:
- `ecef::Array{<:Real, 1}`: Earth-fixed coordinated [m]
- `use_degrees:Bool`: If `true` returns result in units of degrees

Returns:
- `geoc`: Geocentric coordinates (lon, lat, altitude) [rad] / [deg]
"""
function sECEFtoGEOC(ecef::Array{<:Real, 1} ; use_degrees::Bool=false)
    # Expand ECEF coordinates
    x, y, z = ecef

    # Compute geocentric coordinates
    lat = atan(z, sqrt(x*x + y*y))
    lon = atan(y, x)
    alt = sqrt(x*x + y*y + z*z) - WGS84_a

    # Convert output to degrees
    if use_degrees
        lat = lat*180.0/pi
        lon = lon*180.0/pi
    end

    return [lon, lat, alt]
end

########################
# Geodetic Convertions #
########################

export sGEODtoECEF
"""
Convert geodetic position to equivalent Earth-fixed position.

Arguments:
- `geod::Array{<:Real, 1}`: Geodetic coordinates (lon, lat, altitude) [rad] / [deg]
- `use_degrees:Bool`: If `true` interpret input as being in degrees.

Returns:
- `ecef::Array{<:Real, 1}`: Earth-fixed coordinates [m]
"""
function sGEODtoECEF(geod::Array{<:Real, 1} ; use_degrees::Bool=false)
    # Extract lat and lon
    lon = geod[1]
    lat = geod[2]
    
    # Handle non-explict use-degrees
    if length(geod) == 3
        alt = geod[3]
    else
        alt = 0.0
    end

    # Convert input to radians
    if use_degrees
        lat = lat*pi/180.0
        lon = lon*pi/180.0
    end

    # Check validity of input
    if lat < -pi/2 || lat > pi/2
        throw(ArgumentError("Lattiude, $lat, out of range. Must be between -90 and 90 degrees."))
    end

    # Compute Earth-fixed position vector
    N = WGS84_a / sqrt(1.0 - ECC2*sin(lat)^2)
    x =           (N+alt)*cos(lat)*cos(lon)
    y =           (N+alt)*cos(lat)*sin(lon)
    z =  ((1.0-ECC2)*N+alt)*sin(lat)
    
    return [x, y, z]
end

export sECEFtoGEOD
"""
Convert geodetic coordinaties to Earth-fixed position

Arguments:
- `ecef::Array{<:Real, 1}`: Earth-fixed position [m]
- `use_degrees:Bool`: If `true` returns result in units of degrees

Returns:
- `geod::Array{<:Real, 1}`: Geocentric coordinates (lon, lat, altitude) [rad] / [deg]
"""
function sECEFtoGEOD(ecef::Array{<:Real, 1} ; use_degrees::Bool=false)
    # Expand ECEF coordinates
    x, y, z = ecef


    # Compute intermediate quantities
    epsilon  = eps(Float64) * 1.0e3 * WGS84_a # Convergence requirement as function of machine precision
    rho2 = x^2 + y^2                      # Square of the distance from the z-axis
    dz   = ECC2 * z
    N    = 0.0

    # Iteratively compute refine coordinates
    while true
        zdz    = z + dz
        Nh     = sqrt(rho2 + zdz^2)
        sinphi = zdz / Nh
        N      = WGS84_a / sqrt(1.0 - ECC2 * sinphi^2)
        dz_new = N * ECC2 * sinphi

        # Check convergence requirement
        if abs(dz - dz_new) < epsilon
            break
        end

        dz = dz_new
    end

    # Extract geodetic coordinates
    zdz = z + dz
    lat = atan(zdz, sqrt(rho2))
    lon = atan(y, x)
    alt = sqrt(rho2 + zdz^2) - N

    # Convert output to degrees
    if use_degrees
        lat = lat*180.0/pi
        lon = lon*180.0/pi
    end

    return [lon, lat, alt]
end

#######
# ENZ #
#######

export rECEFtoENZ
"""
Compute the rotation matrix from the Earth-fixed to the East-North-Up
coorindate basis.

Arguments:
- `station_ecef::Array{<:Real, 1}`: Earth-fixed cartesian station coordinates
- `conversion::Bool`: Conversion type to use. Can be "geocentric" or "geodetic"

Returns:
- `E::Array{Real, 2}`: Topocentric rotation matrix
"""
function rECEFtoENZ(ecef::Array{<:Real, 1} ; conversion::String="geodetic")   
    if length(ecef) < 3
        throw(ArgumentError("Input coordinates must be length 3."))
    end

    # Compute Station Lat-Lon-Altitude
    if conversion == "geodetic"
        lat, lon, = sECEFtoGEOD(ecef, use_degrees=true)
    elseif conversion == "geocentric"
        lat, lon, = sECEFtoGEOC(ecef, use_degrees=true)
    else
        throw(ArgumentError("Unknown conversion method: $conversion"))
    end

    # Compute ENZ basis vectors
    eE = [-sin(lon) ; cos(lon) ; 0]
    eN = [-sin(lat)*cos(lon) ; -sin(lat)*sin(lon) ; cos(lat)]
    eZ = [cos(lat)*cos(lon) ; cos(lat)*sin(lon) ; sin(lat)]

    # Construct Rotation matrix
    E = hcat(eE, eN, eZ)'

    # Return Result
    return E
end

export rENZtoECEF
"""
Compute the rotation matrix from the Earth-fixed to the South-East-Zenith 
coorindate basis.

Arguments:
- `station_ecef::Array{<:Real, 1}`: Earth-fixed cartesian station coordinates
- `conversion::Bool`: Conversion type to use. Can be "geocentric" or "geodetic"

Returns:
- `E::Array{Float64, 2}`: Topocentric rotation matrix
"""
function rENZtoECEF(ecef::Array{<:Real, 1} ; conversion::String="geodetic")
    # Check input coordinates
    if length(ecef) < 3
        throw(ArgumentError("Input coordinates must be length 3."))
    end

    return rECEFtoENZ(ecef, conversion=conversion)'
end

export sECEFtoENZ
"""
Compute the coordinates of an object in the topocentric frame of an
Earth-fixed frame

Arguments:
- `station_ecef::Array{<:Real, 1}`: Earth-fixed cartesian station coordinates
- `ecef::Array{<:Real, 1}`: Coordinates of the object in Earth-fixed station
- `conversion::Bool`: Conversion type to use. Can be "geocentric" or "geodetic"

Returns:
- `E::Array{Float64, 2}`: Topocentric rotation matrix
"""
function sECEFtoENZ(station_ecef::Array{<:Real, 1}, ecef::Array{<:Real, 1} ; conversion::String="geodetic")
    # Check input sizes
    if length(ecef) < 3
        throw(ArgumentError("Input ecef state must be at least length 3."))
    end

    if length(station_ecef) < 3
        throw(ArgumentError("Input station coordinates must be length 3."))
    end

    # Compute ENZ Rotation matrix
    E = rECEFtoENZ(station_ecef, conversion=conversion)

    # Transform range
    range_ecef = ecef[1:3] - station_ecef
    range_enz  = E * range_ecef

    # Transform range-rate (if necessary)
    if length(ecef) == 6
        range_rate_ecef = ecef[4:6]
        range_rate_enz  = E * range_rate_ecef
    end

    # Return
    if length(ecef) == 6
        sat_enz = vcat(range_enz, range_rate_enz)
    else
        sat_enz = range_enz
    end
    
    return sat_enz
end

export sENZtoECEF
"""
Compute the coordinates of an object in the topocentric frame of an
Earth-fixed frame

Arguments:
- `station_ecef::Array{<:Real, 1}`: Earth-fixed cartesian station coordinates
- `sez::Array{<:Real, 1}`: SEZ coordinates of the object
- `conversion::Bool`: Conversion type to use. Can be "geocentric" or "geodetic"

Returns:
- `E::Array{Float64, 2}`: Topocentric rotation matrix
"""
function sENZtoECEF(station_ecef::Array{<:Real, 1}, enz::Array{<:Real, 1} ; conversion::String="geodetic")
    # Check input sizes
    if length(enz) < 3
        throw(ArgumentError("Input ENZ state must be at least length 3."))
    end

    if length(station_ecef) < 3
        throw(ArgumentError("Input station coordinates must be length 3."))
    end

    # Compute ENZ Rotation matrix
    E = rENZtoECEF(station_ecef, conversion=conversion)

    # Transform range
    range_enz  = enz[1:3]
    range_ecef = E * range_enz

    # Transform range-rate (if necessary)
    if length(enz) == 6
        range_rate_enz  = enz[4:6]
        range_rate_ecef = E * range_rate_enz
    end

    # Return
    if length(enz) == 6
        sat_ecef = vcat(range_ecef + station_ecef, range_rate_ecef)
    else
        sat_ecef = range_ecef + station_ecef
    end
    
    return sat_ecef
end

#######
# SEZ #
#######

export rECEFtoSEZ
"""
Compute the rotation matrix from the Earth-fixed to the South-East-Zenith 
coorindate basis.

Arguments:
- `station_ecef::Array{<:Real, 1}`: Earth-fixed cartesian station coordinates
- `conversion::Bool`: Conversion type to use. Can be "geocentric" or "geodetic"

Returns:
- `E::Array{Float64, 2}`: Topocentric rotation matrix
"""
function rECEFtoSEZ(ecef::Array{<:Real, 1} ; conversion::String="geodetic")    
    if length(ecef) < 3
        throw(ArgumentError("Input coordinates must be length 3."))
    end

    # Compute Station Lat-Lon-Altitude
    if conversion == "geodetic"
        lat, lon, = sECEFtoGEOD(ecef, use_degrees=true)
    elseif conversion == "geocentric"
        lat, lon, = sECEFtoGEOC(ecef, use_degrees=true)
    else
        throw(ArgumentError("Unknown conversion method: $conversion"))
    end

    # Compute SEZ basis vectors
    eS = [sin(lat)*cos(lon) ; sin(lat)*sin(lon) ; -cos(lat)]
    eE = [-sin(lon) ; cos(lon) ; 0]
    eZ = [cos(lat)*cos(lon) ; cos(lat)*sin(lon) ; sin(lat)]

    # Construct Rotation matrix
    E = hcat(eS, eE, eZ)'

    # Return Result
    return E
end


export rSEZtoECEF
"""
Compute the rotation matrix from the Earth-fixed to the South-East-Zenith 
coorindate basis.

Arguments:
- `station_ecef::Array{<:Real, 1}`: Earth-fixed cartesian station coordinates
- `conversion::Bool`: Conversion type to use. Can be "geocentric" or "geodetic"

Returns:
- `E::Array{Float64, 2}`: Topocentric rotation matrix
"""
function rSEZtoECEF(ecef::Array{<:Real, 1} ; conversion::String="geodetic")
    # Check input coordinates
    if length(ecef) < 3
        throw(ArgumentError("Input coordinates must be length 3."))
    end

    return rECEFtoSEZ(ecef, conversion=conversion)'
end


export sECEFtoSEZ
"""
Compute the coordinates of an object in the topocentric frame of an
Earth-fixed frame

Arguments:
- `station_ecef::Array{<:Real, 1}`: Earth-fixed cartesian station coordinates
- `ecef::Array{<:Real, 1}`: Coordinates of the object in Earth-fixed station
- `conversion::Bool`: Conversion type to use. Can be "geocentric" or "geodetic"

Returns:
- `E::Array{Float64, 2}`: Topocentric rotation matrix
"""
function sECEFtoSEZ(station_ecef::Array{<:Real, 1}, ecef::Array{<:Real, 1} ; conversion::String="geodetic")    
    # Check input sizes
    if length(ecef) < 3
        throw(ArgumentError("Input ecef state must be at least length 3."))
    end

    if length(station_ecef) < 3
        throw(ArgumentError("Input station coordinates must be length 3."))
    end

    # Construct SEZ Rotation matrix
    E = rECEFtoSEZ(station_ecef, conversion=conversion)

    # Transform range
    range_ecef = ecef[1:3] - station_ecef
    range_sez  = E * range_ecef

    # Transform range-rate (if necessary)
    if length(ecef) == 6
        range_rate_ecef = ecef[4:6]
        range_rate_sez  = E * range_rate_ecef
    end

    # Return
    if length(ecef) == 6
        sez = vcat(range_sez, range_rate_sez)
    else
        sez = range_sez
    end
    
    return sez
end

export sSEZtoECEF
"""
Compute the coordinates of an object in the topocentric frame of an
Earth-fixed frame

Arguments:
- `station_ecef::Array{<:Real, 1}`: Earth-fixed cartesian station coordinates
- `sez::Array{<:Real, 1}`: SEZ coordinates of the object
- `conversion::Bool`: Conversion type to use. Can be "geocentric" or "geodetic"

Returns:
- `E::Array{Float64, 2}`: Topocentric rotation matrix
"""
function sSEZtoECEF(station_ecef::Array{<:Real, 1}, sez::Array{<:Real, 1} ; conversion::String="geodetic")
    # Check input sizes
    if length(sez) < 3
        throw(ArgumentError("Input SEZ state must be at least length 3."))
    end

    if length(station_ecef) < 3
        throw(ArgumentError("Input station coordinates must be length 3."))
    end

    # Compute ENZ Rotation matrix
    E = rSEZtoECEF(station_ecef, conversion=conversion)

    # Transform range
    range_sez  = sez[1:3]
    range_ecef = E * range_sez

    # Transform range-rate (if necessary)
    if length(sez) >= 6
        range_rate_sez  = sez[4:6]
        range_rate_ecef = E * range_rate_sez
    end

    # Return
    if length(sez) >= 6
        sat_ecef = vcat(range_ecef + station_ecef, range_rate_ecef)
    else
        sat_ecef = range_ecef + station_ecef
    end
    
    return sat_ecef
end

###############
# Topocentric #
###############

export sENZtoAZEL
"""
Convert East-North-Zenith topocentric state to azimuth, elevation, and range.

Arguments:
- `x::Array{<:Real, 1}`: East-North-Up state
- `use_degrees:Bool`: If `true` returns result in units of degrees

Returns:
- `azel::Array{<:Real, 1}`: Azimuth, elevation and range [rad; rad; m]
"""
function sENZtoAZEL(x::Array{<:Real, 1} ; use_degrees::Bool=false)
    # Check inputs
    if !(length(x) == 3 || length(x) == 6)
        throw(ArgumentError("Input ENZ state must be length 3 or 6."))
    end

    # Expand values
    rE, rN, rZ = x[1], x[2], x[3]        
    
    # Range
    rho = norm(x[1:3])

    # Elevation
    el = atan(rZ, sqrt(rE^2 + rN^2))

    # Azimuth
    az = 0.0
    if el != pi/2 # Non-singular azimuth 
        az = atan(rE, rN)
        if az < 0
            az += 2*pi
        end
    else # Azimuth may be singular for 90 deg elevation
        if length(x) != 6
            az = 0.0
            # @warn "Could not resolve singularity calculating azimuth."
        else
            # Use rate information to get azimuth if there is a singularity
            # in the position
            az = atan(x[4], x[5])
        end
    end

    # Output
    azel = [az ; el ; rho]

    if use_degrees
        azel[1] *= 180.0/pi
        azel[2] *= 180.0/pi
    end

    # Process Rate information
    if length(x) == 6
        rdE, rdN, rdZ = x[4], x[5], x[6]

        # Range-rate
        rhod = dot(x[1:3], x[4:6])/rho

        # Elevation-rate
        eld = (rdZ - norm(x[4:6])*sin(el))/sqrt(rE^2 + rN^2)

        # Azimuth-rate
        azd = (rdE*rN - rdN*rE)/(rE^2 + rN^2)

        # Output
        azel_rate = [azd ; eld ; rhod]
        if use_degrees
            azel_rate[1] *= 180/pi
            azel_rate[2] *= 180/pi
        end
    end

    # Return
    if length(x) == 6
        return vcat(azel, azel_rate)
    else
        return azel
    end
end

export sSEZtoAZEL
"""
Convert South-East-Zenith topocentric state to azimuth, elevation, and range.

Arguments:
- `x::Array{<:Real, 1}`: South-East-Zenith state
- `use_degrees:Bool`: If `true` returns result in units of degrees

Returns:
- `azel::Array{<:Real, 1}`: Azimuth, elevation and range [rad; rad; m]
"""
function sSEZtoAZEL(x::Array{<:Real, 1} ; use_degrees::Bool=false)
    # Check inputs
    if !(length(x) == 3 || length(x) == 6)
        throw(ArgumentError("Input rECEFtoSEZ state must be length 3 or 6."))
    end

    # Expand values
    rS, rE, rZ = x[1], x[2], x[3]        
    
    # Range
    rho = norm(x[1:3])

    # Elevation
    el = atan(rZ, sqrt(rS^2 + rE^2))

    # Azimuth
    az = 0.0
    if el != pi/2
        az = atan(rE, -rS)
        if az < 0
            az += 2*pi
        end
    else
        if length(x) != 6
            az = 0.0
            # @warn "Could not resolve singularity calculating azimuth."
        else
            # Use rate information to get azimuth if there is a singularity
            # in the position
            az = atan(x[6], -x[5])
        end
    end

    # Output
    azel = [az ; el ; rho]

    if use_degrees
        azel[1] *= 180.0/pi
        azel[2] *= 180.0/pi
    end

    # Process Rate information
    if length(x) == 6
        rdS, rdE, rdZ = x[4], x[5], x[6]

        # Range-rate
        rhod = dot(x[1:3], x[4:6])/rho

        # Elevation-rate
        eld = (rdZ - norm(x[4:6])*sin(el))/sqrt(rS^2 + rE^2)

        # Azimuth-rate
        azd = (rdS*rE - rdE*rS)/(rS^2 + rE^2)

        # Output
        azel_rate = [azd ; eld ; rhod]
        if use_degrees
            azel_rate[1] *= 180/pi
            azel_rate[2] *= 180/pi
        end
    end

    # Return
    if length(x) == 6
        return vcat(azel, azel_rate)
    else
        return azel
    end
end

end