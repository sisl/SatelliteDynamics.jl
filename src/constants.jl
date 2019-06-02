__precompile__(true)
module Constants

# Mathematical Constants
export RAD2DEG
"""
Constant to convert radians to degrees. Equal to 360/2pi. [deg/rad]
"""
const RAD2DEG = 360.0/(2.0*pi)

export DEG2RAD
"""
Constant to convert radians to degrees. Equal to 2pi/360. [rad/deg]
"""
const DEG2RAD = 2.0*pi/360.0

export AS2RAD
"""
Constant to convert arcseconds to radians. Equal to 2pi/(360*3600). [rad/as]
"""
const AS2RAD = 2.0*pi/360.0/3600.0

export RAD2AS
"""
Constant to convert radians to arcseconds. Equal to 2pi/(360*3600) [as/ras]
"""
const RAD2AS = 360.0*3600.0/pi/2.0

# Physical Constants
export C_LIGHT
"""
Speed of light in vacuum. [m/s]

D. Vallado, _Fundamentals of Astrodynamics and Applications_ (4th Ed.), 2010
"""
const C_LIGHT     = 299792458.0                 # [m/s]Exact definition Vallado

export AU
"""
Astronomical Unit. Equal to the mean distance of the Earth from the sun.
TDB-compatible value. [m]

P. Gérard and B. Luzum, IERS Technical Note 36, 2010
"""
const AU          = 1.49597870700e11            # [m] Astronomical Unit IAU 2010

# Time Consants

export SECONDS_IN_DAY
"""
Number of seconds in a day.
"""
const SECONDS_IN_DAY = 86400

export MJD_ZERO
"""
Offset of Modified Julian Days representation with respect to Julian Days. For 
a time, t, MJD_ZERO is equal to:

    MJD_ZERO = t_jd - t_mjd

Where t_jd is the epoch represented in Julian Days, and t_mjd is the epoch in
Modified Julian Days.

O. Montenbruck, and E. Gill, _Satellite Orbits: Models, Methods and 
Applications_, 2012.
"""
const MJD_ZERO = 2400000.5

export MJD2000
"""
Modified Julian Date of January 1, 2000 00:00:00. Value is independent of time
scale.

O. Montenbruck, and E. Gill, _Satellite Orbits: Models, Methods and 
Applications_, 2012.
"""
const MJD2000  = 51544.0

export GPS_TAI
"""
Offset of GPS time system with respect to TAI time system.

O. Montenbruck, and E. Gill, _Satellite Orbits: Models, Methods and 
Applications_, 2012.
"""
const GPS_TAI  = -19.0

export TAI_GPS
"""
Offset of TAI time system with respect to GPS time system.

O. Montenbruck, and E. Gill, _Satellite Orbits: Models, Methods and 
Applications_, 2012.
"""
const TAI_GPS  = -GPS_TAI

export TT_TAI
"""
Offset of TT time system with respect to TAI time system.

O. Montenbruck, and E. Gill, _Satellite Orbits: Models, Methods and 
Applications_, 2012.
"""
const TT_TAI   = 32.184

export TAI_TT
"""
Offset of TAI time system with respect to TT time system.

O. Montenbruck, and E. Gill, _Satellite Orbits: Models, Methods and 
Applications_, 2012.
"""
const TAI_TT   = -TT_TAI

export GPS_TT
"""
Offset of GPS time system with respect to TT time system.

O. Montenbruck, and E. Gill, _Satellite Orbits: Models, Methods and 
Applications_, 2012.
"""
const GPS_TT   = GPS_TAI + TAI_TT

export TT_GPS
"""
Offset of TT time system with respect to GPS time system.

O. Montenbruck, and E. Gill, _Satellite Orbits: Models, Methods and 
Applications_, 2012.
"""
const TT_GPS   = -GPS_TT

export GPS_ZERO
"""
Modified Julian Date of the start of the GPS time system in the GPS time system.
This date was January 6, 1980 0H as reckond in the UTC time system.

O. Montenbruck, and E. Gill, _Satellite Orbits: Models, Methods and 
Applications_, 2012.
"""
const GPS_ZERO = 44244.0

# Earth Constants
export R_EARTH
"""
Earth's equatorial radius. [m]

GGM05s Gravity Model
"""
const R_EARTH     = 6.378136300e6               # [m] GGM05s Value

export WGS84_a
"""
Earth's semi-major axis as defined by the WGS84 geodetic system. [m]

NIMA Technical Report TR8350.2
"""
const WGS84_a     = 6378137.0                   # WGS-84 semi-major axis

export WGS84_f
"""
Earth's ellipsoidal flattening.  WGS84 Value.

NIMA Technical Report TR8350.2
"""
const WGS84_f     = 1.0/298.257223563           # WGS-84 flattening

export GM_EARTH
"""
Earth's Gravitational constant [m^3/s^2]

O. Montenbruck, and E. Gill, _Satellite Orbits: Models, Methods and 
Applications_, 2012.
"""
const GM_EARTH    = 3.986004415e14              # [m^3/s^2] GGM05s Value

export e_EARTH
"""
Earth's first eccentricity. WGS84 Value. [dimensionless]

NIMA Technical Report TR8350.2
"""
const e_EARTH     = 8.1819190842622e-2          # [] First Eccentricity WGS84 Value

export J2_EARTH
"""
Earth's first zonal harmonic. [dimensionless]

GGM05s Gravity Model.
"""
const J2_EARTH    = 0.0010826358191967          # [] GGM05s value

export OMEGA_EARTH
"""
Earth axial rotation rate. [rad/s]

D. Vallado, _Fundamentals of Astrodynamics and Applications_ (4th Ed.), p. 222, 2010
"""
const OMEGA_EARTH = 7.292115146706979e-5        # [rad/s] Taken from Vallado 4th Ed page 222

# Sun Constants
export GM_SUN
"""
Gravitational constant of the Sun. [m^3/s^2]

O. Montenbruck, and E. Gill, _Satellite Orbits: Models, Methods and 
Applications_, 2012.
"""
const GM_SUN      = 132712440041.939400*1e9     # Gravitational constant of the Sun

export R_SUN
"""
Nominal solar photospheric radius. [m]

O. Montenbruck, and E. Gill, _Satellite Orbits: Models, Methods and 
Applications_, 2012.
"""
const R_SUN       = 6.957*1e8                   # Nominal solar radius corresponding to photospheric radius

export P_SUN
"""
Nominal solar radiation pressure at 1 AU. [N/m^2]

O. Montenbruck, and E. Gill, _Satellite Orbits: Models, Methods and 
Applications_, 2012.
"""
const P_SUN       = 4.560E-6                    # [N/m^2] (~1367 W/m^2) Solar radiation pressure at 1 AU

# Celestial Constants - from JPL DE430 Ephemerides
export GM_MOON
"""
Gravitational constant of the Moon. [m^3/s^2]

O. Montenbruck, and E. Gill, _Satellite Orbits: Models, Methods and 
Applications_, 2012.
"""
const GM_MOON     = 4902.800066*1e9

export GM_MERCURY
"""
Gravitational constant of the Mercury. [m^3/s^2]

O. Montenbruck, and E. Gill, _Satellite Orbits: Models, Methods and 
Applications_, 2012.
"""
const GM_MERCURY  = 22031.780000*1e9

export GM_VENUS
"""
Gravitational constant of the Venus. [m^3/s^2]

O. Montenbruck, and E. Gill, _Satellite Orbits: Models, Methods and 
Applications_, 2012.
"""
const GM_VENUS    = 324858.592000*1e9

export GM_MARS
"""
Gravitational constant of the Mars. [m^3/s^2]

O. Montenbruck, and E. Gill, _Satellite Orbits: Models, Methods and 
Applications_, 2012.
"""
const GM_MARS     = 42828.37521*1e9

export GM_JUPITER
"""
Gravitational constant of the Jupiter. [m^3/s^2]

O. Montenbruck, and E. Gill, _Satellite Orbits: Models, Methods and 
Applications_, 2012.
"""
const GM_JUPITER  = 126712764.8*1e9

export GM_SATURN
"""
Gravitational constant of the Saturn. [m^3/s^2]

O. Montenbruck, and E. Gill, _Satellite Orbits: Models, Methods and 
Applications_, 2012.
"""
const GM_SATURN   = 37940585.2*1e9

export GM_URANUS
"""
Gravitational constant of the Uranus. [m^3/s^2]

O. Montenbruck, and E. Gill, _Satellite Orbits: Models, Methods and 
Applications_, 2012.
"""
const GM_URANUS   = 5794548.6*1e9

export GM_NEPTUNE
"""
Gravitational constant of the Neptune. [m^3/s^2]

O. Montenbruck, and E. Gill, _Satellite Orbits: Models, Methods and 
Applications_, 2012.
"""
const GM_NEPTUNE  = 6836527.100580*1e9

export GM_PLUTO
"""
Gravitational constant of the Pluto. [m^3/s^2]

O. Montenbruck, and E. Gill, _Satellite Orbits: Models, Methods and 
Applications_, 2012.
"""
const GM_PLUTO    = 977.000000*1e9

end