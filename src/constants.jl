__precompile__(true)
module Constants

# Mathematical Constants
export AS2RAD, RAD2AS
const AS2RAD = 2.0*pi/360.0/3600.0
const RAD2AS = 360.0*3600.0/pi/2.0

# Physical Constants
export C_LIGHT, AU
const C_LIGHT     = 299792458.0                 # [m/s]Exact definition Vallado
const AU          = 1.49597870700e11            # [m] Astronomical Unit IAU 2010

# Time Consants
export MJD_ZERO, MJD2000, GPS_TAI, TAI_GPS, TT_TAI, TAI_TT, GPS_TT, TT_GPS, GPS_ZERO
const MJD_ZERO = 2400000.5
const MJD2000  = 51544.0
const GPS_TAI  = -19.0
const TAI_GPS  = -GPS_TAI
const TT_TAI   = 32.184
const TAI_TT   = -TT_TAI
const GPS_TT   = GPS_TAI + TAI_TT
const TT_GPS   = -GPS_TT
const GPS_ZERO = 44244.0

# Earth Constants
export R_EARTH, WGS84_a, WGS84_f, GM_EARTH, e_EARTH, J2_EARTH, OMEGA_EARTH
const R_EARTH     = 6.378136300e6               # [m] GGM05s Value
const WGS84_a     = 6378137.0                   # WGS-84 semi-major axis
const WGS84_f     = 1.0/298.257223563           # WGS-84 flattening
const GM_EARTH    = 3.986004415e14              # [m^3/s^2] GGM05s Value
const e_EARTH     = 8.1819190842622e-2          # [] First Eccentricity WGS84 Value
const J2_EARTH    = 0.0010826358191967          # [] GGM05s value
const OMEGA_EARTH = 7.292115146706979e-5        # [rad/s] Taken from Vallado 4th Ed page 222

# Sun Constants
export GM_SUN, R_SUN, P_SUN, GM_MOON, GM_MERCURY
const GM_SUN      = 132712440041.939400*1e9     # Gravitational constant of the Sun
const R_SUN       = 6.957*1e8                   # Nominal solar radius corresponding to photospheric radius
const P_SUN       = 4.560E-6                    # [N/m^2] (~1367 W/m^2) Solar radiation pressure at 1 AU

# Celestial Constants - from JPL DE430 Ephemerides
export GM_VENUS, GM_MARS, GM_JUPITER, GM_SATURN, GM_URANUS, GM_NEPTUNE, GM_PLUTO
const GM_MOON     = 4902.800066*1e9
const GM_MERCURY  = 22031.780000*1e9
const GM_VENUS    = 324858.592000*1e9
const GM_MARS     = 42828.37521*1e9
const GM_JUPITER  = 126712764.8*1e9
const GM_SATURN   = 37940585.2*1e9
const GM_URANUS   = 5794548.6*1e9
const GM_NEPTUNE  = 6836527.100580*1e9
const GM_PLUTO    = 977.000000*1e9

end