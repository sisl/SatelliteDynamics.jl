__precompile__(true)
module Time

using Printf
using SOFA
using SatelliteDynamics.Constants

################
# Time Methods #
################

export caldate_to_mjd
"""
Convert a Gregorian calendar date to the equivalent Modified Julian Date representation of that time instant.

# Aguments:
- `year::Int` Year
- `year::Int` Month
- `year::Int` Day
- `hour::Int` Hour
- `minute::Int` Minute 
- `second::Real` Seconds
- `nanoseconds::Real` Microseconds

# Returns:
- `mjd::Float64` Modified Julian Date of Epoch
"""
function caldate_to_mjd(year::Int, month::Int, day::Int, hour=0::Int, minute=0::Int, second=0.0::Real, nanosecond=0.0::Real)
    status, jd, fd = iauDtf2d("TAI", year, month, day, hour, minute, second + nanosecond/1.0e9)

    mjd = (jd - MJD_ZERO) + fd

    return mjd
end

export mjd_to_caldate
"""
Convert a Modified Julian Date to the equivalent Gregorian calendar date representation of the same instant in time.

# Aguments:
- `mjd::Real`: Modified Julian Date of Epoch

# Returns:
- `year::Int32`: Year
- `year::Int32`: Month
- `year::Int32`: Day
- `hour::Int32`: Hour
- `minute::Int32`: Minute 
- `second::Float64`: Seconds
- `nanoseconds::Float64`: nanosecondss
"""
function mjd_to_caldate(mjd::Real)
    status, iy, im, id, ihmsf = iauD2dtf("TAI", 9, Constants.MJD_ZERO, mjd)

    return iy, im, id, ihmsf[1], ihmsf[2], ihmsf[3], ihmsf[4]/1.0e9
end

export caldate_to_jd
"""
Convert a Gregorian calendar date to the equivalent Julian Date representation of that time instant.

# Aguments:
- `year::Int`: Year
- `year::Int`: Month
- `year::Int`: Day
- `hour::Int`: Hour
- `minute::Int`: Minute 
- `second::Real`: Seconds
- `nanoseconds::Real`: nanosecondss

# Returns:
- `mjd::Float64`: Julian Date of Epoch
"""
function caldate_to_jd(year::Int, month::Int, day::Int, hour=0::Int, minute=0::Int, second=0.0::Real, nanoseconds=0.0::Real)
    status, jd, fd = iauDtf2d("TAI", year, month, day, hour, minute, second + nanoseconds/1.0e6)

    mjd = jd + fd

    return mjd
end

export jd_to_caldate
"""
Convert a Julian Date to the equivalent Gregorian calendar date representation of the same instant in time.

# Aguments:
- `mjd::Real`: Julian Date of Epoch

# Returns:
- `year::Int32`: Year
- `year::Int32`: Month
- `year::Int32`: Day
- `hour::Int32`: Hour
- `minute::Int32`: Minute 
- `second::Float64`: Seconds
- `microsecond::Float64`: Microseconds
"""
function jd_to_caldate(mjd::Real)
    status, iy, im, id, ihmsf = iauD2dtf("TAI", 9, MJD_ZERO, mjd - MJD_ZERO)

    return iy, im, id, ihmsf[1], ihmsf[2], ihmsf[3], ihmsf[4]/1.0e9
end


export elapsed_from_epoch
"""
Compute the number of elapsed seconds since a given Epoch from the day number. Can be used to compute the elapsed time since a given Julian or Modified Julian Date.

# Arguments:
- `day_number::Real`: Day number, can contain fractional days. Asummes all days are a uniform 86400.0 seconds in length.
- `day_epoch::Real`: Modified Julian Date of Epoch

# Returns:
- `t::Float`: Number of elapsed seconds between the Provided Modified
    Julian date and the epoch.
"""
function elapsed_from_epoch(day_number::Real, day_epoch=0::Real)
    return (day_number - day_epoch)*86400.0
end

export days_from_elapsed
"""
Computes the day number in a given time scale given the elapsed time since epoch and the epoch itself.

Assumes all days are counted using a uniform 86400.0 seconds over the time span.

# Arguments:
- `t::Real`: Elapsed seconds since the `day_epoch`.
- `day_epoch::Real`: Day number of the epoch. Common values are `SatelliteDynamics.Constants.MJD_ZERO` (to get the Julian Day number) or `SatelliteDynamics.Constants.MJD2000` (to get Modified Julian Days if reckoning time from January 1, 2000 0H)

# Returns:
- `days::Float`: Number of elapsed days in the time scale.
"""
function days_from_elapsed(t::Real, day_epoch=0::Real)
    return (t/86400.0 + day_epoch)
end

# #########################
# # Time System Alignment #
# #########################

# """
# Internal function to align a Julian Date and fractional days so the
# fractional day value is within +/- 1.0.
# """
# function align_jd_fd(jd::Real, fd::Real)
#     if fd < 0
#         days = - floor(fd)
#         # days = 1
#         jda = jd - days
#         fda = fd + days
#     elseif fd >= 1
#         days = floor(fd)
#         jda  = jd + days
#         fda  = fd - days
#     else
#         jda, fda = jd, fd
#     end

#     return jda, fda
# end


# function time_system_offset(mjd::Real, tsys_src::String, tsys_dst::String)
#     # If no transformatio needed return early
#     if tsys_src == tsys_dst
#         return 0.0
#     end
    
#     offset = 0.0

#     # Convert To TAI 
#     if tsys_src == "GPS"
#         offset += TAI_GPS
#     elseif tsys_src == "TT"
#         offset += TAI_TT
#     elseif tsys_src == "UTC"
#         status, iy, im, id, ihmsf = iauD2dtf("UTC", 6, MJD_ZERO, mjd) # Returns TAI-UTC
#         status, dutc = iauDat(iy, im, id, (ihmsf[1]*3600 + ihmsf[2]*60 + ihmsf[3] + ihmsf[4]/1e6)/86400.0)
#         offset += dutc
#     elseif tsys_src == "UT1"
#         # Convert UT1 -> UTC
#         offset -= SatelliteDynamics.Universe.UT1_UTC(mjd)

#         # Convert UTC -> TAI
#         status, iy, im, id, ihmsf = iauD2dtf("UTC", 6, MJD_ZERO, mjd + offset) # Returns TAI-UTC
#         status, dutc = iauDat(iy, im, id, (ihmsf[1]*3600 + ihmsf[2]*60 + ihmsf[3] + ihmsf[4]/1e6)/86400.0)
#         offset += dutc
#     elseif tsys_src == "TAI"
#         # Do nothing in this case
#     end

#     # Covert from TAI to source
#     if tsys_dst == "GPS"
#         offset += GPS_TAI
#     elseif tsys_dst == "TT"
#         offset += TT_TAI
#     elseif tsys_dst == "UTC"
#         # Initial UTC guess
#         u1, u2 = MJD_ZERO, mjd + offset/86400.0

#         # Iterate to get the UTC time
#         for i in 1:3
#             status, d1, d2 = iauUtctai(u1, u2)

#             # Adjust UTC guess
#             u1 += MJD_ZERO - d1
#             u2 += mjd + offset/86400.0 - d2
#         end

#         # Compute Caldate from two-part date
#         status, iy, im, id, ihmsf = iauD2dtf("UTC", 6, u1, u2)

#         status, dutc = iauDat(iy, im, id, (ihmsf[1]*3600 + ihmsf[2]*60 + ihmsf[3] + ihmsf[4]/1e6)/86400.0)
#         offset -= dutc
#     elseif tsys_dst == "UT1"
#         # Initial UTC guess
#         u1, u2 = MJD_ZERO, mjd + offset/86400.0

#         # Iterate to get the UTC time
#         for i in 1:3
#             status, d1, d2 = iauUtctai(u1, u2)

#             # Adjust UTC guess
#             u1 += MJD_ZERO - d1
#             u2 += mjd + offset/86400.0 - d2
#         end

#         # Compute Caldate from two-part date
#         status, iy, im, id, ihmsf = iauD2dtf("UTC", 6, u1, u2)

#         status, dutc = iauDat(iy, im, id, (ihmsf[1]*3600 + ihmsf[2]*60 + ihmsf[3] + ihmsf[4]/1e6)/86400.0)
#         offset -= dutc

#         # Convert UTC to UT1
#         offset += SatelliteDynamics.Universe.UT1_UTC(u1 + u2 + offset/86400.0 - MJD_ZERO)
#     elseif tsys_dst == "TAI"
#         # Do nothing in this case
#     end

#     return offset
# end

# function convert_time_system(jd_src::Real, fd_src::Real, tsys_src::String, tsys_dst::String)

#     # Convert To TAI 
#     if tsys_src == "GPS"
#         jd, fd = jd_src, fd_src + TAI_GPS/86400.0
#     elseif tsys_src == "TT"
#         jd, fd = iauTttai(jd_src, fd_src)
#     elseif tsys_src == "UTC"
#         jd, fd = iauUtctai(jd_src, fd_src)
#     elseif tsys_src == "UT1"
#         # Convert UT1 -> UTC
#         jd, fd = iauUt1utc(jd_src, fd_src, SatelliteDynamics.Universe.UT1_UTC(jd_src +fd_src - MJD_ZERO))

#         # Convert UTC -> TAI
#         jd, fd = iauUtctai(jd, fd)
#     elseif tsys_src == "TAI"
#         jd, fd = jd_src, fd_src
#     end

#     jd, fd = align_jd_fd(jd, fd)

#     # Covert from TAI to source
#     if tsys_dst == "GPS"
#         jd_dst, fd_dst = jd, fd + GPS_TAI/86400.0
#     elseif tsys_dst == "TT"
#         jd_dst, fd_dst = iauTaitt(jd, fd)
#     elseif tsys_dst == "UTC"
#         jd_dst, fd_dst = iauTaiutc(jd, fd)
#     elseif tsys_dst == "UT1"
#         # Convert TAI to UTC
#         jd_dst, fd_dst = iauTaiutc(jd, fd)
        
#         # Convert UTC to UT1
#         jd_dst, fd_dst = iauUtcut1(jd, fd, SatelliteDynamics.Universe.UT1_UTC(jd + fd - MJD_ZERO))
#     elseif tsys_dst == "TAI"
#         jd_dst, fd_dst = jd, fd
#     end

#     return align_jd_fd(jd_dst, fd_dst)
# end


# ###############
# # Epoch Class #
# ###############

# struct Epoch
#     jd::Int32
#     fd::Float64
#     kahan_c::Float64
# end

# function Epoch(year::Real, month::Real, day::Real, hour::Real, minute::Real, second::Real ; tsys="TAI"::String)
#     status, jd, fd = iauDtf2d("TAI", year, month, day, hour, minute, second)

#     jd, fd = convert_time_system(jd, fd, tsys, "TAI")

#     new(jd, fd, 0.0)
# end

# function Base.show(io::IO, epc::Epoch)
#     year, month, day, hour, minute, second = cal(epc, tsys="UTC")
#     s = @sprintf "Epoch(%02d-%02d-%02dT%02d:%02d:%02.3fZ)" year month day hour minute second;
#     print(io, s)
# end

# function jd(epc::Epoch ; tsys="TAI"::String)
#     jd, fd = convert_time_system(epc.jd, epc.fd, "TAI", tsys)

#     return jd + fd
# end

# function mjd(epc::Epoch ; tsys="TAI"::String)
#     jd, fd = convert_time_system(epc.jd, epc.fd, "TAI", tsys)

#     return (jd + fd - MJD_ZERO)
# end

# function cal(epc::Epoch ; tsys="TAI"::String)
#     jd, fd = convert_time_system(epc.jd, epc.fd, "TAI", tsys)

#     year, month, day, hour, minute, second, nanoseconds = jd_to_caldate(jd+fd)

#     return year, month, day, hour, minute, second + nanoseconds
# end

# function day_of_year(epc::Epoch ; tsys="TAI"::String)
#     jd, fd = convert_time_system(epc.jd, epc.fd, "TAI", tsys)

#     year, month, day, hour, minute, second, = jd_to_caldate(jd+fd)

#     mjd0 = caldate_to_mjd(year, 1, 1)
#     mjd  = caldate_to_mjd(year, month, day, hour, minute, second)

#     # Get day of year
#     doy = mjd - mjd0 + 1.0

#     return doy
# end

# function gmst(epc::Epoch ; use_degrees=false::Bool)
#     uta, utb = convert_time_system(epc.jd, epc.fd, "TAI", "UT1")
#     tta, ttb = convert_time_system(epc.jd, epc.fd, "TAI", "TT")

#     st = iauGmst06(uta, utb, tta, ttb)

#     if use_degrees
#         st *= 180.0/pi
#     end

#     return st
# end

# function gast(epc::Epoch ; use_degrees=false::Bool)
#     uta, utb = convert_time_system(epc.jd, epc.fd, "TAI", "UT1")
#     tta, ttb = convert_time_system(epc.jd, epc.fd, "TAI", "TT")

#     st = iauGst06a(uta, utb, tta, ttb)

#     if use_degrees
#         st *= 180.0/pi
#     end

#     st = 0.0
#     return st
# end

# # Comparison operators
# function Base.:(==)(epc_left::Epoch, epc_right::Epoch)
#     return (epc_left.jd == epc_right.jd) && (abs(epc_left.fd - epc_right.fd) < 10^-4)
# end

# function Base.:(!=)(epc_left::Epoch, epc_right::Epoch)
#     return (epc_left.jd != epc_right.jd) || (abs(epc_left.fd - epc_right.fd) < 10^-4)
# end

# function Base.:(<)(epc_left::Epoch, epc_right::Epoch)
#     return (epc_left.jd < epc_right.jd) || ((epc_left.jd == epc_right.jd) && (epc_left.fd < epc_right.fd)) 
# end

# Base.isless(epc_left::Epoch, epc_right::Epoch) = Base.:(<)(epc_left::Epoch, epc_right::Epoch)


# function Base.:(<=)(epc_left::Epoch, epc_right::Epoch)
#     return (epc_left.jd < epc_right.jd) || (epc_left.jd == epc_right.jd) && ((epc_left.fd < epc_right.fd) || (abs(epc_left.fd - epc_right.fd) < 10^-4))
# end

# function Base.:(>)(epc_left::Epoch, epc_right::Epoch)
#     return (epc_left.jd > epc_right.jd) || ((epc_left.jd == epc_right.jd) && (epc_left.fd > epc_right.fd))
# end

# Base.isgreater(epc_left::Epoch, epc_right::Epoch) = Base.:(>)(epc_left::Epoch, epc_right::Epoch)

# function Base.:(>=)(epc_left::Epoch, epc_right::Epoch)
#     return (epc_left.jd > epc_right.jd) || (epc_left.jd == epc_right.jd) && ((epc_left.fd > epc_right.fd) || (abs(epc_left.fd - epc_right.fd) < 10^-4))
# end

# # Arithmetic operators
# function Base.:+(epc::Epoch, t::Real)
#     y        = t/86400.0 - epc.kahan_c
#     _t       = epc.fd + y
#     kahan_c  = (_t - epc.fd) - y
#     fd       = _t
    
#     jd, fd = epc.jd, fd
#     while abs(fd) >= 1.0
#         jd, fd = align_jd_fd(jd, fd)
#     end

#     return Epoch(jd, fd, kahan_c)
# end

# function Base.:+(t::Real, epc::Epoch)
#     return epc + t
# end

# function Base.:-(epc::Epoch, t::Real)
#     return epc + (-t)
# end

# function Base.:-(epc_left::Epoch, epc_right::Epoch)
#     # Return difference in time between two epochs
#     return (epc_left.jd - epc_right.jd)*86400.0 +(epc_left.fd - epc_right.fd)*86400.0
# end

end