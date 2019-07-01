__precompile__(true)
module Time

using Printf
using SOFA
using SatelliteDynamics.Constants
using SatelliteDynamics.Universe: UT1_UTC

#############
# Constants #
#############

"""
Defines valid time systems usable as inputs in the time module, and in functions
pertaining to the `Epoch` class.

Valid systems are: `"GPS"`, `"TAI"`, `"TT"`, `"UTC"`, and `"UT1"`
"""
VALID_TIME_SYSTEMS = ["GPS", "TAI", "TT", "UTC", "UT1"]


####################
# Helper Functions #
####################

string_to_nanoseconds(str::AbstractString) = parse(Float64, str)*10^(9-length(str))

function align_epoch_data(days::Integer, seconds::Integer, nanoseconds::Real)
    # Align nanoseconds to [0, 1.0e9)
    while nanoseconds < 0.0
        nanoseconds += 1.0e9
        seconds     -= 1
    end

    while nanoseconds >= 1.0e9
        nanoseconds -= 1.0e9
        seconds     += 1
    end

    # Align seconds to [0, 86400)
    while seconds < 0
        seconds += 86400
        days    -= 1
    end

    while seconds >= 86400
        seconds -= 86400
        days    += 1
    end

    return days, seconds, nanoseconds
end

################
# Time Methods #
################

export caldate_to_mjd
"""
Convert a Gregorian calendar date to the equivalent Modified Julian Date representation of that time instant.

# Aguments:
- `year::Integer` Year
- `year::Integer` Month
- `year::Integer` Day
- `hour::Integer` Hour
- `minute::Integer` Minute 
- `second::Real` Seconds
- `nanoseconds::Real` Nanoseconds

Returns:
- `mjd::Float64` Modified Julian Date of Epoch
"""
function caldate_to_mjd(year::Integer, month::Integer, day::Integer, hour::Integer=0, minute::Integer=0, second::Real=0.0, nanoseconds::Real=0.0)
    status, jd, fd = iauDtf2d("TAI", year, month, day, hour, minute, second + nanoseconds/1.0e9)

    mjd = (jd - Constants.MJD_ZERO) + fd

    return mjd
end

export mjd_to_caldate
"""
Convert a Modified Julian Date to the equivalent Gregorian calendar date representation of the same instant in time.

# Aguments:
- `mjd::Real`: Modified Julian Date of Epoch

Returns:
- `year::Int`: Year
- `year::Int`: Month
- `year::Int`: Day
- `hour::Int`: Hour
- `minute::Int`: Minute 
- `second::Float64`: Seconds
- `nanoseconds::Float64`: Nanoseconds
"""
function mjd_to_caldate(mjd::Real)
    status, iy, im, id, ihmsf = iauD2dtf("TAI", 9, Constants.MJD_ZERO, mjd)

    return iy, im, id, ihmsf[1], ihmsf[2], ihmsf[3], ihmsf[4]/1.0e9
end

export caldate_to_jd
"""
Convert a Gregorian calendar date to the equivalent Julian Date representation of that time instant.

# Aguments:
- `year::Integer`: Year
- `year::Integer`: Month
- `year::Integer`: Day
- `hour::Integer`: Hour
- `minute::Integer`: Minute 
- `second::Real`: Seconds
- `nanoseconds::Real`: Nanoseconds

Returns:
- `mjd::Float64`: Julian Date of Epoch
"""
function caldate_to_jd(year::Integer, month::Integer, day::Integer, hour::Integer=0, minute::Integer=0, second::Real=0.0, nanoseconds::Real=0.0)
    status, jd, fd = iauDtf2d("TAI", year, month, day, hour, minute, second + nanoseconds/1.0e9)

    jd = jd + fd

    return jd
end

export jd_to_caldate
"""
Convert a Julian Date to the equivalent Gregorian calendar date representation of the same instant in time.

# Aguments:
- `jd::Real`: Julian Date of Epoch

Returns:
- `year::Int`: Year
- `year::Int`: Month
- `year::Int`: Day
- `hour::Int`: Hour
- `minute::Int`: Minute 
- `second::Float64`: Seconds
- `microsecond::Float64`: Nanoseconds
"""
function jd_to_caldate(jd::Real)
    status, iy, im, id, ihmsf = iauD2dtf("TAI", 9, jd, 0)

    return iy, im, id, ihmsf[1], ihmsf[2], ihmsf[3], ihmsf[4]/1.0e9
end


export elapsed_from_epoch
"""
Compute the number of elapsed seconds since a given Epoch from the day number. Can be used to compute the elapsed time since a given Julian or Modified Julian Date.

Arguments:
- `day_number::Real`: Day number, can contain fractional days. Asummes all days are a uniform 86400.0 seconds in length.
- `day_epoch::Real`: Modified Julian Date of Epoch

Returns:
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

Arguments:
- `t::Real`: Elapsed seconds since the `day_epoch`.
- `day_epoch::Real`: Day number of the epoch. Common values are `SatelliteDynamics.Constants.MJD_ZERO` (to get the Julian Day number) or `SatelliteDynamics.Constants.MJD2000` (to get Modified Julian Days if reckoning time from January 1, 2000 0H)

Returns:
- `days::Float`: Number of elapsed days in the time scale.
"""
function days_from_elapsed(t::Real, day_epoch=0::Real)
    return (t/86400.0 + day_epoch)
end

###############
# Epoch Class #
###############

"""
Returns whether a symbol is a valid time system.

Valid systems are: `"GPS"`, `"TAI"`, `"TT"`, `"UTC"`, and `"UT1"`

Arguments:
- `tsys::String`: time system symbol to check for validity.

Returns:
- `valid::Bool`: Return `true` if tsys is a valid time system.
"""
function valid_time_system(tsys::String)
    return tsys in VALID_TIME_SYSTEMS
end

export Epoch
"""
The `Epoch` type represents a single instant in time. It is used throughout the
SatelliteDynamics module. It is meant to provide a clear definition of moments
in time and provide a convenient interface display time in various representations
as well as in differrent time systems. The internal data members are also chosen
such that the representation maintains nanosecond-precision in reprersenation
of time and doesn't accumulate floating-point arithmetic errors larger than
nanoseconds even after centuries.

Supports `+`, `+=`, `-`, and `-=` operators. Two Epoch's can be differenced to
return the time difference between two Epochs. If adding a `Real` number it is
interpreted as an offset in seconds to add to the Epoch.

The class also supports all arithmetic operators: `==`, `!=`, `<`, `<=`, `>`, `>=`

Arguments:
- `year::Int` Year
- `year::Int` Month
- `year::Int` Day
- `hour::Int` Hour (optional)
- `minute::Int` Minute (optional)
- `second::Real` Seconds (optional)
- `nanoseconds::Real` Nanoseconds (optional)
- `tsys::String`: Time system of the epoch at initialization

The Epoch class can be also be initialized from a string. Examples of Valid String constructors are: 

```julia
epc = Epoch("2018-12-20")
epc = Epoch("2018-12-20T16:22:19.0Z")
epc = Epoch("2018-12-20T16:22:19.123Z")
epc = Epoch("2018-12-20T16:22:19.123456789Z")
epc = Epoch("2018-12-20T16:22:19Z")
epc = Epoch("20181220T162219Z")
epc = Epoch("2018-12-01 16:22:19 GPS")
epc = Epoch("2018-12-01 16:22:19.0 GPS")
epc = Epoch("2018-12-01 16:22:19.123 GPS")
epc = Epoch("2018-12-01 16:22:19.123456789 GPS")
```
"""
struct Epoch
    # All days, seconds, and nanoseconds are stored internally in the TAI time scale, conversion to from TAI is done on intpu/output interacting.
    days::Int # Total days [0, ∞)
    seconds::Int # Integer seconds [0, 86400)
    nanoseconds::Float64 # Fractional seconds [0, 1)
    tsys::String # Time system of epoch
end

function Epoch(year::Real, month::Real, day::Real, hour::Real=0, minute::Real=0, 
               second::Real=0, nanosecond::Real=0.0; tsys::String="TAI")
    if !valid_time_system(tsys)
        error("Invalid time system $(String(tsys))")
    end

    # Convert date into days and fractional days in desired time system
    status, jd, fd = iauDtf2d("TAI", year, month, day, hour, minute, 0)

    # Log initial output
    # @debug "iauDtf2d returned ($status, $jd, $fd)"

    if status != 0
        error("Non-zero SOFA status returned initializing Epoch.")
    end

    # Get time system offset based on days and fractional days using SOFA
    tsys_offset     = time_system_offset(jd, fd, tsys, "TAI")
    foffset, offset = modf(tsys_offset)

    # Ensure days is an integer number, add entire fractional component to the
    # fractional days variable
    fdays, days = modf(jd)
    fd          = fd + fdays

    # Convert fractional days into total seconds still retaining fractional part
    seconds = fd * 86400.0
    fsecs, secs = modf(seconds)

    # @debug "seconds: $seconds"
    # @debug "fsecs: $fsecs, secs: $secs"

    # Now trip second of the fractional part
    fsecond, second = modf(second)
    # @debug "fsecond: $fsecond, second: $second"

    # Add the integer parts together
    seconds = secs + second + offset
    # @debug "seconds: $seconds"

    # Convert the fractional parts to nanoseconds
    nanoseconds = nanosecond + fsecs*1e9 + fsecond*1e9 + foffset*1e9
    nanoseconds = Float64(nanoseconds)
    # @debug "nanoseconds: $nanoseconds"

    # @debug "Initializing Epoch: days: $days, seconds: $seconds, nanoseconds: $nanoseconds"
    
    # Ensure type consistency of input:
    days        = Int(days)
    seconds     = Int(seconds)
    nanoseconds = Float64(nanoseconds)

    days, seconds, nanoseconds = align_epoch_data(days, seconds, nanoseconds)
    # @debug "Initializing Epoch: days: $days, seconds: $seconds, nanoseconds: $nanoseconds"

    # Parse inputs into Epoch internal storage
    Epoch(days, seconds, nanoseconds, tsys)
end

# All REGEX used for epoch initialization MUST BE AN EXACT MATCH.
# Otherwise the regex parse improperly match on a partial fit and fail.
VALID_EPOCH_REGEX = [
    r"^(\d{4})\-(\d{2})\-(\d{2})$",
    # r"^(\d{4})\-(\d{2})\-(\d{2})[T](\d{2})\:(\d{2})\:(\d{2})([+-])(\d{2})\:(\d{2})$",
    r"^(\d{4})\-(\d{2})\-(\d{2})[T](\d{2})\:(\d{2})\:(\d{2})[Z]$",
    r"^(\d{4})\-(\d{2})\-(\d{2})[T](\d{2})\:(\d{2})\:(\d{2})[.](\d*)[Z]$",
    r"^(\d{4})(\d{2})(\d{2})[T](\d{2})(\d{2})(\d{2})[Z]$",
    r"^(\d{4})\-(\d{2})\-(\d{2})\s(\d{2})\:(\d{2})\:(\d{2})\.*\s*(\d*)\s*([A-Z]*)$"
]


function Epoch(str::String)
    year        = 0
    month       = 0
    day         = 0
    hour        = 0
    minute      = 0
    second      = 0.0
    nanoseconds = 0.0
    tsys        = "UTC"

    m = nothing
    # Iterate through valid regex string 
    for regex in VALID_EPOCH_REGEX
        m = match(regex, str)
        if !(m === nothing)
            # @debug m
            # Parse date (common to all)
            year  = parse(Int, m[1])
            month = parse(Int, m[2])
            day   = parse(Int, m[3])

            # Parse time (most have this)
            if length(m.captures) >= 6
                hour   = parse(Int, m[4])
                minute = parse(Int, m[5])
                second = parse(Float64, m[6])
            end

            # Parse additional types
            if length(m.captures) == 7
                nanoseconds = string_to_nanoseconds(m[7])
            elseif length(m.captures) == 8
                if m[7] != ""
                    nanoseconds = string_to_nanoseconds(m[7])
                end
                tsys = string(m[8])

                if !(tsys in VALID_TIME_SYSTEMS)
                    throw(ErrorException("Parsed invalid time system: \"$tsys\""))
                end
            end

            # Exit early since a match has been found
            break
        end
    end

    # No valid match found throw error
    if m === nothing
        error("Invalid Epoch string. Must be iso8061 compliant.")
    end

    # @debug "regex parsed: $year $month $day $hour $minute $second $nanoseconds $(string(tsys))"

    return Epoch(year, month, day, hour, minute, second, nanoseconds, tsys=tsys)
end

####################
# Epoch Arithmetic #
####################

function Base.:+(epc::Epoch, t::Real)
    # Immidiately separate seconds and fractional seconds
    fseconds, seconds = modf(t)
    seconds = Int(seconds) # Conver seconds to integer seconds

    # Compute time delta aligned 
    dt_days        = div(seconds, 86400)
    dt_seconds     = seconds % 86400
    dt_nanoseconds = fseconds*1e9

    # Perform additon to get new epoch
    days        = epc.days + dt_days
    seconds     = epc.seconds + dt_seconds
    nanoseconds = epc.nanoseconds + dt_nanoseconds

    # Align to proper ranges in reverse order
    days, seconds, nanoseconds = align_epoch_data(days, seconds, nanoseconds)

    return Epoch(days, seconds, nanoseconds, epc.tsys)
end

function Base.:+(t::Real, epc::Epoch)
    return epc + t
end

function Base.:-(epc::Epoch, t::Real)
    return epc + (-t)
end

function Base.:-(epc_left::Epoch, epc_right::Epoch)
    # Return difference in time between two epochs
    return (epc_left.days - epc_right.days)*86400.0 + (epc_left.seconds - epc_right.seconds) + (epc_left.nanoseconds - epc_right.nanoseconds)/1e9
end

##########
# Output #
##########


function Base.show(io::IO, epc::Epoch)
    year, month, day, hour, minute, second, nanodseconds = caldate(epc, tsys="UTC")

    s = @sprintf "Epoch(%02d-%02d-%02dT%02d:%02d:%02d.%03dZ)" year month day hour minute second floor(nanodseconds/1.0e6);

    print(io, s)
end

export epoch_to_jdfd
"""
Compute the two-part date format used by SOFA.jl functions forr a given Epoch.

Arguments:
- `epc::Epoch`: Epoch
- `tsys::String`: Time system to return output in

Returns:
- `d1::Real`: First part of two part date. [days]
- `d2::Real`: Second part of two part date. [days]
"""
function epoch_to_jdfd(epc::Epoch; tsys::String=epc.tsys)
    offset = time_system_offset(epc, "TAI", tsys)

    return epc.days, (epc.seconds + offset + epc.nanoseconds/1.0e9)/86400.0
end

export caldate
"""
Return the Gregorian calendar date for a specific 

Arguments:
- `epc::Epoch`: Input epoch
- `tsys::String`: Time system to compute output in.

Returns:
- `year::Int`: Year of epoch
- `month::Int`: Month of epoch
- `day::Int`: Day of epoch
- `hour::Int`: Hour of epoch
- `minute::Int`: Minute of epoch
- `second::Int`: Second of epoch
- `nanoseconds::Int`: Year of epoch
"""
function caldate(epc::Epoch; tsys::String=epc.tsys)
    offset = time_system_offset(epc, "TAI", tsys)

    status, iy, im, id, ihmsf = iauD2dtf("TAI", 9, 
                                    epc.days, 
                                    (epc.seconds + offset + epc.nanoseconds/1.0e9)/86400.0)

    return iy, im, id, ihmsf[1], ihmsf[2], ihmsf[3], ihmsf[4]
end

export jd
"""
Compute the Julian Date for a specific epoch

Arguments:
- `epc::Epoch`: Epoch
- `tsys::String`: Time system to return output in

Returns:
- `jd::Real`: Julian date of the epoch in the requested time system
"""
function jd(epc::Epoch; tsys::String=epc.tsys)
    offset = time_system_offset(epc, "TAI", tsys)

    return (epc.days + (epc.seconds + epc.nanoseconds/1.0e9 + offset)/86400.0)
end

export mjd
"""
Compute the Modified Julian Date for a specific epoch

Arguments:
- `epc::Epoch`: Epoch
- `tsys::String`: Time system to return output in

Returns:
- `mjd::Real`: Julian date of the epoch in the requested time system
"""
function mjd(epc::Epoch; tsys::String=epc.tsys)
    offset = time_system_offset(epc, "TAI", tsys)

    return (epc.days + (epc.seconds + epc.nanoseconds/1.0e9 + offset)/86400.0) - Constants.MJD_ZERO
end


export day_of_year
"""
Return the day-of-year number for a given `Epoch`. 

January 1 0h of each year will return 1.

Arguments:
- `epc::Epoch`: Epoch
- `tsys::String`: Time system to return output in

Returns:
- `doy::Real`: Day of year number. 
"""
function day_of_year(epc::Epoch; tsys::String=epc.tsys)
    # Compute MJD of first day of yearr
    year, month, day, hour, minute, second, = caldate(epc, tsys=tsys)
    mjd0 = caldate_to_mjd(year, 1, 1)

    # Compute MJD of current day
    offset = time_system_offset(epc, "TAI", tsys)
    mjd = (epc.days + (epc.seconds + epc.nanoseconds/1.0e9 + offset)/86400.0) - Constants.MJD_ZERO

    # Get day of year is the difference
    doy = mjd - mjd0 + 1.0

    return doy
end

export gmst
"""
Compute the Greenwich Mean Sidereal Time for the given Epoch.

Arguments:
- `epc::Epoch`: Epoch
- `use_degrees::Bool`: Return output in degrees (Default: false)

Returns:
- `gmst::Real`: Greenwich Mean Sidereal Time [rad/deg]
"""
function gmst(epc::Epoch; use_degrees::Bool=false)
    uta, utb = epoch_to_jdfd(epc, tsys="UT1")
    tta, ttb = epoch_to_jdfd(epc, tsys="TT")

    gmst = iauGmst06(uta, utb, tta, ttb)

    
    return use_degrees ? gmst*180.0/pi : gmst
end

export gast
"""
Compute the Greenwich Mean Sidereal Time for the given Epoch.

Arguments:
- `epc::Epoch`: Epoch
- `use_degrees::Bool`: Return output in degrees (Default: false)

Returns:
- `gast::Real`: Greenwich Apparent Sidereal Time [rad/deg]
"""
function gast(epc::Epoch; use_degrees::Bool=false)
    uta, utb = epoch_to_jdfd(epc, tsys="UT1")
    tta, ttb = epoch_to_jdfd(epc, tsys="TT")

    gast = iauGst06a(uta, utb, tta, ttb)

    return use_degrees ? gast*180.0/pi : gast
end

# Hash function
function Base.hash(epc::Epoch, h::UInt)
    return hash(epc.days, hash(epc.seconds, hash(epc.nanoseconds, hash(:Epoch, h))))
    # return hash(epc.days, hash(epc.seconds, hash(:Epoch, h)))
end

# Comparison operators
function Base.:(==)(epc_left::Epoch, epc_right::Epoch)
    return ((epc_left.days == epc_right.days) &&
            (epc_left.seconds == epc_right.seconds) &&
            isapprox(epc_left.nanoseconds, epc_right.nanoseconds, atol=1e-3))
end

function Base.:(!=)(epc_left::Epoch, epc_right::Epoch)
    return ((epc_left.days != epc_right.days) ||
            (epc_left.seconds != epc_right.seconds) ||
            !isapprox(epc_left.nanoseconds, epc_right.nanoseconds, atol=1e-3))
end

function Base.:(<)(epc_left::Epoch, epc_right::Epoch)
    return ((epc_left.days < epc_right.days) ||
            ((epc_left.days == epc_right.days) &&
             (epc_left.seconds < epc_right.seconds)) ||
             ((epc_left.days == epc_right.days) &&
             (epc_left.seconds == epc_right.seconds) &&
             (epc_left.nanoseconds < epc_right.nanoseconds)))
end

Base.isless(epc_left::Epoch, epc_right::Epoch) = Base.:(<)(epc_left::Epoch, epc_right::Epoch)

function Base.:(<=)(epc_left::Epoch, epc_right::Epoch)
    return ((epc_left.days < epc_right.days) ||
            ((epc_left.days == epc_right.days) &&
            (epc_left.seconds < epc_right.seconds)) ||
            ((epc_left.days == epc_right.days) &&
            (epc_left.seconds == epc_right.seconds) &&
            (epc_left.nanoseconds <= epc_right.nanoseconds)))
end

function Base.:(>)(epc_left::Epoch, epc_right::Epoch)
    return ((epc_left.days > epc_right.days) ||
            ((epc_left.days == epc_right.days) &&
             (epc_left.seconds > epc_right.seconds)) ||
             ((epc_left.days == epc_right.days) &&
             (epc_left.seconds == epc_right.seconds) &&
             (epc_left.nanoseconds > epc_right.nanoseconds)))

end

Base.isgreater(epc_left::Epoch, epc_right::Epoch) = Base.:(>)(epc_left::Epoch, epc_right::Epoch)

function Base.:(>=)(epc_left::Epoch, epc_right::Epoch)
    return ((epc_left.days > epc_right.days) ||
            ((epc_left.days == epc_right.days) &&
            (epc_left.seconds > epc_right.seconds)) ||
            ((epc_left.days == epc_right.days) &&
            (epc_left.seconds == epc_right.seconds) &&
            (epc_left.nanoseconds >= epc_right.nanoseconds)))
end

#######################
# Time System Offsets #
#######################

export time_system_offset
"""
Compute the offset between two time systems at a given Epoch.

The offset (in seconds) is computed as:

    time_system_offset = tsys_dest - tsys_src

The value returned is the number of seconds that musted be added to the source time system given the input epoch, to get the equivalent epoch.

Conversions are accomplished using SOFA C library calls.
Epoch.

Arguments:
- `jd::Real`: Part 1 of two-part date (Julian days)
- `fd::Real`: Part 2 of two-part date (Fractional days)
- `tsys_src::String`: Base time system
- `tsys_dest::String`: Destination time system

Returns:
- `offset::Float`: Offset between soruce and destination time systems in seconds.
"""
function time_system_offset(jd, fd, tsys_src::String, tsys_dest::String)
    # If no transformation is needed needed return early
    if tsys_src == tsys_dest
        return 0.0
    end
    
    offset = 0.0

    # Convert To TAI 
    if tsys_src == "GPS"
        offset += TAI_GPS
    elseif tsys_src == "TT"
        offset += TAI_TT
    elseif tsys_src == "UTC"
        status, iy, im, id, ihmsf = iauD2dtf("UTC", 6, jd, fd) # Returns TAI-UTC
        status, dutc = iauDat(iy, im, id, (ihmsf[1]*3600 + ihmsf[2]*60 + ihmsf[3] + ihmsf[4]/1e6)/86400.0)
        offset += dutc
    elseif tsys_src == "UT1"
        # Convert UT1 -> UTC
        offset -= UT1_UTC((jd - Constants.MJD_ZERO) + fd)

        # Convert UTC -> TAI
        status, iy, im, id, ihmsf = iauD2dtf("UTC", 6, jd, fd + offset) # Returns TAI-UTC
        status, dutc = iauDat(iy, im, id, (ihmsf[1]*3600 + ihmsf[2]*60 + ihmsf[3] + ihmsf[4]/1e6)/86400.0)
        offset += dutc
    elseif tsys_src == "TAI"
        # Do nothing in this case
    end

    # Covert from TAI to source
    if tsys_dest == "GPS"
        offset += GPS_TAI
    elseif tsys_dest == "TT"
        offset += TT_TAI
    elseif tsys_dest == "UTC"
        # Initial UTC guess
        u1, u2 = jd, fd + offset/86400.0

        # Iterate to get the UTC time
        for i in 1:3
            status, d1, d2 = iauUtctai(u1, u2)

            # Adjust UTC guess
            u1 += jd - d1
            u2 += fd + offset/86400.0 - d2
        end

        # Compute Caldate from two-part date
        status, iy, im, id, ihmsf = iauD2dtf("UTC", 6, u1, u2)

        status, dutc = iauDat(iy, im, id, (ihmsf[1]*3600 + ihmsf[2]*60 + ihmsf[3] + ihmsf[4]/1e6)/86400.0)
        offset -= dutc
    elseif tsys_dest == "UT1"
        # Initial UTC guess
        u1, u2 = jd, fd + offset/86400.0

        # Iterate to get the UTC time
        for i in 1:3
            status, d1, d2 = iauUtctai(u1, u2)

            # Adjust UTC guess
            u1 += jd - d1
            u2 += fd + offset/86400.0 - d2
        end

        # Compute Caldate from two-part date
        status, iy, im, id, ihmsf = iauD2dtf("UTC", 6, u1, u2)

        status, dutc = iauDat(iy, im, id, (ihmsf[1]*3600 + ihmsf[2]*60 + ihmsf[3] + ihmsf[4]/1e6)/86400.0)
        offset -= dutc

        # Convert UTC to UT1
        offset += UT1_UTC(u1 + u2 + offset/86400.0 - Constants.MJD_ZERO)
    elseif tsys_dest == "TAI"
        # Do nothing in this case
    end

    return offset
end

function time_system_offset(epc::Epoch, tsys_src::String, tsys_dest::String)
    jd = epc.days
    fd = (epc.seconds + epc.nanoseconds/1.0e9)/86400.0
    return time_system_offset(jd, fd, tsys_src, tsys_dest)
end

end # Time