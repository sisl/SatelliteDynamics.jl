function get_kp_decimal(number::String, digit::String)
    if isempty(strip(string(number)))
        number = 0.0
    else
        number = parse(Int, number)
    end
    decimal = 0.0
    if parse(Int, digit) == 3
        decimal = 1.0/3.0
    elseif parse(Int, digit) == 7
        decimal = 2.0/3.0
    end

    return number + decimal
end

export SpaceWeatherData
"""
SpaceWeatherData stores space weather data. The data contains both geomagnetic
index data and solar flux data.

# Attributes
- `geomagnetic::SpaceWeatherData` Geomagnetic index data
"""
struct SpaceWeatherData
    cycle::Dict{Int, Tuple{Int, Int}}
    geomagnetic_data::Dict{Int, Tuple{Array{Float64, 1}, Float64, Array{Float64, 1}, Float64}}
    solarflux_data::Dict{Int, Tuple{Float64, Float64, Float64, Float64}}
end


function SpaceWeatherData(filepath::String)
    cycle_data = Dict{Int, Tuple{Int, Int}}()
    geomag_data = Dict{Int, Tuple{Array{Float64, 1}, Float64, Array{Float64, 1}, Float64}}()
    solarflux_data = Dict{Int, Tuple{Float64, Float64, Float64, Float64}}()

    # Load file and read in data
    open(filepath, "r") do fp
        line = readline(fp)
        while !occursin("BEGIN OBSERVED", line)
            # Skip header
            line = readline(fp)
        end

        # Read first line of observed data
        line = readline(fp)

        while !occursin("END OBSERVED", line)
            # Split the line into individual items
            split_line = split(line)

            # Parse Date Information
            year  = parse(Int, split_line[1])
            month = parse(Int, split_line[2])
            day   = parse(Int, split_line[3])
            mjd   = caldate_to_mjd(year, month, day)

            bsrn  = parse(Int, split_line[4])
            dom   = parse(Int, split_line[5])

            cycle_data[mjd] = (bsrn, dom)

            # Parse Geomagnetic Data

            kp_data = Float64[]
            for i in 6:13
                push!(kp_data, get_kp_decimal(split_line[1], split_line[2]))
            end

            kp_sum = get_kp_decimal(split_line[14][1:2], split_line[14][3])

            ap_data = Float64[]

            for i in 15:22
                push!(ap_data, parse(Float64, split_line[i]))
            end

            ap_avg = parse(Float64, split_line[23])

            geomag_data[mjd] = (kp_data, kp_sum, ap_data, ap_avg)

            # Parse Solar Flux Data
            f107_adj = parse(Float64, split_line[27])
            f107_obs = parse(Float64, split_line[31])
            f107_adj_avg = parse(Float64, split_line[30])
            f107_obs_avg = parse(Float64, split_line[33])

            solarflux_data[mjd] = (f107_obs, f107_adj, f107_obs_avg, f107_adj_avg)

            # Read in next line
            line = readline(fp)
        end

        while !occursin("BEGIN DAILY_PREDICTED", line)
            # Skip header
            line = readline(fp)
        end

        # Read first line of predicted data
        line = readline(fp)

        while !occursin("END DAILY_PREDICTED", line)
            # Split the line into individual items
            split_line = split(line)

            # Parse Date Information
            year  = parse(Int, split_line[1])
            month = parse(Int, split_line[2])
            day   = parse(Int, split_line[3])
            mjd   = caldate_to_mjd(year, month, day)

            bsrn  = parse(Int, split_line[4])
            dom   = parse(Int, split_line[5])

            cycle_data[mjd] = (bsrn, dom)

            # Parse Geomagnetic Data

            kp_data = Float64[]
            for i in 6:13
                push!(kp_data, get_kp_decimal(split_line[1], split_line[2]))
            end

            kp_sum = get_kp_decimal(split_line[14][1:2], split_line[14][3])

            ap_data = Float64[]

            for i in 15:22
                push!(ap_data, parse(Float64, split_line[i]))
            end

            ap_avg = parse(Float64, split_line[23])

            geomag_data[mjd] = (kp_data, kp_sum, ap_data, ap_avg)

            # Parse Solar Flux Data
            f107_adj = parse(Float64, split_line[27])
            f107_obs = parse(Float64, split_line[31])
            f107_adj_avg = parse(Float64, split_line[30])
            f107_obs_avg = parse(Float64, split_line[33])

            solarflux_data[mjd] = (f107_obs, f107_adj, f107_obs_avg, f107_adj_avg)

            # Read in next line
            line = readline(fp)
        end
    end

    return SpaceWeatherData(cycle_data, geomag_data, solarflux_data)
end

export SPACE_WEATHER_DATA
"""
Module-wide global SpaceWeatherData. This data object is used as the
default source of geomagnetic data by dynamics models if no explicit 
SpaceWeatherData file is provided to those transformations.

This value can be overridden in your own code as follows:

```julia
SatelliteDynamics.SpaceWeather.SPACE_WEATHER_DATA = SpaceWeatherData(index_file)
```

This global variable defaults used are pull from the Celestrak Space Weather
data repository: https://celestrak.org/SpaceData/
"""
global SPACE_WEATHER_DATA = SpaceWeatherData(abspath(@__DIR__, "../../data/sw19571001.txt"))

export load_space_weather_data
"""
Load new space weather data into the module-wide global SpaceWeatherData structure.

Arguments:
- `filename::String` Path to the space weather data file.
"""
function load_space_weather_data(filename::String)
    global SPACE_WEATHER_DATA = SpaceWeatherData(filename) 
end

export KpIndices
"""
Retrieve Kp Indices for day in question

Arguments:
- `mjd::Real` Modified Julian date of desired data. Time system of input is UT1
- `epc::Epoch` Epoch of desired input Time system of input is UT1
"""
function KpIndices(swdata::SpaceWeatherData, mjd::Real)
    return swdata.geomagnetic_data[floor(Int, mjd)][1]
end

KpIndices(mjd::Real)  = KpIndices(SpaceWeatherData, mjd)
KpIndices(epc::Epoch) = KpIndices(SpaceWeatherData, mjd(epc, tsys="UT1"))

export KpIndex
"""
Kp index for the given epoch.

Arguments:
- `mjd::Real` Modified Julian date of desired data. Time system of input is UT1
- `epc::Epoch` Epoch of desired input Time system of input is UT1
"""
function KpIndex(swdata::SpaceWeatherData, mjd::Real)
    # Base UT1 day for index
    mjd_ut1  = floor(Int, mjd)

    # Get Index into vector based on 3-hour window
    hour_idx = floor(Int, (mjd - mjd_ut1)*8)
    return swdata.geomagnetic_data[mjd_ut1][1][hour_idx+1]
end

KpIndex(mjd::Real)  = KpIndex(SpaceWeatherData, mjd)
KpIndex(epc::Epoch) = KpIndex(SpaceWeatherData, mjd(epc, tsys="UT1"))

export KpDailyIndex
"""
Retrieve daily Kp index for the day in question.

Arguments:
- `mjd::Real` Modified Julian date of desired data. Time system of input is UT1
- `epc::Epoch` Epoch of desired input Time system of input is UT1
"""
function KpDailyIndex(swdata::SpaceWeatherData, mjd::Real)
    return swdata.geomagnetic_data[floor(Int, mjd)][2]
end

KpDailyIndex(mjd::Real)  = KpDailyIndex(SpaceWeatherData, mjd)
KpDailyIndex(epc::Epoch) = KpDailyIndex(SpaceWeatherData, mjd(epc, tsys="UT1"))

export ApIndices
"""
Retrieve Ap Indices for day in question

Arguments:
- `mjd::Real` Modified Julian date of desired data. Time system of input is UT1
- `epc::Epoch` Epoch of desired input Time system of input is UT1
"""
function ApIndices(swdata::SpaceWeatherData, mjd::Real)
    return swdata.geomagnetic_data[floor(Int, mjd)][3]
end

ApIndices(mjd::Real)  = ApIndices(SpaceWeatherData, mjd)
ApIndices(epc::Epoch) = ApIndices(SpaceWeatherData, mjd(epc, tsys="UT1"))

export ApIndex
"""
Ap index for the given epoch.

Arguments:
- `mjd::Real` Modified Julian date of desired data. Time system of input is UT1
- `epc::Epoch` Epoch of desired input Time system of input is UT1
"""
function ApIndex(swdata::SpaceWeatherData, mjd::Real)
    # Base UT1 day for index
    mjd_ut1  = floor(Int, mjd)

    # Get Index into vector based on 3-hour window
    hour_idx = floor(Int, (mjd - mjd_ut1)*8)
    return swdata.geomagnetic_data[mjd_ut1][3][hour_idx+1]
end

ApIndex(mjd::Real)  = ApIndex(SpaceWeatherData, mjd)
ApIndex(epc::Epoch) = ApIndex(SpaceWeatherData, mjd(epc, tsys="UT1"))

export ApDailyIndex
"""
Retrieve daily Ap index for the day in question.

Arguments:
- `mjd::Real` Modified Julian date of desired data. Time system of input is UT1
- `epc::Epoch` Epoch of desired input Time system of input is UT1
"""
function ApDailyIndex(swdata::SpaceWeatherData, mjd::Real)
    return swdata.geomagnetic_data[floor(Int, mjd)][4]
end

ApDailyIndex(mjd::Real)  = ApDailyIndex(SpaceWeatherData, mjd)
ApDailyIndex(epc::Epoch) = ApDailyIndex(SpaceWeatherData, mjd(epc, tsys="UT1"))

export f107Data
"""
F10.7 cm solar flux data on the day in question.

Arguments:
- `mjd::Real` Modified Julian date of desired data. Time system of input is UT1
- `epc::Epoch` Epoch of desired input Time system of input is UT1

Returns:
- `data::Tuple{Float64, Float64, Float64, Float64}` Flux data on the day in question
Elements are observed, adjusted, 81-observed average, 81-adjusted average.
"""
function f107Data(solarflux::SpaceWeatherData, mjd::Real)
    return solarflux.solarflux_data[floor(Int, mjd)]
end

f107Data(mjd::Real)  = f107Data(SPACE_WEATHER_DATA, mjd)
f107Data(epc::Epoch) = f107Data(SPACE_WEATHER_DATA, mjd(epc, tsys="UT1"))

export f107Observed
"""
Retrieve the F10.7 cm solar flux data observed at the day in question.

Arguments:
- `mjd::Real` Modified Julian date of desired data. Time system of input is UT1
- `epc::Epoch` Epoch of desired input Time system of input is UT1
"""
function f107Observed(solarflux::SpaceWeatherData, mjd::Real)
    return solarflux.solarflux_data[floor(Int, mjd)][1]
end

f107Observed(mjd::Real)  = f107Observed(SPACE_WEATHER_DATA, mjd)
f107Observed(epc::Epoch) = f107Observed(SPACE_WEATHER_DATA, mjd(epc, tsys="UT1"))

export f107Adjusted
"""
Retrieve the F10.7 cm solar flux data adjusted to 1 AU distance on the day in 
question.

Arguments:
- `mjd::Real` Modified Julian date of desired data. Time system of input is UT1
- `epc::Epoch` Epoch of desired input Time system of input is UT1
"""
function f107Adjusted(solarflux::SpaceWeatherData, mjd::Real)
    return solarflux.solarflux_data[floor(Int, mjd)][2]
end

f107Adjusted(mjd::Real)  = f107Adjusted(SPACE_WEATHER_DATA, mjd)
f107Adjusted(epc::Epoch) = f107Adjusted(SPACE_WEATHER_DATA, mjd(epc, tsys="UT1"))

export f107ObservedAvg
"""
Retrieve the 81-day average of the F10.7 cm solar flux data observed at the 
observatory over the previous 81-days.

Arguments:
- `mjd::Real` Modified Julian date of desired data. Time system of input is UT1
- `epc::Epoch` Epoch of desired input Time system of input is UT1
"""
function f107ObservedAvg(solarflux::SpaceWeatherData, mjd::Real)
    return solarflux.solarflux_data[floor(Int, mjd)][3]
end

f107ObservedAvg(mjd::Real)  = f107ObservedAvg(SPACE_WEATHER_DATA, mjd)
f107ObservedAvg(epc::Epoch) = f107ObservedAvg(SPACE_WEATHER_DATA, mjd(epc, tsys="UT1"))

export f107AdjustedAvg
"""
Retrieve the 81-day average of the F10.7 cm solar flux data adjusted to 1 AU 
over the previous 81-days.

Arguments:
- `mjd::Real` Modified Julian date of desired data. Time system of input is UT1
- `epc::Epoch` Epoch of desired input Time system of input is UT1
"""
function f107AdjustedAvg(solarflux::SpaceWeatherData, mjd::Real)
    return solarflux.solarflux_data[floor(Int, mjd)][4]
end

f107AdjustedAvg(mjd::Real)  = f107AdjustedAvg(SPACE_WEATHER_DATA, mjd)
f107AdjustedAvg(epc::Epoch) = f107AdjustedAvg(SPACE_WEATHER_DATA, mjd(epc, tsys="UT1"))