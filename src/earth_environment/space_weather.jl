##########################
# Geomagnetic Index Data #
##########################

export GeomagneticIndexData
"""
GeomagneticIndexData stores geomagnetic index data stores the geomagnetic
activity Kp index and the dervied Ap indices.

# Attributes
- `data::Dict`
"""
struct GeomagneticIndexData
    data::Dict{Int, Tuple{Array{Float64, 1}, Float64, Array{Float64, 1}, Float64}}
end

function get_kp_decimal(number, digit)
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

function GeomagneticIndexData(filepath::String)
    data = Dict{Int, Tuple{Array{Float64, 1}, Float64, Array{Float64, 1}, Float64}}()

    # Load file and read in data
    open(filepath, "r") do fp
        for line in readlines(fp)
            if !isempty(line)
                # Day of observation
                year  = 2000 + parse(Int, line[1:2])
                month = parse(Int, line[3:4])
                day   = parse(Int, line[5:6])

                # Data
                kp_data = Float64[]

                for i in 13:2:27
                    push!(kp_data, get_kp_decimal(line[i], line[i+1]))
                end

                kp_daily = get_kp_decimal(line[29:30], line[31])

                # Data
                ap_data = Float64[]

                for i in 32:3:55
                    push!(ap_data, parse(Float64, line[i:(i+2)]))
                end

                ap_daily = parse(Float64, line[56:58])
                
                mjd = round(Int, caldate_to_mjd(year, month, day))
                data[mjd] = (kp_data, kp_daily, ap_data, ap_daily)
            end
        end
    end

    GeomagneticIndexData(data)
end

export GEOMAGNETIC_DATA
"""
Module-wide global GeomagneticIndexData. This data object is used as the
default source of geomagnetic data by dynamics models if no explicit 
GeomagneticIndexData file is provided to those transformations.

This value can be overridden in your own code as follows:

```julia
SatelliteDynamics.SpaceWeather.GEOMAGNETIC_DATA = GeomagneticIndexData(index_file)
```

This global variable defaults used are pull from GFZ Potsdam FTP server:
ftp://ftp.gfz-potsdam.de/pub/home/obs/kp-ap/wdc
"""
global GEOMAGNETIC_DATA = GeomagneticIndexData(abspath(@__DIR__, "../../data/kpall.wdc"))

export KpIndices
"""
Retrieve Kp Indices for day in question

Arguments:
- `mjd::Real` Modified Julian date of desired data. Time system of input is UT1
- `epc::Epoch` Epoch of desired input Time system of input is UT1
"""
function KpIndices(geoindices::GeomagneticIndexData, mjd::Real)
    return geoindices.data[floor(Int, mjd)][1]
end

KpIndices(mjd::Real)  = KpIndices(GEOMAGNETIC_DATA, mjd)
KpIndices(epc::Epoch) = KpIndices(GEOMAGNETIC_DATA, mjd(epc, tsys="UT1"))

export KpIndex
"""
Kp index for the given epoch.

Arguments:
- `mjd::Real` Modified Julian date of desired data. Time system of input is UT1
- `epc::Epoch` Epoch of desired input Time system of input is UT1
"""
function KpIndex(geoindices::GeomagneticIndexData, mjd::Real)
    # Base UT1 day for index
    mjd_ut1  = floor(Int, mjd)

    # Get Index into vector based on 3-hour window
    hour_idx = floor(Int, (mjd - mjd_ut1)*8)
    return geoindices.data[mjd_ut1][1][hour_idx+1]
end

KpIndex(mjd::Real)  = KpIndex(GEOMAGNETIC_DATA, mjd)
KpIndex(epc::Epoch) = KpIndex(GEOMAGNETIC_DATA, mjd(epc, tsys="UT1"))

export KpDailyIndex
"""
Retrieve daily Kp index for the day in question.

Arguments:
- `mjd::Real` Modified Julian date of desired data. Time system of input is UT1
- `epc::Epoch` Epoch of desired input Time system of input is UT1
"""
function KpDailyIndex(geoindices::GeomagneticIndexData, mjd::Real)
    return geoindices.data[floor(Int, mjd)][2]
end

KpDailyIndex(mjd::Real)  = KpDailyIndex(GEOMAGNETIC_DATA, mjd)
KpDailyIndex(epc::Epoch) = KpDailyIndex(GEOMAGNETIC_DATA, mjd(epc, tsys="UT1"))

export ApIndices
"""
Retrieve Ap Indices for day in question

Arguments:
- `mjd::Real` Modified Julian date of desired data. Time system of input is UT1
- `epc::Epoch` Epoch of desired input Time system of input is UT1
"""
function ApIndices(geoindices::GeomagneticIndexData, mjd::Real)
    return geoindices.data[floor(Int, mjd)][3]
end

ApIndices(mjd::Real)  = ApIndices(GEOMAGNETIC_DATA, mjd)
ApIndices(epc::Epoch) = ApIndices(GEOMAGNETIC_DATA, mjd(epc, tsys="UT1"))

export ApIndex
"""
Ap index for the given epoch.

Arguments:
- `mjd::Real` Modified Julian date of desired data. Time system of input is UT1
- `epc::Epoch` Epoch of desired input Time system of input is UT1
"""
function ApIndex(geoindices::GeomagneticIndexData, mjd::Real)
    # Base UT1 day for index
    mjd_ut1  = floor(Int, mjd)

    # Get Index into vector based on 3-hour window
    hour_idx = floor(Int, (mjd - mjd_ut1)*8)
    return geoindices.data[mjd_ut1][3][hour_idx+1]
end

ApIndex(mjd::Real)  = ApIndex(GEOMAGNETIC_DATA, mjd)
ApIndex(epc::Epoch) = ApIndex(GEOMAGNETIC_DATA, mjd(epc, tsys="UT1"))

export ApDailyIndex
"""
Retrieve daily Ap index for the day in question.

Arguments:
- `mjd::Real` Modified Julian date of desired data. Time system of input is UT1
- `epc::Epoch` Epoch of desired input Time system of input is UT1
"""
function ApDailyIndex(geoindices::GeomagneticIndexData, mjd::Real)
    return geoindices.data[floor(Int, mjd)][4]
end

ApDailyIndex(mjd::Real)  = ApDailyIndex(GEOMAGNETIC_DATA, mjd)
ApDailyIndex(epc::Epoch) = ApDailyIndex(GEOMAGNETIC_DATA, mjd(epc, tsys="UT1"))

###################
# Solar Flux Data #
###################

export SolarFluxData
"""
SolarFluxData stores solar flux index data.

The data is a dictionary keyed to the modified julian date in UT of the day of 
observation. All observations are taken at 20:00 UT at the of question. The values
are a tuple storing the solar radio flux the F10.7cm solar radio flux. The first 
entry is the flux observed directly at the Ottowa radio observatory station,
while the second entry is the flux value adjusted to 1 AU. The final two entries
is the 81-day average, cenerted on the day entry, of the observed and adjusted 
solar radio flux.

# Attributes
- `data::Dict{Int, Tuple{Float64, Float64, Float64}` Stores observed F10.7cm
solar radio flux, the adjusted solar radio flux, and the 81-day average centered
on the day in question.
"""
struct SolarFluxData
    data::Dict{Int, Tuple{Float64, Float64, Float64, Float64}}
end

function SolarFluxData(filepath::String)
    data = Dict{Int, Tuple{Float64, Float64, Float64, Float64}}()

    # Used to calculate the average solar radio flux for the days in question
    avg  = []

    # Load file and read in data
    open(filepath, "r") do fp
        for i in 1:2
            # Skip first two lines 
            readline(fp)
        end

        for line in readlines(fp)
            split_line = split(line)

            # Day of observation
            year  = parse(Int, split_line[1][1:4])
            month = parse(Int, split_line[1][5:6])
            day   = parse(Int, split_line[1][7:8])

            # Data
            f107_obs = parse(Float64, split_line[5])
            f107_adj = parse(Float64, split_line[6])

            # Push data to array for computing 81 day average
            push!(avg, [year month day f107_obs f107_adj])
        end
    end

    # Format data in matrix to make it possible to compute 81-day average
    avg = vcat(avg...)
    
    day_avg = 81            # Number of days in average
    day_os  = round(Int, (day_avg-1)/2) # One sided average direction
    n_data  = size(avg)[1]
    for i in 1:n_data
        min_index = max(1, i-day_os)
        max_index = min(n_data, i+day_os)

        n_avg   = 0 # Number of data points in average
        obs_avg = 0 # Observed flux average
        adj_avg = 0 # Adjusted flux average
        for idx in min_index:max_index
            n_avg += 1

            obs_avg += avg[idx, 4]
            adj_avg += avg[idx, 5]
        end

        # Calculate data average
        obs_avg = obs_avg/n_avg
        adj_avg = adj_avg/n_avg

        # Insert entry into dictionary
        year      = round(Int, avg[i, 1])
        month     = round(Int, avg[i, 2])
        day       = round(Int, avg[i, 3])
        mjd       = round(Int, caldate_to_mjd(year, month, day))
        data[mjd] = (avg[i, 4], avg[i, 5], obs_avg, adj_avg)
    end

    SolarFluxData(data)
end

export SOLAR_FLUX_DATA
"""
Module-wide global SolarFluxData object. This data object is used as the
default source of Solar Flux Data by dynamics models if no explicit 
SolarFluxData file is provided to those transformations.

This value can be overridden in your own code as follows:

```julia
SatelliteDynamics.SpaceWeather.SOLAR_FLUX_DATA = SolarFluxData(flux_file)
```

This global variable defaults uses flux data provide by National Resources Canada: 
ftp://ftp.seismo.nrcan.gc.ca/spaceweather/solar_flux/daily_flux_values/fluxtable.txt
"""
global SOLAR_FLUX_DATA = SolarFluxData(abspath(@__DIR__, "../../data/fluxtable.txt"))

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
function f107Data(solarflux::SolarFluxData, mjd::Real)
    return solarflux.data[floor(Int, mjd)]
end

f107Data(mjd::Real)  = f107Data(SOLAR_FLUX_DATA, mjd)
f107Data(epc::Epoch) = f107Data(SOLAR_FLUX_DATA, mjd(epc, tsys="UT1"))

export f107Observed
"""
Retrieve the F10.7 cm solar flux data observed at the day in question.

Arguments:
- `mjd::Real` Modified Julian date of desired data. Time system of input is UT1
- `epc::Epoch` Epoch of desired input Time system of input is UT1
"""
function f107Observed(solarflux::SolarFluxData, mjd::Real)
    return solarflux.data[floor(Int, mjd)][1]
end

f107Observed(mjd::Real)  = f107Observed(SOLAR_FLUX_DATA, mjd)
f107Observed(epc::Epoch) = f107Observed(SOLAR_FLUX_DATA, mjd(epc, tsys="UT1"))

export f107Adjusted
"""
Retrieve the F10.7 cm solar flux data adjusted to 1 AU distance on the day in 
question.

Arguments:
- `mjd::Real` Modified Julian date of desired data. Time system of input is UT1
- `epc::Epoch` Epoch of desired input Time system of input is UT1
"""
function f107Adjusted(solarflux::SolarFluxData, mjd::Real)
    return solarflux.data[floor(Int, mjd)][2]
end

f107Adjusted(mjd::Real)  = f107Adjusted(SOLAR_FLUX_DATA, mjd)
f107Adjusted(epc::Epoch) = f107Adjusted(SOLAR_FLUX_DATA, mjd(epc, tsys="UT1"))

export f107ObservedAvg
"""
Retrieve the 81-day average of the F10.7 cm solar flux data observed at the 
observatory centered on the day in question.

Arguments:
- `mjd::Real` Modified Julian date of desired data. Time system of input is UT1
- `epc::Epoch` Epoch of desired input Time system of input is UT1
"""
function f107ObservedAvg(solarflux::SolarFluxData, mjd::Real)
    return solarflux.data[floor(Int, mjd)][3]
end

f107ObservedAvg(mjd::Real)  = f107ObservedAvg(SOLAR_FLUX_DATA, mjd)
f107ObservedAvg(epc::Epoch) = f107ObservedAvg(SOLAR_FLUX_DATA, mjd(epc, tsys="UT1"))

export f107AdjustedAvg
"""
Retrieve the 81-day average of the F10.7 cm solar flux data adjusted to 1 AU 
centered on the day in question.

Arguments:
- `mjd::Real` Modified Julian date of desired data. Time system of input is UT1
- `epc::Epoch` Epoch of desired input Time system of input is UT1
"""
function f107AdjustedAvg(solarflux::SolarFluxData, mjd::Real)
    return solarflux.data[floor(Int, mjd)][4]
end

f107AdjustedAvg(mjd::Real)  = f107AdjustedAvg(SOLAR_FLUX_DATA, mjd)
f107AdjustedAvg(epc::Epoch) = f107AdjustedAvg(SOLAR_FLUX_DATA, mjd(epc, tsys="UT1"))