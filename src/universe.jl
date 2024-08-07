###############
# Remote Data #
###############

DATA_DIR = abspath(joinpath(abspath(string(@__DIR__)), "../data"))

function download_file(url, file)
    filepath = abspath(joinpath(DATA_DIR, file))
    tempfilepath = filepath * ".tmp"

    @debug("Downloading datafile. URL: $url DESTINATION: $filepath")
    
    # Remove any temp files
    rm(tempfilepath, force=true)

    # Attempt to download data
    download(url, tempfilepath)

    # Move temporary file into permanent location
    rm(filepath, force=true)
    mv(tempfilepath, filepath)

    @debug("Downloaded all datafiles.")
end

export download_kp
"""
Download geomagnetic indices.

Notes:
1. Data source is Celestract: https://celestrak.com/SpaceData/sw19571001.txt
"""
function download_kp()
    download_file("https://celestrak.com/SpaceData/sw19571001.txt", "sw19571001.txt")
end

export download_solar_flux
"""
Download F10.7cm Solar Flux data.

10.7cm solar flux is the standard measure of solar activity in space weather models.

Notes:
1. Data source is NRC Canada solar flux tables: https://www.spaceweather.gc.ca/solar_flux_data/daily_flux_values/fluxtable.txt
"""
function download_solar_flux()
    @suppress download_file("https://www.spaceweather.gc.ca/solar_flux_data/daily_flux_values/fluxtable.txt", "fluxtable.txt")
end

export download_all_data
"""
Downloads package datafiles into folders `\$PACKAGE_ROOT/DIR`

Downloads the following files:
- IERS C04 IAU2000A Earth Orientation Data
- IERS C04 IAU1980 Earth Orientation Data
- IERS Bulletin A/B IAU2000 Earth Orientation Data
"""
function download_all_data()
    # Earth Orientation Data
    for (product_name, product_data) in eop_products
        (url, dest) = product_data
        download_file(url, basename(dest))
    end

    # Geomagnetic inidies
    download_kp()

    # Solar Flux
    download_solar_flux()
end

##########################
# Earth Orientation Data #
##########################

# Define product dictionary
eop_products = Dict(
    "C04_20" => ("https://datacenter.iers.org/data/latestVersion/EOP_20_C04_one_file_1962-now.txt", abspath(string(@__DIR__), "../data/EOP_20_C04_one_file_1962-now.txt")),
    "FINALS_2000" => ("https://datacenter.iers.org/data/latestVersion/finals.all.iau2000.txt", abspath(string(@__DIR__), "../data/finals.all.iau2000.txt"))
)

export EarthOrientationData
"""
Earth Orientation Data object. This object stores Earth Orientation Parameters (EOP) and Celestial Intermediate Pole (CIP) offsets. 
The object is used to compute the Earth Orientation Parameters for a given Modified Julian Date (MJD).

Attributes:
- `data::Dict{Int, Tuple{Float64, Float64, Float64}}` EOP data stored as a dictionary with MJD as the key and a tuple of (UT1-UTC [s], xp [rad], yp [rad]) as the value.
- `last_eop_data::Int` Last MJD value in the EOP data.
- `dXdYData::Dict{Int, Tuple{Float64, Float64}}` CIP offsets stored as a dictionary with MJD as the key and a tuple of (dX [rad], dY [rad]) as the value.
- `last_dxdy_mjd::Int` Last MJD value in the dX/dY data.
"""
struct EarthOrientationData
    data::Dict{Int, Tuple{Float64, Float64, Float64}}
    last_eop_data::Int
    dXdYData::Dict{Int, Tuple{Float64, Float64}}
    last_dxdy_mjd::Int
end

"""
Load Earth Orientation Data from a file. The caller must specify the product type to interpret the file as.
If the product is not recognized, an error will be thrown.

Arguments:
- `product::Symbol` The IERS product type can be `:C04_20` or `:FINALS_2000`

Returns:
- `EarthOrientationData` Earth Orientation Data object.
"""
function EarthOrientationData(product::Symbol) 
    product_string = string(product)


    # Load in Data from filepath
    if product_string == "FINALS_2000"
        return EarthOrientationDataFromFinals(eop_products[product_string][2])
    elseif product_string == "C04_20"
        return EarthOrientationDataFromC04_20(eop_products[product_string][2])
    else
        error("Unknown symbol $(String(product))")
    end

    return EarthOrientationData(eop_data, last_eop_data, dxdy_data, last_dxdy_mjd)
end

"""
Load Earth Orientation Data from a file. The caller must specify the product type to interpret the file as.
If the product is not recognized, an error will be thrown.

Arguments:
- `filepath::String` Path to the Earth Orientation Data file.
- `product::Symbol` The IERS product type can be `:C04_20` or `:FINALS_2000`

Returns:
- `EarthOrientationData` Earth Orientation Data object.
"""
function EarthOrientationData(filepath::String, product::Symbol)
    if product == :C04_20
        return EarthOrientationDataFromC04_20(filepath)
    elseif product == :FINALS_2000
        return EarthOrientationDataFromFinals(filepath)
    else
        error("Unknown symbol $(String(:product))")
    end
end


"""
Parse the Earth Orientation Data from a C04_20 file.

Arguments:
- `filepath::String` Path to the C04_20 Earth Orientation Data file.

Returns:
- `EarthOrientationData` Earth Orientation Data object.
"""
function EarthOrientationDataFromC04_20(filepath::String)
    eop_data = Dict{Int, Tuple{Float64, Float64, Float64}}()
    last_eop_data = 0

    dxdy_data = Dict{Int, Tuple{Float64, Float64}}()
    last_dxdy_mjd = 0

    try 
        open(filepath, "r") do product_file
            for i in 1:6
                # Read first 6 lines to skip to data
                readline(product_file)
            end

            for line in readlines(product_file)
                split_line = split(line)
                mjd_utc = parse(Float64, split_line[5])

                # Parse always-present EOP data
                ut1_utc = parse(Float64, split_line[8])
                xp      = parse(Float64, split_line[6])*AS2RAD
                yp      = parse(Float64, split_line[7])*AS2RAD

                eop_data[mjd_utc] = (ut1_utc, xp, yp)

                if mjd_utc > last_eop_data
                    last_eop_data = mjd_utc
                end

                # Parse dX and dY data - Guaranteed to be present for C04_20

                dX = parse(Float64, split_line[9]) * AS2RAD
                dY = parse(Float64, split_line[10]) * AS2RAD

                dxdy_data[mjd_utc] = (dX, dY)

                if mjd_utc > last_dxdy_mjd
                    last_dxdy_mjd = mjd_utc
                end
                
            end
        end
    catch e
        error("Failed to parse \"$filepath\" as C04_20 EOP File. Error: $e")
    end

    return EarthOrientationData(eop_data, last_eop_data, dxdy_data, last_dxdy_mjd)
end

"""
Parse the Earth Orientation Data from a Finals 2000 file.

Arguments:
- `filepath::String` Path to the Finals 2000 Earth Orientation Data file.

Returns:
- `EarthOrientationData` Earth Orientation Data object.
"""
function EarthOrientationDataFromFinals(filepath::String)
    # Initialize data storage
    eop_data = Dict{Int, Tuple{Float64, Float64, Float64}}()
    last_eop_data = 0

    dxdy_data = Dict{Int, Tuple{Float64, Float64}}()
    last_dxdy_mjd = 0

    # Load in Data from filepath
    try
        for line in readlines(filepath)
            if line[17] == 'P' || line[17] == 'I'
                mjd_utc = parse(Int, line[8:12])            # MJD (UTC)
                ut1_utc = parse(Float64, line[59:68])         # UT1-UTC [s]
                xp      = parse(Float64, line[19:27])*AS2RAD  # xp [rad]
                yp      = parse(Float64, line[38:46])*AS2RAD  # yp [rad]
                

                eop_data[mjd_utc] = (ut1_utc, xp, yp)

                if mjd_utc > last_eop_data
                    last_eop_data = mjd_utc
                end

                # Parse dX and dY data - Not always present in FINALS_2000

                if length(line) >= 125 && (line[96] == 'I' || line[96] == 'P')
                    dX = parse(Float64, line[100:106]) * AS2RAD * 1.0e-3 # Apply 1.0e-3 scale factor for file
                    dY = parse(Float64, line[119:125]) * AS2RAD * 1.0e-3

                    dxdy_data[mjd_utc] = (dX, dY)

                    if mjd_utc > last_dxdy_mjd
                        last_dxdy_mjd = mjd_utc
                    end
                end
            end
        end
    catch e
        error("Failed to parse \"$filepath\" as Finals 2000 EOP File: $e")
    end

    return EarthOrientationData(eop_data, last_eop_data, dxdy_data, last_dxdy_mjd)
end

# Declare global Earth Orientation Data Object used by Reference System Calls
export EOP
"""
Module-wide global EarthOrientationData object. This data object is used as the
default source of Earth Orientation Data by reference system transformations if
no explicit EarthOrientationData file is provided to those transformations.

This value can be overridden in your own code as follows:

```julia
SatelliteDynamics.EOP = EarthOrientationData(:FINALS_2000)
```

This global variable defaults to use the module's internal version of `"FINALS_2000"` 
if it is not otherwise set/provided.
"""
global EOP = EarthOrientationData(:FINALS_2000)

# Access Methods
export UT1_UTC
"""
Compute the offset between the UT1 and UTC time systems in seconds. If the EarthOrientationData argument is ommitted the function will use the default module-global value.

Arguments:
- `eop::EarthOrientationData` EarthOrientationData object to use to compute the offset
- `mjd::Real` Modified Julian Date in UTC of the Epoch for which the UT1-UTC offset is desired.
- `interp::Bool` Whether to linearly interpolate the parameter data to the input MJD.

Returns:
- `ut1_utc::Float` UT1 - UTC offset. [s] 
"""
function UT1_UTC(eop::EarthOrientationData, mjd::Real; interp::Bool=false)

    if mjd > eop.last_eop_data
        error("Requested time (MJD: $mjd) is beyond the last EOP data point (MJD: $(eop.last_eop_data)). Consider updating downloaded EOP data files (e.g. using `download_all_data()`).")
    end

    if interp
        x1 = floor(mjd)
        x2 = floor(mjd) + 1
        y1 = eop.data[convert(Int, floor(mjd))][1]
        y2 = eop.data[convert(Int, floor(mjd)+1)][1]
        x  = (y2 - y1)/(x2 - x1) * (mjd - x1) + y1
        return x
    else
        return eop.data[convert(Int, floor(mjd))][1]
    end
end

UT1_UTC(mjd::Real; interp::Bool=false) = UT1_UTC(EOP, mjd, interp=interp)

export POLE_LOCATOR
"""
Compute the location of the pole. Returns x- and y- components as a tuple with the units of [radians].  If the EarthOrientationData argument is ommitted the function will use the default module-global value.

Arguments:
- `eop::EarthOrientationData` EarthOrientationData object to use to compute the offset
- `mjd::Real` Modified Julian Date in UTC of the Epoch for which the pole locator is desired.
- `interp::Bool` Whether to linearly interpolate the parameter data to the input MJD.

Returns:
- `pole_locator::Tuple{Float, Float}` (x, y) pole location in radians.
"""
function POLE_LOCATOR(eop::EarthOrientationData, mjd::Real; interp::Bool=false)

    if mjd > eop.last_eop_data
        error("Requested time (MJD: $mjd) is beyond the last EOP data point (MJD: $(eop.last_eop_data)). Consider updating downloaded EOP data files (e.g. using `download_all_data()`).")
    end

    if interp
        x1 = floor(mjd)
        x2 = floor(mjd) + 1

        # Get values converted to array for interpolation
        y1 = [v for v in eop.data[convert(Int, floor(mjd))][2:3]]
        y2 = [v for v in eop.data[convert(Int, floor(mjd)+1)][2:3]]
        x  = (y2 - y1)/(x2 - x1) * (mjd - x1) + y1
        return x
    else
        return eop.data[convert(Int, floor(mjd))][2:3]
    end
end

POLE_LOCATOR(mjd::Real; interp::Bool=false) = POLE_LOCATOR(EOP, mjd, interp=interp)

export XP
"""
Compute the x-component of the pole locator in [radians]. If the first EarthOrientationData argument is ommitted the function will use the default module-global value.

Arguments:
- `eop::EarthOrientationData` EarthOrientationData object to use to compute the offset
- `mjd::Real` Modified Julian Date in UTC of the Epoch for which the xp value is desired.
- `interp::Bool` Whether to linearly interpolate the parameter data to the input MJD.

Returns:
- `xp::Float` x-component of pole locator in radians.
"""
function XP(eop::EarthOrientationData, mjd::Real; interp=false)
    if mjd > eop.last_eop_data
        error("Requested time (MJD: $mjd) is beyond the last EOP data point (MJD: $(eop.last_eop_data)). Consider updating downloaded EOP data files (e.g. using `download_all_data()`).")
    end
    
    if interp
        x1 = floor(mjd)
        x2 = floor(mjd) + 1
        y1 = eop.data[convert(Int, floor(mjd))][2]
        y2 = eop.data[convert(Int, floor(mjd)+1)][2]
        x  = (y2 - y1)/(x2 - x1) * (mjd - x1) + y1
        return x
    else
        return eop.data[convert(Int, floor(mjd))][2]
    end
end

XP(mjd::Real; interp::Bool=false) = XP(EOP, mjd, interp=interp)

export YP
"""
Compute the y-component of the pole locator in [radians]. If the first EarthOrientationData argument is ommitted the function will use the default module-global value.

Arguments:
- `eop::EarthOrientationData` EarthOrientationData object to use to compute the offset
- `mjd::Real` Modified Julian Date in UTC of the Epoch for which the yp value is desired.
- `interp::Bool` Whether to linearly interpolate the parameter data to the input MJD.

Returns:
- `yp::Float` y-component of pole locator in radians.
"""
function YP(eop::EarthOrientationData, mjd::Real; interp::Bool=false)

    if mjd > eop.last_eop_data
        error("Requested time (MJD: $mjd) is beyond the last EOP data point (MJD: $(eop.last_eop_data)). Consider updating downloaded EOP data files (e.g. using `download_all_data()`).")
    end

    if interp
        x1 = floor(mjd)
        x2 = floor(mjd) + 1
        y1 = eop.data[convert(Int, floor(mjd))][3]
        y2 = eop.data[convert(Int, floor(mjd)+1)][3]
        x  = (y2 - y1)/(x2 - x1) * (mjd - x1) + y1
        return x
    else
        return eop.data[convert(Int, floor(mjd))][3]
    end
end

YP(mjd::Real; interp::Bool=false) = YP(EOP, mjd, interp=interp)

export POLE_OFFSETS
"""
Compute the Celestial Intermediate Pole (CIP) offsets, `(xp, xp)`.  If the EarthOrientationData argument is ommitted the function will use the default module-global value.

Arguments:
- `eop::EarthOrientationData` EarthOrientationData object to use to compute the offset
- `mjd::Real` Modified Julian Date in UTC of the Epoch for which the pole locator is desired.
- `interp::Bool` Whether to linearly interpolate the parameter data to the input MJD.

Returns:
- `pole_offsets::Tuple{Float, Float}` (x, y) pole location in radians.
"""
function POLE_OFFSETS(eop::EarthOrientationData, mjd::Real; interp::Bool=false)

    if mjd > eop.last_dxdy_mjd
        error("Requested time (MJD: $mjd) is beyond the last dX/dY data point (MJD: $(eop.last_dxdy_mjd)). Consider updating downloaded EOP data files (e.g. using `download_all_data()`).")
    end

    if interp
        x1 = floor(mjd)
        x2 = floor(mjd) + 1

        # Get values converted to array for interpolation
        y1 = [v for v in eop.data[convert(Int, floor(mjd))][1:2]]
        y2 = [v for v in eop.data[convert(Int, floor(mjd)+1)][1:2]]
        x  = (y2 - y1)/(x2 - x1) * (mjd - x1) + y1
        return x
    else
        return eop.dXdYData[convert(Int, floor(mjd))][1:2]
    end
end

POLE_OFFSETS(mjd::Real; interp::Bool=false) = POLE_OFFSETS(EOP, mjd, interp=interp)

export DX
"""
Compute x-component of the Celestial Intermediate Pole offset (`dX`) in [radians]. If the first EarthOrientationData argument is ommitted the function will use the default module-global value.

Arguments:
- `eop::EarthOrientationData` EarthOrientationData object to use to compute the offset
- `mjd::Real` Modified Julian Date in UTC of the Epoch for which the `dX` value is desired.
- `interp::Bool` Whether to linearly interpolate the parameter data to the input MJD.

Returns:
- `dX::Float` x-component of pole offset in radians.
"""
function DX(eop::EarthOrientationData, mjd::Real; interp=false)

    if mjd > eop.last_dxdy_mjd
        error("Requested time (MJD: $mjd) is beyond the last dX/dY data point (MJD: $(eop.last_dxdy_mjd)). Consider updating downloaded EOP data files (e.g. using `download_all_data()`).")
    end

    if interp
        x1 = floor(mjd)
        x2 = floor(mjd) + 1
        y1 = eop.dXdYData[convert(Int, floor(mjd))][1]
        y2 = eop.dXdYData[convert(Int, floor(mjd)+1)][1]
        x  = (y2 - y1)/(x2 - x1) * (mjd - x1) + y1
        return x
    else
        return eop.dXdYData[convert(Int, floor(mjd))][1]
    end
end

DX(mjd::Real; interp::Bool=false) = DX(EOP, mjd, interp=interp)

export DY
"""
Compute y-component of the Celestial Intermediate Pole offset (`dY`) in [radians]. If the first EarthOrientationData argument is ommitted the function will use the default module-global value.

Arguments:
- `eop::EarthOrientationData` EarthOrientationData object to use to compute the offset
- `mjd::Real` Modified Julian Date in UTC of the Epoch for which the `dY` value is desired.
- `interp::Bool` Whether to linearly interpolate the parameter data to the input MJD.

Returns:
- `dY::Float` y-component of pole offset in radians.
"""
function DY(eop::EarthOrientationData, mjd::Real; interp::Bool=false)

    if mjd > eop.last_dxdy_mjd
        error("Requested time (MJD: $mjd) is beyond the last dX/dY data point (MJD: $(eop.last_dxdy_mjd)). Consider updating downloaded EOP data files (e.g. using `download_all_data()`).")
    end

    if interp
        x1 = floor(mjd)
        x2 = floor(mjd) + 1
        y1 = eop.dXdYData[convert(Int, floor(mjd))][2]
        y2 = eop.dXdYData[convert(Int, floor(mjd)+1)][2]
        x  = (y2 - y1)/(x2 - x1) * (mjd - x1) + y1
        return x
    else
        return eop.dXdYData[convert(Int, floor(mjd))][2]
    end
end

DY(mjd::Real; interp::Bool=false) = DY(EOP, mjd, interp=interp)

export set_eop
"""
Set Earth orientation data values for a specific date in the module global EarthOrientationData object.

Arguments:
- `mjd::Real` Modified Julian Date in UTC of the Epoch for which the Earth orientation data is aligned to.
- `ut1_utc::Real` Offset between UT1 and UTC in seconds.
- `xp::Real` x-component of the pole locator in radians.
- `yp::Real` y-component of the pole locator in radians.
"""
function set_eop(mjd::Real, ut1_utc::Real, xp::Real, yp::Real)
    EOP.data[convert(Int, floor(mjd))] = (ut1_utc, xp, yp)
end

export load_eop
"""
Load new Earth orientation data into the module global EarthOrientationData object. The product can be one of the symbols: `:C04_20` or `:FINALS_2000`.

Arguments:
- `product::Symbol` Loads a different set of EarthOrientationData values into the module-wide global EarthOrientationData parameters.
"""
function load_eop(product::Symbol)
    global EOP = EarthOrientationData(product::Symbol) 
end

function load_eop(filepath::String, product::Symbol)
    global EOP = EarthOrientationData(filepath, product) 
end

#################
# Gravity Model #
#################

grav_products = Dict(
    :EGM2008_20 => abspath(string(@__DIR__), "../data/EGM2008_20.gfc"),
    :EGM2008_90 => abspath(string(@__DIR__), "../data/EGM2008_90.gfc"),
    :GGM01S => abspath(string(@__DIR__), "../data/GGM01S.gfc"),
    :GGM01S => abspath(string(@__DIR__), "../data/GGM05S.gfc"),
)


function line_starts_with(line::String, str::String)
    if length(line) > length(str) && line[1:length(str)] == str
        return true
    else
        return false
    end
end

export GravModel
"""
GravModel stores a spherical harmonic gravity field in memory. Can store normalized or denomalized coefficients. Package contains EGM2008, GGM01S, and GGM0S gravity models, as well as the default gravity model of EGM2008 truncated to degree and order 90.

Additional gravity field models can be downloaded from: <http://icgem.gfz-potsdam.de/home>

Arguments:
- `filepath::string` Path to spherical harmonic gravity model file.
"""
struct GravModel
    name::String
    normalized::Bool
    R::Float64
    GM::Float64
    n_max::Int64
    m_max::Int64
    data::AbstractArray{Float64, 2}
end

function GravModel(filepath::String)
    model_name = ""
    normalized = false
    R          = 0.0
    GM         = 0.0
    n_max      = 0.0
    m_max      = 0.0
    data       = Array{Float64, 2}(undef, 1, 1)

    # Parse File
    for line in readlines(filepath)
        # Replace non-standard float formatting in GFC files
        line = replace(line, "D+" => "e+")
        line = replace(line, "D-" => "e-")

        if line_starts_with(line, "modelname")
            model_name = split(line)[2]
        elseif line_starts_with(line, "max_degree")
            n_max = parse(Int16, split(line)[2])
            m_max = n_max
            data = zeros(Float64, n_max+1, m_max+1)
        elseif line_starts_with(line, "earth_gravity_constant")
            GM = parse(Float64, split(line)[2])
        elseif line_starts_with(line, "radius")
            R = parse(Float64, split(line)[2])
        elseif line_starts_with(line, "norm")
            if split(line)[2] == "fully_normalized"
                normalized = true
            else
                normalized = false
            end
        elseif line_starts_with(line, "gfc")
            line_split = split(line)
            n = parse(Int16, line_split[2])
            m = parse(Int16, line_split[3])
            C = parse(Float64, line_split[4])
            S = parse(Float64, line_split[5])

            data[n+1, m+1] = C
            if m != 0
                data[m+1-1, n+1] = S
            end
        end
    end

    # Pre-Allocate 
    return GravModel(model_name, normalized, R, GM, n_max, m_max, data)
end

function GravModel(product_name::Symbol)
    return GravModel(grav_products[product_name])
end


# Declare glrobal Gravity Model used by dynamics model calls
export GRAVITY_MODEL
"""
Module-wide global GravityModel object. This data object is used as the
default spherical harmonic gravity field unless one is otherwise provided.

This value can be overridden in your own code as follows:

```julia
SatelliteDynamics.GravityModel = GravityModel(PATH_TO_YOUR_GRAVITY_MODEL)
```

This global variable defaults to use the module's internal version of the EGM2008 model truncated to order and degree 90, if it is not otherwise set.
"""
global GRAVITY_MODEL = GravModel(abspath(@__DIR__, "../data/EGM2008_90.gfc"))

export load_gravity_model
"""
Load new gravity model into module global EarthOrientationData object.

Arguments:
- `gfc_file::String` File path of gravity field model
"""
function load_gravity_model_file(gfc_file::String)
    global GRAVITY_MODEL = GravModel(gfc_file::String) 
end

"""
Load new gravity model into module global EarthOrientationData object. The product can be one of the symbols: `:EGM2008_20`, `:EGM2008_90`, `:GGM01S`, `:GGM05S`, or the filepath to a text-encoded gravity model file.

Arguments:
- `product_name::Symbol` Symbol of a known gravity field product. Valid ones are: `:EGM2008_20`, `:EGM2008_90`, `:GGM01S`, `:GGM05S`
"""
function load_gravity_model(product_name::Symbol)
    global GRAVITY_MODEL = GravModel(product_name::Symbol) 
end

export GRAV_COEF
function GRAV_COEF(i::Int, j::Int)
    # Offset into matrix to deal with julia indexing
    return GRAVITY_MODEL.data[i+1, j+1]
end