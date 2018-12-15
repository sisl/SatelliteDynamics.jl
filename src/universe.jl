__precompile__(true)
module Universe

#################
# Julia Imports #
#################

using HTTP

###################
# Package Imports #
###################

using SatelliteDynamics.Constants

##########################
# Earth Orientation Data #
##########################

valid_products = [:C04_14, :C04_80, :FINALS_2000]

# Define product dictionary
eop_products = Dict(
    :C04_14 => ("https://datacenter.iers.org/data/latestVersion/224_EOP_C04_14.62-NOW.IAU2000A224.txt", abspath(string(@__DIR__), "../data/EOP_C04_14.62-NOW.IAU2000A.txt")),
    :C04_80 => ("https://datacenter.iers.org/data/latestVersion/223_EOP_C04_14.62-NOW.IAU1980223.txt", abspath(string(@__DIR__), "../data/EOP_C04_80.62-NOW.IAU2000A.txt")),
    :FINALS_2000 => ("https://datacenter.iers.org/data/latestVersion/9_FINALS.ALL_IAU2000_V2013_019.txt", abspath(string(@__DIR__), "../data/FINALS.ALL_IAU2000.txt"))
)

export EarthOrientationData
struct EarthOrientationData
    data::Dict
end

function EarthOrientationData(product::Symbol) 
    # Initialize Data Array
    eop_data = Dict{Int32, Tuple{Float64, Float64, Float64}}()

    # Load in Data from filepath
    if product == :FINALS_2000
        for line in readlines(eop_products[product][2])
            if line[17] == 'P' || line[17] == 'I'
                mjd_utc = parse(Int32, line[8:12])            # MJD (UTC)
                ut1_utc = parse(Float64, line[59:68])         # UT1-UTC [s]
                xp      = parse(Float64, line[19:27])*AS2RAD  # xp [rad]
                yp      = parse(Float64, line[38:46])*AS2RAD  # yp [rad]

                eop_data[mjd_utc] = (ut1_utc, xp, yp)
            end
        end
    elseif product == :C04_14 || product == :C04_80
        open(eop_products[product][2], "r") do product_file
            for i in 1:14
                # Read first 14 lines to skip to data
                readline(product_file)
            end

            for line in readlines(product_file)
                split_line = split(line)
                mjd_utc = parse(Int32, split_line[4])
                ut1_utc = parse(Float64, split_line[7])
                xp      = parse(Float64, split_line[5])*AS2RAD
                yp      = parse(Float64, split_line[6])*AS2RAD

                eop_data[mjd_utc] = (ut1_utc, xp, yp)
            end
        end
    else
        @error "Unknown symbol $(String(:product))"
    end

    return EarthOrientationData(eop_data)
end

# Declare global Earth Orientation Data Object used by Reference System Calls
export EOP
global EOP = EarthOrientationData(:FINALS_2000)

# Access Methods
export UT1_UTC
function UT1_UTC(eop::EarthOrientationData, mjd::Real; interp::Bool=false)
    if interp
        x1 = floor(mjd)
        x2 = floor(mjd) + 1
        y1 = eop.data[convert(Int32, floor(mjd))][1]
        y2 = eop.data[convert(Int32, floor(mjd)+1)][1]
        x  = (y2 - y1)/(x2 - x1) * (mjd - x1) + y1
        return x
    else
        return eop.data[convert(Int32, floor(mjd))][1]
    end
end

UT1_UTC(mjd::Real; interp::Bool=false) = UT1_UTC(EOP, mjd, interp=interp)

export POLE_LOCATOR
function POLE_LOCATOR(eop::EarthOrientationData, mjd::Real; interp::Bool=false)
    if interp
        x1 = floor(mjd)
        x2 = floor(mjd) + 1

        # Get values converted to array for interpolation
        y1 = [v for v in eop.data[convert(Int32, floor(mjd))][2:3]]
        y2 = [v for v in eop.data[convert(Int32, floor(mjd)+1)][2:3]]
        x  = (y2 - y1)/(x2 - x1) * (mjd - x1) + y1
        return x
    else
        return eop.data[convert(Int32, floor(mjd))][2:3]
    end
end

POLE_LOCATOR(mjd::Real; interp::Bool=false) = POLE_LOCATOR(EOP, mjd, interp=interp)

export XP
function XP(eop::EarthOrientationData, mjd::Real; interp=false)
    if interp
        x1 = floor(mjd)
        x2 = floor(mjd) + 1
        y1 = eop.data[convert(Int32, floor(mjd))][2]
        y2 = eop.data[convert(Int32, floor(mjd)+1)][2]
        x  = (y2 - y1)/(x2 - x1) * (mjd - x1) + y1
        return x
    else
        return eop.data[convert(Int32, floor(mjd))][2]
    end
end

XP(mjd::Real; interp::Bool=false) = XP(EOP, mjd, interp=interp)

export YP
function YP(eop::EarthOrientationData, mjd::Real; interp::Bool=false)
    if interp
        x1 = floor(mjd)
        x2 = floor(mjd) + 1
        y1 = eop.data[convert(Int32, floor(mjd))][3]
        y2 = eop.data[convert(Int32, floor(mjd)+1)][3]
        x  = (y2 - y1)/(x2 - x1) * (mjd - x1) + y1
        return x
    else
        return eop.data[convert(Int32, floor(mjd))][3]
    end
end

YP(mjd::Real; interp::Bool=false) = YP(EOP, mjd, interp=interp)

export set_eop
function set_eop(mjd::Real, ut1_utc::Float64, xp::Float64, yp::Float64)
    EOP.data[convert(Int32, floor(mjd))] = (ut1_utc, xp*AS2RAD, yp*AS2RAD)
end

export load_eop
function load_eop(product::Symbol)
    global EOP = EarthOrientationData(product::Symbol) 
end

#################
# Gravity Model #
#################

function line_starts_with(line::String, str::String)
    if length(line) > length(str) && line[1:length(str)] == str
        return true
    else
        return false
    end
end

struct GravModel
    name::String
    normalized::Bool
    n_max::Int64
    m_max::Int64
    R::Float64
    GM::Float64
    data::Array{Float64, 2}

    function GravModel(filepath::String)
        model_name = ""
        normalized = false
        R          = 0.0
        GM         = 0.0
        data       = Array{Float64, 2}(undef, 1, 1)
        n_max      = 0.0
        m_max      = 0.0

        # Parse File
        for line in readlines(filepath)
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
        new(model_name, normalized, n_max, m_max, R, GM, data)
    end
end

# Declare glrobal Gravity Model used by dynamics model calls
export GRAVITY_MODEL
global GRAVITY_MODEL = GravModel(abspath(@__DIR__, "../data/EGM2008_90.gfc"))

export load_gravity_model
function load_gravity_model(gfc_file::String)
    global GRAVITY_MODEL = GravModel(gfc_file::String) 
end

export GRAV_COEF
function GRAV_COEF(i::Int, j::Int)
    # Offset into matrix to deal with julia indexing
    return GRAVITY_MODEL.data[i+1, j+1]
end

##########
# Update #
##########

export update_eop
function update_eop(product_name::Symbol)
    if (product_name != :C04_14) && (product_name != :C04_80) && (product_name != :FINALS_2000)
        error("Unknown product type \"$(String(product_name))\"")
    end

    @debug "Getting product: \"$(String(product_name))\""

    product_url  = eop_products[product_name][1]
    product_file = eop_products[product_name][2]

    @debug "IERS Product Server URL: \"$product_url\""
    @debug "Local IERS file location: \"$product_file\""

    HTTP.open("GET", product_url) do http
        open(product_file, "w") do file
            write(file, http)
        end
    end
end

end