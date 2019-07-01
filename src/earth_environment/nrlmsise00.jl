__precompile__(true)
module NRLMSISE00
"""
Julia NRLMSISE00 implementation is Dominik Brodowski's C implmenation adapted
for Julia.

Ref: https://www.brodo.de/space/nrlmsise/
"""

using Printf
using SatelliteDynamics.Time: Epoch, mjd, day_of_year
using SatelliteDynamics.Coordinates: sECEFtoGEOD
using SatelliteDynamics.EarthEnvironment.SpaceWeather

##############
# Model Data #
##############

include("nrlmsise00_data.jl")

###################
# Model Structure #
###################

"""
Internal data structure used for setting model parameters
"""
mutable struct NRLMSISE_Flags
    switches::Array{<:Int, 1}
    sw::Array{<:Real, 1}
    swc::Array{<:Real, 1}

    function NRLMSISE_Flags(switches::Array{<:Int, 1}=zeros(Int, 24),
            sw::Array{<:Real, 1}=zeros(Float32, 24),
            swc::Array{<:Real, 1}=zeros(Float32, 24),)
        new(switches, sw, swc)
    end
end

"""
Internal structure for storing model input

ap_array contains the following magnetic values:
- 1 : daily AP
- 2 : 3 hr AP index for current time 
- 3 : 3 hr AP index for 3 hrs before current time
- 4 : 3 hr AP index for 6 hrs before current time
- 5 : 3 hr AP index for 9 hrs before current time 
- 6 : Average of eight 3 hr AP indices from 12 to 33 hrs prior to current time 
- 7 : Average of eight 3 hr AP indices from 36 to 57 hrs prior to current time 
"""
mutable struct NRLMSISE_Input
    year::Int       # Year, currently ignored
    doy::Int        # Day of year
    sec::Float64    # Seconds in day (UT)
    alt::Float64    # Altitude in [km]
    g_lat::Float64  # Geodetic latitude
    g_lon::Float64  # Geodetic longitude
    lst::Float64    # Local apparent solar time (hours)
    f107A::Float64  # 81 day average of F10.7cm flux (centered on day)
    f107::Float64   # Daily F10.7cm flux for previous day
    ap::Float64     # Magnetic index (daily)
    ap_array::Array{<:Real, 1} # Magnetic index array

    function NRLMSISE_Input(year::Int=2000, doy::Int=1, sec::Float64=0.0, 
                alt::Float64=0.0, g_lat::Float64=0.0, g_lon::Float64=0.0,
                lst::Float64=0.0, f107A::Float64=0.0, f107::Float64=0.0, 
                ap::Float64=0.0, ap_array::Array{<:Real, 1}=zeros(Float64, 7))
        new(year, doy, sec, alt, g_lat, g_lon, lst, f107A, f107, ap, ap_array)
    end
end

"""
# Output variables
d = zeros(Float64, 9)
t = zeros(Float64, 2)
"""
mutable struct NRLMSISE_Output
    d::Array{Float64, 1}
    t::Array{Float64, 1}

    function NRLMSISE_Output(d::Array{<:Real, 1}=zeros(Float64, 9), t::Array{<:Real, 1}=zeros(Float64, 2))
        new(d, t)
    end
end

#####################
# Utility Functions #
#####################

"""
Select internal model flags 
"""
function tselec!(flags::NRLMSISE_Flags)
    for i in 1:24
        if i != 10
            flags.switches[i] == 1 ? flags.sw[i] = 1 : flags.sw[i] = 0
            flags.switches[i] >  0 ? flags.swc[i] = 1 : flags.swc[i] = 0

        else
            flags.sw[i]  = flags.switches[i]
            flags.swc[i] = flags.switches[i]
        end
    end
end


function glatf(lat::Real)
    dgtr = 1.74533E-2

    c2  = cos(2.0*dgtr*lat)

    gv   = 980.616 * (1.0 - 0.0026373 * c2)
    reff = 2.0 * (gv) / (3.085462E-6 + 2.27E-9 * c2) * 1.0E-5

    return gv, reff
end

"""
Chemistry/Dissociation correction for MSIS models

Arguments:
- `alt::Real` Altitude
- `r::Real` Target ratio
- `h1::Real` Transition scale length
- `zh::Real` Altitude of 1/2 R

Returns:
- `e::Real` Correction coefficient
"""
function ccor(alt::Real, r::Real, h1::Real, zh::Real)
    e = (alt - zh)/h1

    if e > 70
        return exp(0)
    elseif e < -70
        return exp(r)
    end

    ex = exp(e)
    e  = r/(1.0+ex)
    return exp(e)
end

"""
Chemistry/Dissociation correction for MSIS models

Arguments:
- `alt::Real` Altitude
- `r::Real` Target ratio
- `h1::Real` Transition scale length
- `zh::Real` Altitude of 1/2 R
- `h2::Real` Transition scale length 2

Returns:
- `e::Real` Correction coefficient
"""
function ccor2(alt::Real, r::Real, h1::Real, zh::Real, h2::Real)
    e1 = (alt - zh)/h1
    e2 = (alt - zh)/h2

    if (e1 > 70) || (e2 > 70)
        return exp(0)
    elseif (e1 < -70) && (e2 < -70) # TODO: Should this be && not ||?
        return exp(r)
    end

    ex1 = exp(e1)
    ex2 = exp(e2)

    ccor2v = r/(1.0+0.5*(ex1+ex2))
    return exp(ccor2v)
end

"""
Compute scale height
"""
function scaleh(alt::Real, xm::Real, temp::Real; gsurf::Real, re::Real)
    rgas = 831.4
    g = gsurf/(1.0 + alt/re)^2
    g = rgas*temp/(g*xm)

    return g
end

"""
Turbopause correction for MSIS models

Arguments:
- `dd::Real` diffusive density
- `dm::Real` full mixed density density
- `zhm::Real` transition scale length
- `xmm::Real` full mixed molecular weight
- `xm::Real` species molecular weight
- `dnet::Real` combined density
"""
function dnet(dd::Real, dm::Real, zhm::Real, xmm::Real, xm::Real)
    a = zhm/(xmm-xm)

    if !((dm>0) && (dd>0))
        @error "dnet log error $dm $dd $xm"
        if (dd==0) && (dm==0)
            dd = 1
        end

        if dm == 0
            return dd
        end

        if dd == 0
            return dm
        end
    end

    ylog = a * log(dm/dd)

    if ylog < -10
        return dd
    end

    if ylog > 10
        return dm
    end

    a = dd*(1.0 + exp(ylog))^(1.0/a)

    return a
end

"""
Integrate cubic spline from xa[1] to x

Arguments:
-`xa::Array{Real, 1}` Array of tabulated function inputs in ascending order by x
-`ya::Array{Real, 1}` Array of tabulated function outputs in ascending order by x
-`y2a::Array{Real, 1}` Array of second derivatives
-`x::Real` Array output point

Returns:
-`y::Real` Array output value
"""
function splini(xa::Array{<:Real, 1}, ya::Array{<:Real, 1}, y2a::Array{<:Real, 1}, n::Int, x::Real)
    yi  = 0
    klo = 1
    khi = 2

    while (x > xa[klo]) && (khi <= n)
        xx = x
        if (khi < n)
            if x < xa[khi]
                xx = x
            else
                xx = xa[khi]
            end
        end
        h    = xa[khi] - xa[klo]
        a    = (xa[khi] - xx)/h
        b    = (xx - xa[klo])/h
        a2   = a^2
        b2   = b^2
        yi  += ((1.0 - a2) * ya[klo] / 2.0 + b2 * ya[khi] / 2.0 + ((-(1.0+a2*a2)/4.0 + a2/2.0) * y2a[klo] + (b2*b2/4.0 - b2/2.0) * y2a[khi]) * h * h / 6.0) * h
        klo += 1
        khi += 1
    end

    return yi
end

"""
Interpolate cubic spline from to x

Arguments:
-`xa::Array{Real, 1}` Array of tabulated function inputs in ascending order by x
-`ya::Array{Real, 1}` Array of tabulated function outputs in ascending order by x
-`y2a::Array{Real, 1}` Array of second derivates
-`x::Real` Array output point

Returns:
-`y::Real` Array output value
"""
function splint(xa::Array{<:Real, 1}, ya::Array{<:Real, 1}, y2a::Array{<:Real, 1}, n::Int, x::Real)
    klo = 1
    khi = n
    k   = 0
    while (khi - klo) > 1
        k = floor(Int, (khi + klo)/2)
        if xa[k] > x
            khi = k
        else
            klo = k
        end
    end

    h = xa[khi] - xa[klo]
    if h == 0
        @error "Bad array input to splint. xa[$khi]: $(xa[khi]), xa[$klo]: $(xa[klo])"
        @error "Array Input: $xa, $ya"
    end

    a = (xa[khi] - x)/h
    b = (x - xa[klo])/h

    y = a * ya[klo] + b * ya[khi] + ((a*a*a - a) * y2a[klo] + (b*b*b - b) * y2a[khi]) * h * h/6.0

    return y
end

"""
Interpolate cubic spline from to x

Arguments:
-`x::Array{Real, 1}` Array of tabulated function inputs in ascending order by x
-`y::Array{Real, 1}` Array of tabulated function outputs in ascending order by x
-`yp1::Real` Specified derivative at x[1]   (yp1 >= 1e30 signals second derivative is zero)
-`ypn::Real` Specified derivative at x[end] (ypn >= 1e30 signals second derivative is zero)
-`x::Real` Array output point

Returns:
-`y::Array{Float64, 1}` Array second derivatives
"""
function spline(x::Array{<:Real, 1}, y::Array{<:Real, 1}, n::Int, yp1::Real, ypn::Real)
    u  = zeros(Float64, n) # Need the + 1 here because C implementation already allocates extra
    y2 = zeros(Float64, n)

    if yp1 > 0.99e30
        y2[1] = 0
        u[1]  = 0
    else
        y2[1] = -0.5
        u[1]  = (3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1)
    end

    # println("u: $u")
    # println("y2: $y2")

    for i in 2:n-1
        sig = (x[i]-x[i-1])/(x[i+1] - x[i-1])
        p     = sig*y2[i-1] + 2.0
        y2[i] = (sig - 1.0)/p
        u[i]  = (6.0 * ((y[i+1] - y[i])/(x[i+1] - x[i]) -(y[i] - y[i-1]) / (x[i] - x[i-1]))/(x[i+1] - x[i-1]) - sig * u[i-1])/p
    end

    # println("u: $u")
    # println("y2: $y2")
    
	if ypn>0.99E30
		qn = 0
		un = 0
    else
		qn = 0.5
		un = (3.0 / (x[n] - x[n-1])) * (ypn - (y[n] - y[n-1])/(x[n] - x[n-1]))
    end

    # println("un: $un")
    
    y2[n] = (un - qn * u[n-1]) / (qn * y2[n-1] + 1.0)

    # println("y2: $y2")
    
    for k in (n-1):-1:1
        y2[k] = y2[k] * y2[k+1] + u[k]
    end

    # println("y2: $y2")

    return y2
end

@inline function zeta(zz::Real, zl::Real; re::Real)
    return ((zz-zl)*(re+zl)/(re+zz))
end

"""
Calculate temperature and density profiles for lower atmosphere

Arguments:
- `alt::Real`
- `dm28m`
- `xnm`
- `mn3`
- `zn3`
- `mn2`
- `zn2`

Returns:
- `n2` N2 number density
- `tz`
- `meso_tn3`
- `meso_tgn3`
- `meso_tn2`
- `meso_tgn2`
"""
function densm!(alt::Real, d0::Real, xm::Real, tz::Real, 
                mn3::Int, zn3::Array{<:Real, 1}, tn3::Array{<:Real, 1}, tgn3::Array{<:Real, 1}, 
                mn2::Int, zn2::Array{<:Real, 1}, tn2::Array{<:Real, 1}, tgn2::Array{<:Real, 1};
                gsurf::Real, re::Real)
    

    # @printf "alt: %.10f xm: %f tz: %.10f\n" alt xm tz

    xs = zeros(Float64, 10)
    ys = zeros(Float64, 10)

    rgas = 831.4

    densm_tmp = d0

    if alt > zn2[1]
        if xm == 0.0
            return tz, tz
        else
            # println("densm return 1")
            return tz, d0
        end
    end

    # Straosphere/Mesosphere Temperature
    if alt > zn2[mn2]
        z = alt
    else
        z = zn2[mn2]
    end
    mn = mn2
    z1 = zn2[1]
    z2 = zn2[mn]
    t1 = tn2[1]
    t2 = tn2[mn]
    zg = zeta(z, z1, re=re)
    zgdif = zeta(z2, z1, re=re)

    # Setup spline nodes
    for k in 1:mn
        xs[k] = zeta(zn2[k], z1, re=re)/zgdif
        ys[k] = 1.0/tn2[k]
    end
    yd1 = -tgn2[1] / (t1*t1) * zgdif
    yd2 = -tgn2[2] / (t2*t2) * zgdif * ((re+z2)/(re+z1))^2

    # Calculate spline coefficients
    y2out = spline(xs, ys, mn, yd1, yd2)
    x     = zg/zgdif
    y     = splint(xs, ys, y2out, mn, x)

    # Temperature at altitude
    tz = 1.0 / y

    # @printf "tz: %.10f densm_tmp: %.10f\n" tz densm_tmp

    if xm != 0
        # Calculate stratosphere/mesosphere density
        glb  = gsurf / (1.0 + z1/re)^2
        gamm = xm * glb * zgdif / rgas

        # Integrate temperature profiles
        yi = splini(xs, ys,y2out, mn, x)
        expl = gamm*yi
        if expl > 50
            expl = 50
        end

        # Density at altitude
        densm_tmp = densm_tmp * (t1 / tz) * exp(-expl)
    end

    if alt > zn3[1]
        if xm == 0.0
            return tz, tz
        else
            # println("densm return 2")
            # @printf "tz: %.10f densm_tmp: %.10f\n" tz densm_tmp
            return tz, densm_tmp
        end
    end

    # Troposphere / stratosphere temperatures 
    z  = alt
    mn = mn3
    z1 = zn3[1]
    z2 = zn3[mn]
    t1 = tn3[1]
    t2 = tn3[mn]
    zg = zeta(z, z1, re=re)
    zgdif = zeta(z2, z1, re=re)

    # @printf "mn: %d z1: %f z2: %f t1: %f t2: %f zg: %f zgdif: %f\n" mn z1 z2 t1 t2 zg zgdif

    # Setup spline nodes
    for k in 1:mn
        xs[k] = zeta(zn3[k], z1, re=re)/zgdif
        ys[k] = 1.0/tn3[k]
    end
    yd1 = -tgn3[1] / (t1*t1) * zgdif
    yd2 = -tgn3[2] / (t2*t2) * zgdif * ((re+z2)/(re+z1))^2

    # Calculate spline coefficients
    # println("mn: $mn")
    # println("xs: $xs")
    # println("ys: $ys")
    # println("yd1: $yd1")
    # println("yd2: $yd2")
    y2out = spline(xs, ys, mn, yd1, yd2)
    # println("y2out: $y2out")
    x = zg/zgdif
    y = splint(xs, ys, y2out, mn, x)
    # println("x: $x")
    # println("y: $y")

    # @printf "tz: %.10f densm_tmp: %.10f\n" tz densm_tmp

    # Temperature at altitude
    tz = 1.0 / y

    # @printf "tz: %.10f densm_tmp: %.10f\n" tz densm_tmp

    if xm != 0
        # Calculate stratosphere/mesosphere density
        glb  = gsurf / (1.0 + z1/re)^2
        gamm = xm * glb * zgdif / rgas

        # Integrate temperature profiles
        yi = splini(xs, ys, y2out, mn, x)
        expl = gamm*yi
        if expl > 50
            expl = 50
        end

        # Density at altitude
        densm_tmp = densm_tmp * (t1 / tz) * exp(-expl)
    end

    if xm == 0.0
        return tz, tz
    else
        # println("densm return 3")
        # @printf "tz: %.10f densm_tmp: %.10f\n" tz densm_tmp
        return tz, densm_tmp
    end    
end

function densu!(alt::Real, dlb::Real, tinf::Real, tlb::Real, xm::Real, 
               alpha::Real, tz::Real, zlb::Real, s2::Real, 
               mn1::Int, zn1::Array{<:Real, 1}, tn1::Array{<:Real, 1}, tgn1::Array{<:Real, 1};
               gsurf::Real, re::Real)

    x          = 0
    rgas       = 831.4
    densu_temp = 1.0
    z1         = 0
    t1         = 0
    zgdif      = 0

    xs    = zeros(Float64, 5)
    ys    = zeros(Float64, 5)
    y2out = zeros(Float64, 5)
    
    # Join altitudes of Bates and spline
    za = zn1[1]
    if alt > za
        z = alt
    else
        z = za
    end

    # println("z: $z, zlb: $zlb, tinf: $tinf, tlb: $tlb")

    # Geopotential altitude difference from ZLB
    zg2 = zeta(z, zlb, re=re)

    # Bates temperatures
    tt = tinf - (tinf - tlb) * exp(-s2*zg2)
    ta = tt
    tz = tt
    densu_temp = tz

    # println("alt: $alt, za: $za")
    if alt < za 
        # Calculate temperature below za
        # Get temperaturre gradient at za from Bates profile
        dta = (tinf - ta) * s2 * ((re+zlb)/(re+za))^2

        tgn1[1] = dta
        tn1[1]  = ta

        if alt > zn1[end]
            z = alt
        else
            z = zn1[end]
        end

        mn = mn1
        z1 = zn1[1]
        z2 = zn1[end]
        t1 = tn1[1]
        t2 = tn1[end]

        # Geopotential difference from z1
        zg    = zeta(z, z1, re=re)
        zgdif = zeta(z2, z1, re=re)

        # Setup spline nodes
        for k in 1:mn
            xs[k] = zeta(zn1[k], z1, re=re)/zgdif
            ys[k] = 1.0/tn1[k]
        end

        # End node derivatives
        yd1 = -tgn1[1] / (t1*t1) * zgdif
        yd2 = -tgn1[2] / (t2*t2) * zgdif * ((re+z2)/(re+z1))^2

        # Calculate spline coefficients
        y2out[5] = 0.0
        # println("xs: $xs")
        # println("ys: $ys")
        # println("mn: $mn")
        # println("yd1: $yd1")
        # println("yd2: $yd2")
        # println("before spline call")
        y2out = spline(xs, ys, mn, yd1, yd2)
        x     = zg/zgdif
        y     = splint(xs, ys, y2out, mn, x)

        # Temperature at altitude
        tz = 1.0 / y
        densu_temp = tz 
        # println("tz: $tz, densu_temp: $densu_temp")
    end

    # println("break 1")
    # println("tz: $tz, densu_temp: $densu_temp")

    if xm == 0
        # println("return 1")
        # println("tz: $tz, densu_temp: $densu_temp")
        return tz, densu_temp
    end

    # Calculate density above za
    glb = gsurf / (1.0+zlb/re)^2
    gamma = xm *glb / (s2 * rgas * tinf)
    expl = exp(-s2 * gamma * zg2)

    if expl > 50
        expl = 50
    end

    if tt <= 0
        expl = 50
    end

    # println("break 2")
    # println("tz: $tz, densu_temp: $densu_temp")

    # @printf "dlb: %.10f tlb: %.10f tt: %.10f alpha: %.10f gamma: %.10f expl: %.10f\n" dlb tlb tt alpha gamma expl

    # density at altitude
    densa      = dlb * (tlb/tt)^(1.0+alpha+gamma)*expl
    densu_temp = densa

    # println("break 3")
    # println("tz: $tz, densu_temp: $densu_temp")


    if alt >= za
        # println("return 2")
        # println("tz: $tz, densu_temp: $densu_temp")
        return tz, densu_temp
    end

    # println("branch here")
    # println("tz: $tz, densu_temp: $densu_temp")

    # Calculate density below za
    glb  = gsurf / (1.0 + z1/re)^2.0
    gamm = xm * glb * zgdif / rgas 

    # Integrate spline temperatures
    yi = splini(xs, ys, y2out, mn, x)
    # println("xs: $xs")
    # println("ys: $ys")
    # println("y2out: $y2out")
    # println("mn: $mn")
    # println("x: $x")

    expl = gamm * yi
    # @printf "%.10f\n" yi
    # println("glb: $glb, gamm: $gamm, yi: $yi, expl: $expl")

    if expl > 50.0
        expl = 50.0
    end

    if tz <= 0
        expl = 50.0
    end

    # density at altitude
    # println("densu_temp: $densu_temp, t1: $t1, tz: $tz, alpha: $alpha, expl: $expl")
    densu_temp = densu_temp * (t1/tz)^(1.0+alpha) * exp(-expl)
    # println("return 3")
    # println("tz: $tz, densu_temp: $densu_temp")
    return tz, densu_temp
end

# Equation A24a
@inline function g0(a::Real, p::Array{<:Real, 1})
    return (a - 4.0 + (p[26] - 1.0) * (a - 4.0 + (exp(-sqrt(p[25]*p[25]) * (a - 4.0)) - 1.0) / sqrt(p[25]*p[25])))
end

# Equation A24c
@inline function sumex(ex::Real)
    return (1.0 + (1.0 - ex^19.0) / (1.0 - ex) * ex^0.5)
end

# Equation A24a
@inline function sg0(ex::Real, p::Array{<:Real, 1}, ap::Array{<:Real, 1})
    return (g0(ap[2],p) + (g0(ap[3],p)*ex + g0(ap[4],p)*ex*ex + 
            g0(ap[5],p)*ex^3.0	+ (g0(ap[6],p)*ex^4.0 + 
            g0(ap[7],p)*ex^12.0)*(1.0- ex^8.0)/(1.0-ex)))/sumex(ex)
end

##################
# NRLMSISE Model #
##################

"""
Calculate G(L) function
"""
function globe7(p::Array{<:Real, 1}, input::NRLMSISE_Input, flags::NRLMSISE_Flags)
    # Working variables
    t = zeros(Float64, 15)

    # LPOLY Variables
    dfa = 0.0
    plg = Array{Float64, 1}[zeros(Float64, 9) for i in 1:4]
    ctloc,  stloc  = 0.0, 0.0
    c2tloc, s2tloc = 0.0, 0.0
    s3tloc, c3tloc = 0.0, 0.0
    apdf = 0.0
    apt  = zeros(Float64, 4)

    # Upper therrmosphere parameters
    sr   = 7.2722E-5
    dgtr = 1.74533E-2
    dr   = 1.72142E-2
    hr   = 0.2618

    tloc = input.lst
    for j in 1:14
        t[j] = 0
    end

    # Calculate legendre polynomials
    c  = sin(input.g_lat * dgtr)
    s  = cos(input.g_lat * dgtr)
    c2 = c*c
    c4 = c2*c2
    s2 = s*s

    plg[1][2] = c
	plg[1][3] = 0.5*(3.0*c2 -1.0)
	plg[1][4] = 0.5*(5.0*c*c2-3.0*c)
	plg[1][5] = (35.0*c4 - 30.0*c2 + 3.0)/8.0
	plg[1][6] = (63.0*c2*c2*c - 70.0*c2*c + 15.0*c)/8.0
	plg[1][7] = (11.0*c*plg[1][6] - 5.0*plg[1][5])/6.0
    # plg[1][8] = (13.0*c*plg[1][7] - 6.0*plg[1][6])/7.0
	plg[2][2] = s
	plg[2][3] = 3.0*c*s
	plg[2][4] = 1.5*(5.0*c2-1.0)*s
	plg[2][5] = 2.5*(7.0*c2*c-3.0*c)*s
	plg[2][6] = 1.875*(21.0*c4 - 14.0*c2 +1.0)*s
	plg[2][7] = (11.0*c*plg[2][6]-6.0*plg[2][5])/5.0
    # plg[2][8] = (13.0*c*plg[2][7]-7.0*plg[2][6])/6.0
    # plg[2][9] = (15.0*c*plg[2][8]-8.0*plg[2][7])/7.0
	plg[3][3] = 3.0*s2
	plg[3][4] = 15.0*s2*c
	plg[3][5] = 7.5*(7.0*c2 -1.0)*s2
	plg[3][6] = 3.0*c*plg[3][5]-2.0*plg[3][4]
	plg[3][7] = (11.0*c*plg[3][6]-7.0*plg[3][5])/4.0
	plg[3][8] = (13.0*c*plg[3][7]-8.0*plg[3][6])/5.0
	plg[4][4] = 15.0*s2*s
	plg[4][5] = 105.0*s2*s*c 
	plg[4][6] = (9.0*c*plg[4][5]-7.0*plg[4][4])/2.0
	plg[4][7] = (11.0*c*plg[4][6]-8.0*plg[4][5])/3.0

    if !(((flags.sw[8] == 0) && (flags.sw[9] == 0)) && (flags.sw[15] == 0))
        stloc  = sin(hr*tloc)
        ctloc  = cos(hr*tloc)
        s2tloc = sin(2.0*hr*tloc)
        c2tloc = cos(2.0*hr*tloc)
        s3tloc = sin(3.0*hr*tloc)
        c3tloc = cos(3.0*hr*tloc)
    end

    cd32 = cos(dr*(input.doy - p[32]))
    cd18 = cos(2.0*dr*(input.doy - p[18]))
    cd14 = cos(dr*(input.doy - p[14]))
    cd39 = cos(2.0*dr*(input.doy - p[39]))

    # F10.7 Effect
    df   = input.f107 - input.f107A
    dfa = input.f107A - 150.0
    t[1] = p[20]*df*(1.0+p[60]*dfa) + p[21]*df*df + p[22]*dfa + p[30]*dfa^2.0
    f1   = 1.0 + (p[48]*dfa +p[20]*df+p[21]*df*df)*flags.swc[2]
    f2   = 1.0 + (p[50]*dfa+p[20]*df+p[21]*df*df)*flags.swc[2]

    # Time independent
    t[2] = (p[2]*plg[1][3]+ p[3]*plg[1][5]+p[23]*plg[1][7]) +
           (p[15]*plg[1][3])*dfa*flags.swc[2] +p[27]*plg[1][2]
      
    # Symmetrical annual
    t[3] = p[19]*cd32

    # Symmetrical semiannual
    t[4] =  (p[16]+p[17]*plg[1][3])*cd18

    # Asymmetrical annual
    t[5] =  f1*(p[10]*plg[1][2]+p[11]*plg[1][4])*cd14

    # Asymmetrical semiannual
    t[6] = p[38]*plg[1][2]*cd39

    # Diurnal
    if flags.sw[7] != 0
        t71  = (p[12]*plg[2][3])*cd14*flags.swc[6]
        t72  = (p[13]*plg[2][3])*cd14*flags.swc[6]
        t[7] = (f2*((p[4]*plg[2][2] + p[5]*plg[2][4] + p[28]*plg[2][6] + t71) * 
               ctloc + (p[7]*plg[2][2] + p[8]*plg[2][4] + p[29]*plg[2][6] + t72)*stloc))
    end     

    # Semidiurnal
    if flags.sw[9] != 0
        t81  = (p[24]*plg[3][4]+p[36]*plg[3][6])*cd14*flags.swc[6]
		t82  = (p[34]*plg[3][4]+p[37]*plg[3][6])*cd14*flags.swc[6]
        t[8] = (f2*((p[6]*plg[3][3]+ p[42]*plg[3][5] + t81)*c2tloc + 
               (p[9]*plg[3][3] + p[43]*plg[3][5] + t82)*s2tloc))
    end

    # Terdiurnal
    if flags.sw[15] != 0
        t[14] = (f2*((p[40]*plg[4][4]+(p[94]*plg[4][5]+p[47]*plg[4][7])*cd14*flags.swc[6])*s3tloc
                 + (p[41]*plg[4][4]+(p[95]*plg[4][5]+p[49]*plg[4][7])*cd14*flags.swc[6])*c3tloc))
    end

    # Magnetic activity based on daily AP
    if flags.sw[10] == -1
        if p[52] != 0
            exp1 = exp(-10800.0*sqrt(p[52]*p[52])/(1.0+p[139]*(45.0-sqrt(input.g_lat*input.g_lat))))
        end

        if exp1 > 0.99999
            exp1 = 0.99999
        end

        if p[25] < 1.0e-4
            p[25] = 1.0e-4
        end

        apt[1] = sg0(exp1, p, input.ap_array)

        if flags.sw[10] != 0
            t[8] = (apt[1]*(p[51]+p[97]*plg[1][3]+p[55]*plg[1][5] + 
                    (p[126]*plg[1][2]+p[127]*plg[1][4]+p[128]*plg[1][6])*cd14*flags.swc[6] +
                    (p[129]*plg[2][2]+p[130]*plg[2][4]+p[131]*plg[2][6])*flags.swc[8] * 
                    cos(hr*(tloc-p[132]))))
        end
    else
        apd = input.ap - 4.0
        p44 = p[44]
        p45 = p[45]
        if p44 < 0
            p44 = 1.0e-5
        end
        apdf = apd + (p45-1.0)*(apd + (exp(-p44 * apd) - 1.0)/p44)
        if flags.sw[10] != 0
            t[9] = (apdf*(p[33]+p[46]*plg[1][3]+p[35]*plg[1][5] +
                    (p[101]*plg[1][2]+p[102]*plg[1][4]+p[103]*plg[1][6])*cd14*flags.swc[6] +
                    (p[122]*plg[2][2]+p[123]*plg[2][4]+p[124]*plg[2][6])*flags.swc[8]*
                    cos(hr*(tloc-p[125]))))
        end
    end

    if flags.sw[11] != 0 && (input.g_lon > - 1000.0) != 0
        # Longitudinal
        if flags.sw[12] != 0
			t[11] = ((1.0 + p[81]*dfa*flags.swc[2])*
                    ((p[65]*plg[2][3]+p[66]*plg[2][5]+p[67]*plg[2][7]
                    + p[104]*plg[2][2]+p[105]*plg[2][4]+p[106]*plg[2][6]
                    + flags.swc[6]*(p[110]*plg[2][2]+p[111]*plg[2][4]+p[112]*plg[2][6])*cd14)*
                        cos(dgtr*input.g_lon)
                    + (p[91]*plg[2][3]+p[92]*plg[2][5]+p[93]*plg[2][7]
                    + p[107]*plg[2][2]+p[108]*plg[2][4]+p[109]*plg[2][6]
                    + flags.swc[6]*(p[113]*plg[2][2]+p[114]*plg[2][4]+p[115]*plg[2][6])*cd14)*
                    sin(dgtr*input.g_lon)))
        end

        # UT and mixed UT, longitude
        if flags.sw[13] != 0
            t[12] = ((1.0+p[96]*plg[1][2])*(1.0+p[82]*dfa*flags.swc[2])*
				    (1.0+p[120]*plg[1][2]*flags.swc[6]*cd14)*
				    ((p[69]*plg[1][2]+p[70]*plg[1][4]+p[71]*plg[1][6])*
                    cos(sr*(input.sec-p[72]))))
                    
			t[12] += (flags.swc[12]*(p[77]*plg[3][4]+p[78]*plg[3][6]+p[79]*plg[3][8])*
				    cos(sr*(input.sec-p[80])+2.0*dgtr*input.g_lon)*(1.0+p[138]*dfa*flags.swc[2]))
        end

        # UT longitude magnetic activity
        if flags.sw[14] != 0
            if flags.sw[10] == -1
                if p[52] != 0.0
                    t[13] = (apt[1]*flags.swc[12]*(1.0+p[133]*plg[1][2])*
                            ((p[53]*plg[2][3]+p[99]*plg[2][5]+p[68]*plg[2][7])*
                            cos(dgtr*(input.g_lon-p[98])))
                            +apt[1]*flags.swc[12]*flags.swc[6]*
                            (p[134]*plg[2][2]+p[135]*plg[2][4]+p[136]*plg[2][6])*
                            cd14*cos(dgtr*(input.g_lon-p[137])) 
                            +apt[1]*flags.swc[13]* 
                            (p[56]*plg[1][2]+p[57]*plg[1][4]+p[58]*plg[1][6])*
                            cos(sr*(input.sec-p[59])))
                end
            else
                t[13] = (apdf*flags.swc[12]*(1.0+p[121]*plg[1][2])*
					((p[61]*plg[2][3]+p[62]*plg[2][5]+p[63]*plg[2][7])*
					cos(dgtr*(input.g_lon-p[64])))
					+ apdf*flags.swc[12]*flags.swc[6] * 
					(p[116]*plg[2][2]+p[117]*plg[2][4]+p[118]*plg[2][6])* 
					cd14*cos(dgtr*(input.g_lon-p[119])) 
					+ apdf*flags.swc[13]* 
					(p[84]*plg[1][2]+p[85]*plg[1][4]+p[86]*plg[1][6])* 
                    cos(sr*(input.sec-p[76])))
            end
        end
        

    end

    # Params not used: 82, 89, 99, 139-149
    tinf = p[31]
    for i in 1:14
        tinf = tinf + abs(flags.sw[i+1])*t[i]
    end

    # println("t: $t")
    # println("tinf: $tinf dfa: $dfa ctloc: $ctloc stloc: $stloc c2tloc: $c2tloc s2tloc: $s2tloc c3tloc: $c3tloc s3tloc: $s3tloc apdf: $apdf")

    return tinf, dfa, plg, ctloc, stloc, c2tloc, s2tloc, c3tloc, s3tloc, apdf, apt
end

"""
Calculate G(L) function for lower atmosphere
"""
function glob7s(p::Array{<:Real, 1}, input::NRLMSISE_Input, flags::NRLMSISE_Flags; 
                dfa::Real, plg::Array{<:Array{<:Real,1},1}, ctloc::Real, stloc::Real, c2tloc::Real, s2tloc::Real, s3tloc::Real, c3tloc::Real)
    # Working variables
    pset = 2.0
    t    = zeros(Float64, 14)
    dr   = 1.72142E-2
    dgtr = 1.74533E-2

    # Confirm parameter set 
    if p[100] == 0
        p[100] = pset
    end

    if p[100] != pset
        @error "Wrong parameter set for glob7s"
        return -1
    end

    for j in 1:14
        t[j] = 0.0
    end

	cd32 = cos(dr*(input.doy - p[32]))
	cd18 = cos(2.0*dr*(input.doy - p[18]))
	cd14 = cos(dr*(input.doy - p[14]))
	cd39 = cos(2.0*dr*(input.doy - p[39]))

    # F10.7
    t[1] = p[22]*dfa

    # Time independent
    t[2] = p[2]*plg[1][3] + p[3]*plg[1][5] + p[23]*plg[1][7] + p[27]*plg[1][2] + p[15]*plg[1][4] + p[60]*plg[1][6]

    # Symmetrical annual
    t[3] = (p[19]+p[48]*plg[1][3]+p[30]*plg[1][5])*cd32

    # Symmetrical semi annual
    t[4]=(p[16]+p[17]*plg[1][3]+p[31]*plg[1][5])*cd18

    # Asymmetrical annual
    t[5] = (p[10]*plg[1][2]+p[11]*plg[1][4]+p[21]*plg[1][6])*cd14

    # Asymmetrical semiannual
    t[6] = (p[38]*plg[1][2])*cd39

    # Diurnal
    if flags.sw[8] != 0
        t71  = p[12]*plg[2][3]*cd14*flags.swc[6]
		t72  = p[13]*plg[2][3]*cd14*flags.swc[6]
		t[7] = ((p[4]*plg[2][2] + p[5]*plg[2][4] + t71) * ctloc + (p[7]*plg[2][2] + p[8]*plg[2][4] + t72) * stloc)
    end

    # Semidiurnal
    if flags.sw[9] != 0
        t81  = (p[24]*plg[3][4]+p[36]*plg[3][6])*cd14*flags.swc[6]
		t82  = (p[34]*plg[3][4]+p[37]*plg[3][6])*cd14*flags.swc[6]
        t[8] = ((p[6]*plg[3][3] + p[42]*plg[3][5] + t81) * c2tloc + (p[9]*plg[3][3] + p[43]*plg[3][5] + t82) * s2tloc)
    end

    # Terdiurnal
    if flags.sw[15] != 0
        t[14] = p[40] * plg[4][4] * s3tloc + p[41] * plg[4][4] * c3tloc
    end

    # Magnetic activity
    if !((flags.sw[11]==0) || (flags.sw[12]==0) || (input.g_lon <= -1000.0))
        t[11] = (1.0 + plg[1][2]*(p[81]*flags.swc[6]*cos(dr*(input.doy-p[82]))
		        + p[86]*flags.swc[7]*cos(2.0*dr*(input.doy-p[87])))
			    + p[84]*flags.swc[4]*cos(dr*(input.doy-p[85]))
			    + p[88]*flags.swc[5]*cos(2.0*dr*(input.doy-p[89])))*((p[65]*plg[2][3]+p[66]*plg[2][5]+p[67]*plg[2][7]
			    + p[75]*plg[2][2]+p[76]*plg[2][4]+p[77]*plg[2][6])*cos(dgtr*input.g_lon)
			    + (p[91]*plg[2][3]+p[92]*plg[2][5]+p[93]*plg[2][7]
			    + p[78]*plg[2][2]+p[79]*plg[2][4]+p[80]*plg[2][6])*sin(dgtr*input.g_lon))
    end


    tt = 0
    for i in 1:14
        tt += abs(flags.sw[i+1])*t[i]
    end

    return tt
end
!
function gtd7!(input::NRLMSISE_Input, flags::NRLMSISE_Flags, output::NRLMSISE_Output)
    # Working output varriable 
    soutput = NRLMSISE_Output()
    tz = 0.0
    dm28m = 0.0

    # MESO7
    meso_tn1  = zeros(Float64, 5)
    meso_tn2  = zeros(Float64, 4)
    meso_tn3  = zeros(Float64, 5)
    meso_tgn1 = zeros(Float64, 2)
    meso_tgn2 = zeros(Float64, 2)
    meso_tgn3 = zeros(Float64, 2)

    # Working variables
    mn3  = 5
    mn2  = 4
    zn3  = [32.5, 20.0, 15.0, 10.0, 0.0]
    zn2  = [72.5, 55.0, 45.0, 32.5]
    zmix = 62.5

    # Set Input 
    tselec!(flags)

    # Lattitude variation of gravity (none for flags.sw[3] = 0)
    xlat = input.g_lat
    if flags.sw[3] == 0
        xlat = 45.0
    end

    gsurf, re = glatf(xlat)

    xmm = pdm[3][5]

    # Thermosphere / mesosphere (above zn2[1])
    altt = 0.0
    if input.alt > zn2[1]
        altt = input.alt
    else
        altt = zn2[1]
    end

    tmp = input.alt
    input.alt = altt

    dm28, meso_tn1, meso_tgn1, dfa, plg, ctloc, stloc, c2tloc, s2tloc, s3tloc, c3tloc = gts7!(input, flags, soutput, gsurf=gsurf, re=re)

    # println(soutput.d)
    # println(soutput.t)

    altt = input.alt 
    input.alt = tmp

    if flags.sw[1] != 0 # Metric adjustment != 0
        dm28m = dm28*1.0e6
    else
        dm28m = dm28
    end

    output.t[1] = soutput.t[1]
    output.t[2] = soutput.t[2]

    if input.alt >= zn2[1]
        for i in 1:9
            output.d[i] = soutput.d[i]
        end
        # for i in 1:8
        #     @printf "output.d[%d]: %.10f\n" i output.d[i]
        # end
        return
    end

    # Low Mesosphere/Upper stratosphere between zn3[1] and zn2[1]
    meso_tgn2[1] = meso_tgn1[2]
    meso_tn2[1]  = meso_tn1[5]
    meso_tn2[2]  = pma[1][1]*pavgm[1]/(1.0-flags.sw[21]*glob7s(pma[1], input, flags, dfa=dfa, plg=plg, ctloc=ctloc, stloc=stloc, c2tloc=c2tloc, s2tloc=s2tloc, s3tloc=s3tloc, c3tloc=c3tloc))
    meso_tn2[3]  = pma[2][1]*pavgm[2]/(1.0-flags.sw[21]*glob7s(pma[2], input, flags, dfa=dfa, plg=plg, ctloc=ctloc, stloc=stloc, c2tloc=c2tloc, s2tloc=s2tloc, s3tloc=s3tloc, c3tloc=c3tloc))
    meso_tn2[4]  = pma[3][1]*pavgm[3]/(1.0-flags.sw[21]*flags.sw[23]*glob7s(pma[3], input, flags, dfa=dfa, plg=plg, ctloc=ctloc, stloc=stloc, c2tloc=c2tloc, s2tloc=s2tloc, s3tloc=s3tloc, c3tloc=c3tloc))
    meso_tgn2[2] = pavgm[9]*pma[10][1]*(1.0+flags.sw[21]*flags.sw[23]*glob7s(pma[10], input, flags, dfa=dfa, plg=plg, ctloc=ctloc, stloc=stloc, c2tloc=c2tloc, s2tloc=s2tloc, s3tloc=s3tloc, c3tloc=c3tloc))*meso_tn2[4]*meso_tn2[4]/((pma[3][1]*pavgm[3])^2.0)
    meso_tn3[1]  = meso_tn2[4]
    
    # Lower stratrosphere and troposphere below (zn3[1])
    if input.alt < zn3[1]
        meso_tgn3[1] = meso_tgn2[2]
        meso_tn3[2]  = pma[4][1]*pavgm[4]/(1.0-flags.sw[23]*glob7s(pma[4], input, flags, dfa=dfa, plg=plg, ctloc=ctloc, stloc=stloc, c2tloc=c2tloc, s2tloc=s2tloc, s3tloc=s3tloc, c3tloc=c3tloc))
        meso_tn3[3]  = pma[5][1]*pavgm[5]/(1.0-flags.sw[23]*glob7s(pma[5], input, flags, dfa=dfa, plg=plg, ctloc=ctloc, stloc=stloc, c2tloc=c2tloc, s2tloc=s2tloc, s3tloc=s3tloc, c3tloc=c3tloc))
        meso_tn3[4]  = pma[6][1]*pavgm[6]/(1.0-flags.sw[23]*glob7s(pma[6], input, flags, dfa=dfa, plg=plg, ctloc=ctloc, stloc=stloc, c2tloc=c2tloc, s2tloc=s2tloc, s3tloc=s3tloc, c3tloc=c3tloc))
        meso_tn3[5]  = pma[7][1]*pavgm[7]/(1.0-flags.sw[23]*glob7s(pma[7], input, flags, dfa=dfa, plg=plg, ctloc=ctloc, stloc=stloc, c2tloc=c2tloc, s2tloc=s2tloc, s3tloc=s3tloc, c3tloc=c3tloc))
        meso_tgn3[2] = pma[8][1]*pavgm[8]*(1.0+flags.sw[23]*glob7s(pma[8], input, flags, dfa=dfa, plg=plg, ctloc=ctloc, stloc=stloc, c2tloc=c2tloc, s2tloc=s2tloc, s3tloc=s3tloc, c3tloc=c3tloc))*meso_tn3[5]*meso_tn3[5]/((pma[7][1]*pavgm[7])^2.0)
    end

    # println("meso_tn2: $meso_tn2")
    # println("meso_tn3: $meso_tn3")
    # println("meso_tgn2: $meso_tgn2")
    # println("meso_tgn3: $meso_tgn3")

    # Linear transition to full mixing below zn2[0]
    dmc = 0
    if input.alt > zmix
        dmc = 1.0 - (zn2[1]-input.alt)/(zn2[1] - zmix)
    end
    dz28 = soutput.d[3]

    # N2 Density
    dmr = soutput.d[3] / dm28m - 1.0
    tz, output.d[3] = densm!(input.alt, dm28m, xmm, tz, 
                        mn3, zn3, meso_tn3, meso_tgn3, 
                        mn2, zn2, meso_tn2, meso_tgn2,
                        gsurf=gsurf, re=re)
    output.d[3] = output.d[3] * (1.0 + dmr*dmc)

    # HE density
    dmr         = soutput.d[1] / (dz28 * pdm[1][2]) - 1.0
    output.d[1] = output.d[3] * pdm[1][2] * (1.0 + dmr*dmc)

    # O density 
    output.d[2] = 0
    output.d[9] = 0

    # O2 density
    dmr = soutput.d[4] / (dz28 * pdm[4][2]) - 1.0
    output.d[4] = output.d[3] * pdm[4][2] * (1.0 + dmr*dmc)

    # AR density
    dmr = soutput.d[5] / (dz28 * pdm[5][2]) - 1.0
    output.d[5] = output.d[3] * pdm[5][2] * (1.0 + dmr*dmc)

    # H density
    output.d[7] = 0

    # Atomic N density
    output.d[8] = 0

    # Total mass density
    output.d[6] = 1.66e-24 * (  4.0 * output.d[1] 
                             + 16.0 * output.d[2] 
                             + 28.0 * output.d[3]
                             + 32.0 * output.d[4]
                             + 40.0 * output.d[5]
                             +  1.0 * output.d[7]
                             + 14.0 * output.d[8])

    # Correct units
    if flags.sw[1] != 0
        output.d[6] = output.d[6]/1000
    end

    # for i in 1:8
    #     @printf "output.d[%d]: %.10f\n" i output.d[i]
    # end

    # Temperature at altitude
    tz, dd = densm!(input.alt, 1.0, 0, tz, 
                    mn3, zn3, meso_tn3, meso_tgn3, 
                    mn2, zn2, meso_tn2, meso_tgn2,
                    gsurf=gsurf, re=re)
    output.t[2] = tz
end

function gtd7d!(input::NRLMSISE_Input, flags::NRLMSISE_Flags, output::NRLMSISE_Output)
    gtd7!(input, flags, output)

    output.d[6] = 1.66e-24 * (  4.0 * output.d[1] 
                             + 16.0 * output.d[2] 
                             + 28.0 * output.d[3]
                             + 32.0 * output.d[4]
                             + 40.0 * output.d[5]
                             +  1.0 * output.d[7]
                             + 14.0 * output.d[8])
    
    if flags.sw[1] != 0
        output.d[6] = output.d[6]/1000
    end
end

"""
Thermospheric portion of NRLMSISE-00.
See GTD7 for more extensive comments
For alt > 72.5 km
"""
function gts7!(input::NRLMSISE_Input, flags::NRLMSISE_Flags, output::NRLMSISE_Output; gsurf::Real, re::Real)
    dm28  = 0.0
    tz    = 0.0
    mn1   = 5
    zn1   = [120.0, 110.0, 100.0, 90.0, 72.5]
    dgtr  = 1.74533E-2
    dr    = 1.72142E-2
    alpha = [-0.38, 0.0, 0.0, 0.0, 0.17, 0.0, -0.38, 0.0, 0.0]
    altl  = [200.0, 300.0, 160.0, 250.0, 240.0, 450.0, 320.0, 450.0]

    # MESO7
    meso_tn1  = zeros(Float64, 5)
    meso_tgn1 = zeros(Float64, 2)

    za     = pdl[2][16]
    zn1[1] = za

    for j in 1:9
        output.d[j] = 0.0
    end

    # Tinf variations not important below za or zn[1]
    if input.alt > zn1[1]
        tinf_g7, dfa, plg, ctloc, stloc, c2tloc, s2tloc, c3tloc, s3tloc, apdf, apt = globe7(pt,input,flags)
        # println("globe7 temp: $tinf_g7")
        tinf = ptm[1]*pt[1] * (1.0+flags.sw[17]*tinf_g7)
    else
        tinf = ptm[1]*pt[1]
        # println("branch 2")
    end
    output.t[1] = tinf

    # println("tinf: $tinf")

    # Gradient variations not important below zn1[6]
    if input.alt > zn1[5]
        tinf_g7, dfa, plg, ctloc, stloc, c2tloc, s2tloc, c3tloc, s3tloc, apdf, apt = globe7(ps, input, flags)
        g0 = ptm[4]*ps[1] * (1.0 + flags.sw[20]*tinf_g7)
    else
        g0 = ptm[4]*ps[1]
    end

    tinf_g7, dfa, plg, ctloc, stloc, c2tloc, s2tloc, c3tloc, s3tloc, apdf, apt = globe7(pd[4], input, flags)
    tlb = ptm[2] * (1.0 + flags.sw[18]*tinf_g7)*pd[4][1]
    s   = g0 / (tinf - tlb)

    # Lower thermosphere temp variations now significant for density above 300km
    if input.alt < 300.0
        meso_tn1[2]  = ptm[7]*ptl[1][1]/(1.0-flags.sw[19]*glob7s(ptl[1], input, flags, dfa=dfa, plg=plg, ctloc=ctloc, stloc=stloc, c2tloc=c2tloc, s2tloc=s2tloc, s3tloc=s3tloc, c3tloc=c3tloc))
        meso_tn1[3]  = ptm[3]*ptl[2][1]/(1.0-flags.sw[19]*glob7s(ptl[2], input, flags, dfa=dfa, plg=plg, ctloc=ctloc, stloc=stloc, c2tloc=c2tloc, s2tloc=s2tloc, s3tloc=s3tloc, c3tloc=c3tloc))
        meso_tn1[4]  = ptm[8]*ptl[3][1]/(1.0-flags.sw[19]*glob7s(ptl[3], input, flags, dfa=dfa, plg=plg, ctloc=ctloc, stloc=stloc, c2tloc=c2tloc, s2tloc=s2tloc, s3tloc=s3tloc, c3tloc=c3tloc))
        meso_tn1[5]  = ptm[5]*ptl[4][1]/(1.0-flags.sw[19]*flags.sw[21]*glob7s(ptl[4], input, flags, dfa=dfa, plg=plg, ctloc=ctloc, stloc=stloc, c2tloc=c2tloc, s2tloc=s2tloc, s3tloc=s3tloc, c3tloc=c3tloc))
        meso_tgn1[2] = ptm[9]*pma[9][1]*(1.0+flags.sw[19]*flags.sw[21]*glob7s(pma[9], input, flags, dfa=dfa, plg=plg, ctloc=ctloc, stloc=stloc, c2tloc=c2tloc, s2tloc=s2tloc, s3tloc=s3tloc, c3tloc=c3tloc))*meso_tn1[5]*meso_tn1[5]/((ptm[5]*ptl[4][1])^2.0)
    else
        meso_tn1[2]  = ptm[7]*ptl[1][1]
        meso_tn1[3]  = ptm[3]*ptl[2][1]
        meso_tn1[4]  = ptm[8]*ptl[3][1]
        meso_tn1[5]  = ptm[5]*ptl[4][1]
        meso_tgn1[2] = ptm[9]*pma[9][1]*meso_tn1[5]*meso_tn1[5]/((ptm[5]*ptl[4][1])^2.0)
    end

    # N2 variations factor at Zlb
    tinf_g7, dfa, plg, ctloc, stloc, c2tloc, s2tloc, c3tloc, s3tloc, apdf, apt = globe7(pd[3], input, flags)
    g28 = flags.sw[22]*tinf_g7
    
    # Variation of turbopause height
    zhf         = pdl[2][25]*(1.0+flags.sw[6]*pdl[1][25]*sin(dgtr*input.g_lat)*cos(dr*(input.doy-pt[14])))
	output.t[1] = tinf
	xmm         = pdm[3][5]
    z           = input.alt

    # N2 Density 
    # Diffusive density at Zlb
    db28 = pdm[3][1]*exp(g28)*pd[3][1]


    # Diffusive density at Alt
    output.t[2], output.d[3] = densu!(z, db28, tinf, tlb, 28.0, alpha[3], output.t[2], ptm[6], s, mn1, zn1, meso_tn1, meso_tgn1, gsurf=gsurf, re=re)
    dd = output.d[3]

    # println("output.d[3]: $(output.d[3])")
    # return 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0

    # Turbopause
    zh28  = pdm[3][3]*zhf
    zhm28 = pdm[3][4]*pdl[2][6]
    xmd   = 28.0 - xmm

    # Mixed density at Zlb
    tz, b28 = densu!(zh28, db28, tinf, tlb, xmd, (alpha[3]-1.0), tz,ptm[6], s, mn1, zn1, meso_tn1, meso_tgn1, gsurf=gsurf, re=re)
    if flags.sw[16] != 0 && (z <= altl[3]) != 0
        # Mixed density at alt
        tz, dm28 = densu!(z, b28, tinf, tlb, xmm, alpha[3], tz, ptm[6], s, mn1, zn1, meso_tn1, meso_tgn1, gsurf=gsurf, re=re)
        # Net density at alt
        output.d[3] = dnet(output.d[3], dm28, zhm28, xmm, 28.0)
    end

    # println("output.d[3]: $(output.d[3])")

    # HE density
    # Density variation factor at Zlb
    tinf_g7, dfa, plg, ctloc, stloc, c2tloc, s2tloc, c3tloc, s3tloc, apdf, apt = globe7(pd[1], input, flags)
    g4 = flags.sw[22]*tinf_g7
    
    # Diffusive denity at Zlb
    db04 = pdm[1][1]*exp(g4)*pd[1][1]
    
    # Diffusive density at alt
    tz, output.d[1] = densu!(z, db04, tinf, tlb, 4.0, alpha[1], output.t[2], ptm[6], s, mn1, zn1, meso_tn1, meso_tgn1, gsurf=gsurf, re=re)
    dd = output.d[1]
    if flags.sw[16] != 0 && (z < altl[1]) != 0
        # Turbo pause
        zh04 = pdm[1][3]
        
        # Mixed density at Zlb
        tz, b04 = densu!(zh04, db04, tinf, tlb, 4.0-xmm, alpha[1]-1.0, output.t[2], ptm[6], s, mn1, zn1, meso_tn1, meso_tgn1, gsurf=gsurf, re=re)

        # Mixed density at alt
        tz, dm04  = densu!(z, b04, tinf, tlb, xmm, 0.0, output.t[2], ptm[6], s, mn1, zn1, meso_tn1, meso_tgn1, gsurf=gsurf, re=re)
        zhm04 = zhm28

        # Net density at alt
        output.d[1] = dnet(output.d[1], dm04, zhm04, xmm, 4.0)

        # Correction to specified mixing ratio at ground
        rl = log(b28*pdm[1][2]/b04)
        zc04 = pdm[1][5]*pdl[2][1]
        hc04 = pdm[1][6]*pdl[2][2]

        # Net density corrected at alt
        output.d[1] = output.d[1]*ccor(z, rl, hc04, zc04)
    end

    # O density
    
    #  Density variation factor at Zlb
    tinf_g7, dfa, plg, ctloc, stloc, c2tloc, s2tloc, c3tloc, s3tloc, apdf, apt = globe7(pd[2], input, flags)
	g16 = flags.sw[22]*tinf_g7
    
    #  Diffusive density at Zlb
	db16 = pdm[2][1]*exp(g16)*pd[2][1]
    
    # Diffusive density at Alt
	tz, output.d[2] = densu!(z, db16, tinf, tlb, 16.0, alpha[2],output.t[2], ptm[6], s, mn1, zn1, meso_tn1, meso_tgn1, gsurf=gsurf, re=re)
	dd=output.d[2]
	if flags.sw[16] != 0 && (z <= altl[2])
		# Turbopause
        zh16=pdm[2][3]
        
		# Mixed density at Zlb
        tz, b16 = densu!(zh16, db16, tinf, tlb, 16.0-xmm, (alpha[2]-1.0), output.t[2], ptm[6], s, mn1, zn1, meso_tn1, meso_tgn1, gsurf=gsurf, re=re)
        
		# Mixed density at Alt
		tz, dm16  = densu!(z, b16, tinf, tlb, xmm, 0.0, output.t[2], ptm[6], s, mn1, zn1, meso_tn1, meso_tgn1, gsurf=gsurf, re=re)
        zhm16 = zhm28
        
		# Net density at Alt
		output.d[2]=dnet(output.d[2], dm16, zhm16, xmm, 16.0)
        rl    = pdm[2][2]*pdl[2][17]*(1.0+flags.sw[2]*pdl[1][24]*(input.f107A - 150.0))
        hc16  = pdm[2][6]*pdl[2][4]
        zc16  = pdm[2][5]*pdl[2][3]
        hc216 = pdm[2][6]*pdl[2][5]
		output.d[2] = output.d[2]*ccor2(z, rl, hc16, zc16, hc216)
        
		# Chemistry correction
		hcc16 = pdm[2][8]*pdl[2][14]
		zcc16 = pdm[2][7]*pdl[2][13]
		rc16  = pdm[2][4]*pdl[2][15]
        
		# Net density corrected at Alt
		output.d[2] = output.d[2]*ccor(z, rc16, hcc16, zcc16)
    end

    # O2 Density 
    
    # Density variation factor at Zlb
    tinf_g7, dfa, plg, ctloc, stloc, c2tloc, s2tloc, c3tloc, s3tloc, apdf, apt = globe7(pd[5], input, flags)
    # @printf "globe7: %.18f\n" tinf_g7
	g32 = flags.sw[22]*tinf_g7
    
    # Diffusive density at Zlb
    db32 = pdm[4][1]*exp(g32)*pd[5][1]
    
    # Diffusive density at Alt
	tz, output.d[4] = densu!(z, db32, tinf, tlb, 32.0, alpha[4], output.t[2], ptm[6], s, mn1, zn1, meso_tn1, meso_tgn1, gsurf=gsurf, re=re)
    dd=output.d[4]

    
    if flags.sw[16] != 0
		if z <= altl[4]
			# Turbopause
            zh32 = pdm[4][3]
            
			# Mixed density at Zlb
            tz, b32 = densu!(zh32, db32, tinf, tlb, 32.0-xmm, alpha[4]-1.0, output.t[2], ptm[6], s, mn1, zn1, meso_tn1, meso_tgn1, gsurf=gsurf, re=re)
            
			# Mixed density at Alt
			tz, dm32  = densu!(z, b32, tinf, tlb, xmm, 0.0, output.t[2], ptm[6], s, mn1, zn1, meso_tn1, meso_tgn1, gsurf=gsurf, re=re)
            zhm32 = zhm28
            
			# Net density at Alt
            output.d[4] = dnet(output.d[4], dm32, zhm32, xmm, 32.0)
            
			# Correction to specified mixing ratio at ground
			rl = log(b28*pdm[4][2]/b32)
			hc32 = pdm[4][6]*pdl[2][8]
			zc32 = pdm[4][5]*pdl[2][7]
			output.d[4] = output.d[4]*ccor(z, rl, hc32, zc32)
        end

        
        # Correction for general departure from diffusive equilibrium above Zlb
		hcc32  = pdm[4][8]*pdl[2][23]
		hcc232 = pdm[4][8]*pdl[1][23]
		zcc32  = pdm[4][7]*pdl[2][22]
		rc32   = pdm[4][4]*pdl[2][24]*(1.0+flags.sw[2]*pdl[1][24]*(input.f107A - 150.0))
        
        # Net density corrected at Alt
		output.d[4]=output.d[4]*ccor2(z, rc32, hcc32, zcc32, hcc232)
    end


    # AR density
    # Density variation factor at Zlb
    tinf_g7, dfa, plg, ctloc, stloc, c2tloc, s2tloc, c3tloc, s3tloc, apdf, apt = globe7(pd[6], input, flags)
	g40 = flags.sw[22]*tinf_g7
    
    # Diffusive density at Zlb
	db40 = pdm[5][1]*exp(g40)*pd[6][1]
    
    # Diffusive density at Alt
	tz, output.d[5] = densu!(z, db40, tinf, tlb, 40.0, alpha[5], output.t[2], ptm[6], s, mn1, zn1, meso_tn1, meso_tgn1, gsurf=gsurf, re=re)
    dd = output.d[5]
    
    # println("AR calc begin")
    # println("output.d[5]: $(output.d[5])")
    
    if flags.sw[16] != 0 && (z <= altl[5]) != 0
		# Turbopause
		zh40 = pdm[5][3]
        
        # Mixed density at Zlb
		tz, b40 = densu!(zh40, db40, tinf, tlb, 40.0-xmm, alpha[5]-1.0, output.t[2], ptm[6], s, mn1, zn1, meso_tn1, meso_tgn1, gsurf=gsurf, re=re)
        
        # Mixed density at Alt
		tz, dm40  = densu!(z, b40, tinf, tlb, xmm, 0.0, output.t[2], ptm[6], s, mn1, zn1, meso_tn1, meso_tgn1, gsurf=gsurf, re=re)
		zhm40 = zhm28
        
        # Net density at Alt
        output.d[5] = dnet(output.d[5], dm40, zhm40, xmm,40.0)
        # println("output.d[5]: $(output.d[5])")
        
        # Correction to specified mixing ratio at ground
		rl   = log(b28*pdm[5][2]/b40)
		hc40 = pdm[5][6]*pdl[2][10]
		zc40 = pdm[5][5]*pdl[2][9]
        
        # Net density corrected at Alt
        output.d[5] = output.d[5]*ccor(z, rl, hc40, zc40)
        # println("output.d[5]: $(output.d[5])")
    end

    # println("AR calc end")

    # Hydrogen Density
    # Density variation factor at Zlb
    tinf_g7, dfa, plg, ctloc, stloc, c2tloc, s2tloc, c3tloc, s3tloc, apdf, apt = globe7(pd[7], input, flags)
    g1 = flags.sw[22]*tinf_g7
    
    # Diffusive density at Zlb
    db01 = pdm[6][1]*exp(g1)*pd[7][1]
    
    # Diffusive density at Alt
	tz, output.d[7] = densu!(z, db01, tinf, tlb, 1.0, alpha[7], output.t[2], ptm[6], s, mn1, zn1, meso_tn1, meso_tgn1, gsurf=gsurf, re=re)
    dd = output.d[7]
    
	if flags.sw[16] != 0 && ( z <= altl[7])
		# Turbopause
		zh01 = pdm[6][3]
        
        # Mixed density at Zlb
		tz, b01 = densu!(zh01, db01, tinf, tlb, 1.0-xmm, alpha[7]-1.0, output.t[2], ptm[6], s, mn1, zn1, meso_tn1, meso_tgn1, gsurf=gsurf, re=re)
        
        # Mixed density at Alt
		tz, dm01  = densu!(z, b01, tinf, tlb, xmm, 0.0, output.t[2], ptm[6], s, mn1, zn1, meso_tn1, meso_tgn1, gsurf=gsurf, re=re)
		zhm01 = zhm28
        
        # Net density at Alt
		output.d[7] = dnet(output.d[7], dm01, zhm01, xmm, 1.0)
        
        # Correction to specified mixing ratio at ground
		rl   = log(b28*pdm[6][2]*sqrt(pdl[2][18]*pdl[2][18])/b01)
		hc01 = pdm[6][6]*pdl[2][12]
		zc01 = pdm[6][5]*pdl[2][11]
		output.d[7] = output.d[7]*ccor(z, rl, hc01, zc01)
        
        # Chemistry correction
		hcc01 = pdm[6][8]*pdl[2][20]
		zcc01 = pdm[6][7]*pdl[2][19]
		rc01  = pdm[6][4]*pdl[2][21]
        
        # Net density corrected at Alt
		output.d[7] = output.d[7]*ccor(z, rc01, hcc01, zcc01)
    end

    # Atomic Nitrogen density
    # Density variation factor at Zlb
    tinf_g7, dfa, plg, ctloc, stloc, c2tloc, s2tloc, c3tloc, s3tloc, apdf, apt = globe7(pd[8], input, flags)
	g14 = flags.sw[22]*tinf_g7
    
    # Diffusive density at Zlb
	db14 = pdm[7][1]*exp(g14)*pd[8][1]
    
    # Diffusive density at Alt
	tz, output.d[8] = densu!(z, db14, tinf, tlb, 14.0, alpha[8], output.t[2], ptm[6], s, mn1, zn1, meso_tn1, meso_tgn1, gsurf=gsurf, re=re)
	dd = output.d[8]
    
    if flags.sw[16] != 0 && (z <= altl[8]) != 0
		# Turbopause
		zh14 = pdm[7][3]
        
        # Mixed density at Zlb
		tz, b14 = densu!(zh14, db14, tinf, tlb, 14.0-xmm, alpha[8]-1.0, output.t[2], ptm[6], s, mn1, zn1, meso_tn1, meso_tgn1, gsurf=gsurf, re=re)
        
        # Mixed density at Alt
		tz, dm14  = densu!(z, b14, tinf, tlb, xmm, 0.0, output.t[2], ptm[6], s, mn1, zn1, meso_tn1, meso_tgn1, gsurf=gsurf, re=re)
		zhm14 = zhm28
        
        # Net density at Alt
        output.d[8] = dnet(output.d[8], dm14, zhm14, xmm, 14.0)
        
        # Correction to specified mixing ratio at ground
        rl   = log(b28*pdm[7][2]*sqrt(pdl[1][3]*pdl[1][3])/b14)
		hc14 = pdm[7][6]*pdl[1][2]
        zc14 = pdm[7][5]*pdl[1][1]
        output.d[8] = output.d[8]*ccor(z, rl, hc14, zc14)
        
        # Chemistry correction
		hcc14 = pdm[7][8]*pdl[1][5]
		zcc14 = pdm[7][7]*pdl[1][4]
		rc14  = pdm[7][4]*pdl[1][6]
        
        # Net density corrected at Alt
        output.d[8] = output.d[8]*ccor(z, rc14, hcc14, zcc14)
        
    end
    
    # Anomalous Oxygen density
    tinf_g7, dfa, plg, ctloc, stloc, c2tloc, s2tloc, c3tloc, s3tloc, apdf, apt = globe7(pd[9], input, flags)
    g16h  = flags.sw[22]*tinf_g7
    db16h = pdm[8][1]*exp(g16h)*pd[9][1]
    tho   = pdm[8][10]*pdl[1][7]
    output.t[2], dd  = densu!(z, db16h, tho, tho, 16.0, alpha[9], output.t[2], ptm[6], s, mn1, zn1, meso_tn1, meso_tgn1, gsurf=gsurf, re=re)
    zsht  = pdm[8][6]
    zmho  = pdm[8][5]
    zsho  = scaleh(zmho, 16.0, tho, gsurf=gsurf, re=re)
    output.d[9] = dd*exp(-zsht/zsho*(exp(-(z-zmho)/zsht)-1.0))

    # Total mass density
    output.d[6] = 1.66e-24 * (  4.0 * output.d[1] 
                             + 16.0 * output.d[2] 
                             + 28.0 * output.d[3]
                             + 32.0 * output.d[4]
                             + 40.0 * output.d[5]
                             +  1.0 * output.d[7]
                             + 14.0 * output.d[8])

    # Temperature
    z = abs(input.alt)
    output.t[2], ddum = densu!(z, 1.0, tinf, tlb, 0.0, 0.0, output.t[2], ptm[6], s, mn1, zn1, meso_tn1, meso_tgn1, gsurf=gsurf, re=re)

    if flags.sw[1] != 0
        for i in 1:9
            output.d[i] = output.d[i]*1.0e6
            # @printf "output.d[%d]: %.10f\n" i output.d[i]
        end
        output.d[6] = output.d[6]/1000
        # @printf "output.d[%d]: %.10f\n" 6 output.d[6]
    end

    return dm28, meso_tn1, meso_tgn1, dfa, plg, ctloc, stloc, c2tloc, s2tloc, s3tloc, c3tloc
end

export density_nrlmsise00
"""
Computes the local atmospheric density using the NRLMSISE00 atmosphere model.

Arguments:
- `epc::Epoch`: Epoch of computation. Used to lookup space weather data
- `x::Array{<:Real, 1}`: Satellite state in geodetic coordinates [lon, lat, alt]
- `use_degrees:Bool`: If `true` interprets geodetic inputs as being in degrees

Returns:
- `rho:Float64`: Local atmospheric density [kg/m^3]

Notes:
1. Uses the gtd7d subroutine of the NRLMSISE00 to compute local atmospheric density. This subroutine includes the contribution of anomalous Oxygen in local density which is important for satellites above 500 km altitude.

References:
1. _Picone, JM, et al._ NRLMSISE-00 empirical model of the atmosphere: Statistical comparisons and scientific issues _Journal of Geophysical Research: Space Physics_
"""
function density_nrlmsise00(epc::Epoch, x::Array{<:Real, 1}; use_degrees::Bool=false)
    # Create NRLMSISE00 model input
    input  = NRLMSISE_Input()
    flags  = NRLMSISE_Flags()
    output = NRLMSISE_Output()

    # Convert geodetic position
    lon = x[1]
    lat = x[2]
    alt = x[3]

    if !use_degrees
        lon *= 180.0/pi
        lat *= 180.0/pi
    end

    # Look-up space weather data
    mjd_ut1 = mjd(epc, tsys="UT1")
    doy     = floor(Int, day_of_year(epc, tsys="UT1"))
    seconds = (mjd_ut1 - floor(mjd_ut1))*86400.0

    ap_array    = zeros(Float64, 7)
    ap_array[1] = ApDailyIndex(mjd_ut1)
    ap_array[2] = ApIndex(epc)
    ap_array[3] = ApIndex(epc - 3*3600) 
    ap_array[4] = ApIndex(epc - 6*3600)
    ap_array[5] = ApIndex(epc - 9*3600)

    # Ap average for 8 3-hour segments before current time
    n = 0
    for i in 12:3:33
        n += 1
        ap_array[6] += ApIndex(epc - i*3600)
    end
    ap_array[6] = ap_array[6]/n

    # Ap average for 8 3-hour segments before current time
    n = 0
    for i in 36:3:57
        n += 1
        ap_array[7] += ApIndex(epc - i*3600)
    end
    ap_array[7] = ap_array[7]/n

    # Set model flags to default
    flags.switches = ones(Int, 24)
    
    # Set input values
    input.year     = 0         # Unused in model 
    input.doy      = doy       # Day of Year
    input.sec      = seconds   # UT
    input.alt      = alt / 1000.0      # [km]
    input.g_lat    = lat       # [deg]
    input.g_lon    = lon       # [deg]
    input.lst      = (input.sec/3600.0 + input.g_lon/15.0) # Set LST to be self consistent
    input.f107A    = f107ObservedAvg(mjd_ut1)  # 81 day average of F10.7 flux, centered on day
    input.f107     = f107Observed(mjd_ut1)     # Daily F10.7 flux for previous day
    input.ap       = ap_array[1] # Daily magnetic index
    input.ap_array = ap_array

    # Run NRLMSISE00 Atmospheric model - including anomalous oxygen
    gtd7d!(input, flags, output)

    # Return local density
    rho = output.d[6]
    return rho
end

end # NRLMSISE00 Module