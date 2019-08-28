########################
# Two Line Element Set #
########################

@enum gravconsttype begin
    wgs72old
    wgs72
    wgs84
end

# Constants of SGP Models 
const twopi  = 2.0 * pi
const xpdotp = 1440.0 / twopi
# const mu     = 398600.5            # in km3 / s2
# const radiusearthkm = 6378.137     # km
# const xke    = 60.0 / sqrt(radiusearthkm^3 / mu)
# const tumin  = 1.0 / xke
# const j2     =  0.00108262998905
# const j3     = -0.00000253215306
# const j4     = -0.00000161098761
# const j3oj2  = j3 / j2

export TLE
"""
Structure for storing a NORAD Two-Line Element 

Atributes:
- `line1::String` First line of Two-Line Element set including checksum
- `line2::String` Second line of Two-Line Element set including checksum
- `epoch::Epoch` Epoch of element set. 
"""
mutable struct TLE
    # String Representation
    line1::String
    line2::String

    # Propagator settings
    whichconst::gravconsttype
    afspc_mode::Bool
    init::Char

    # Generic values
    satnum::Int
    classification::Char
    intldesg::String
    epochyr::Int
    epoch::Epoch
    elnum::Int
    revnum::Int
    method::Char

    # SGP4 Propagator variables #
    # Common 
    a::Float64
    altp::Float64
    alta::Float64
    epochdays::Float64
    jdsatepoch::Float64
    nddot::Float64
    ndot::Float64
    bstar::Float64
    rcse::Float64
    inclo::Float64
    nodeo::Float64
    ecco::Float64
    argpo::Float64
    mo::Float64
    no::Float64
    
    # Near Earth
    isimp::Int
    aycof::Float64
    con41::Float64
    cc1::Float64
    cc4::Float64
    cc5::Float64
    d2::Float64
    d3::Float64
    d4::Float64
    delmo::Float64
    eta::Float64
    argpdot::Float64
    omgcof::Float64
    sinmao::Float64
    t::Float64
    t2cof::Float64
    t3cof::Float64
    t4cof::Float64
    t5cof::Float64
    x1mth2::Float64
    x7thm1::Float64
    mdot::Float64
    nodedot::Float64
    xlcof::Float64
    xmcof::Float64
    nodecf::Float64

    # Deep Space
    irez::Int
    d2201::Float64
    d2211::Float64
    d3210::Float64
    d3222::Float64
    d4410::Float64
    d4422::Float64
    d5220::Float64
    d5232::Float64
    d5421::Float64
    d5433::Float64
    dedt::Float64
    del1::Float64
    del2::Float64
    del3::Float64
    didt::Float64
    dmdt::Float64
    dnodt::Float64
    domdt::Float64
    e3::Float64
    ee2::Float64
    peo::Float64
    pgho::Float64
    pho::Float64
    pinco::Float64
    plo::Float64
    se2::Float64
    se3::Float64
    sgh2::Float64
    sgh3::Float64
    sgh4::Float64
    sh2::Float64
    sh3::Float64
    si2::Float64
    si3::Float64
    sl2::Float64
    sl3::Float64
    sl4::Float64
    gsto::Float64
    xfact::Float64
    xgh2::Float64
    xgh3::Float64
    xgh4::Float64
    xh2::Float64
    xh3::Float64
    xi2::Float64
    xi3::Float64
    xl2::Float64
    xl3::Float64
    xl4::Float64
    xlamo::Float64
    zmol::Float64
    zmos::Float64
    atime::Float64
    xli::Float64
    xni::Float64

    # Error Status Code
    error::UInt8
end

function Base.show(io::IO, tle::TLE)

    s = @sprintf "TLE(epoch: %s, line1: %s,  line2: %s)" tle.epoch tle.line1 tle.line2

    print(io, s)
end

##################
# SGP4 Internals #
##################


export tle_checksum
"""
Compute Two-Line Element checksum for a given line string.

Arguments:
- `line::String` Input line

Returns:
- `sum::Int` 
"""
function tle_checksum(line::String)
    sum = 0

    for c in line
        if c in "0123456789"
            sum += parse(Int, c)
        elseif c == '-'
            sum += 1
        end
    end

    return sum % 10
end

"""
Internal SGP method to convert year and days to month and hours, minutes, 
seconds.
"""
function days2mdhms(year::Int, days::Real)

    lmonth = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    dayofyr = floor(Int, days)
    # ----------------- find month and day of month ---------------- */
    if ( (year % 4) == 0 )
        lmonth[1] = 29
    end

    i = 1
    inttemp = 0
    while ((dayofyr > inttemp + lmonth[i]) && (i < 12))
        inttemp = inttemp + lmonth[i]
        i = i + 1
    end

    mon = i
    day = dayofyr - inttemp

    # ----------------- find hours minutes and seconds ------------- */
    temp = (days - dayofyr) * 24.0
    hr   = floor(Int, temp)
    temp = (temp - hr) * 60.0
    minute  = floor(Int, temp)
    sec  = (temp - minute) * 60.0

    return mon, day, hr, minute, sec
end

"""
Internal SGP method to convert year, month and hours, minutes, seconds to a
Julian Date.
"""
function jday(year::Int, mon::Int, day::Int, hr::Int, minute::Int, sec::Real)
    jd = 367.0 * year -
        floor((7 * (year + floor((mon + 9) / 12.0))) * 0.25) +
        floor( 275 * mon / 9.0 ) +
        day + 1721013.5 +
        ((sec / 60.0 + minute) / 60.0 + hr) / 24.0
    
    return jd
end

"""
Function finds the Greenwich sidereal time consistent with SGP4 propagator
"""
function gstime(jdut1::Real)
    tut1 = (jdut1 - 2451545.0) / 36525.0
    temp = -6.2e-6* tut1 * tut1 * tut1 + 0.093104 * tut1 * tut1 +
            (876600.0*3600 + 8640184.812866) * tut1 + 67310.54841  # sec
    temp = (temp * DEG2RAD / 240.0) % 2*pi # 360/86400 = 1/240, to deg, to rad

    # Check Quadrants
    if (temp < 0.0)
        temp += 2*pi
    end

    return temp
end

"""
/* -----------------------------------------------------------------------------
*
*                           function getgravconst
*
*  this function gets constants for the propagator. note that mu is identified to
*    facilitiate comparisons with newer models. the common useage is wgs72.
*
*  author        : david vallado                  719-573-2600   21 jul 2006
*
*  inputs        :
*    whichconst  - which set of constants to use  wgs72old, wgs72, wgs84
*
*  outputs       :
*    tumin       - minutes in one time unit
*    mu          - earth gravitational parameter
*    radiusearthkm - radius of the earth in km
*    xke         - reciprocal of tumin
*    j2, j3, j4  - un-normalized zonal harmonic values
*    j3oj2       - j3 divided by j2
*
*  locals        :
*
*  coupling      :
*    none
*
*  references    :
*    norad spacetrack report #3
*    vallado, crawford, hujsak, kelso  2006
  --------------------------------------------------------------------------- */
"""
function getgravconst(whichconst::gravconsttype)

    if whichconst == wgs72old::gravconsttype
        # -- wgs-72 low precision str#3 constants --
        mu     = 398600.79964        # in km3 / s2
        radiusearthkm = 6378.135     # km
        xke    = 0.0743669161
        tumin  = 1.0 / xke
        j2     =   0.001082616
        j3     =  -0.00000253881
        j4     =  -0.00000165597
        j3oj2  =  j3 / j2

        return tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2

    elseif whichconst == wgs72::gravconsttype
        # ------------ wgs-72 constants ------------
        mu     = 398600.8            # in km3 / s2
        radiusearthkm = 6378.135     # km
        xke    = 60.0 / sqrt(radiusearthkm*radiusearthkm*radiusearthkm/mu)
        tumin  = 1.0 / xke
        j2     =   0.001082616
        j3     =  -0.00000253881
        j4     =  -0.00000165597
        j3oj2  =  j3 / j2
        
        return tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2

    elseif whichconst == wgs84::gravconsttype
        # ------------ wgs-84 constants ------------
        mu     = 398600.5            # in km3 / s2
        radiusearthkm = 6378.137     # km
        xke    = 60.0 / sqrt(radiusearthkm*radiusearthkm*radiusearthkm/mu)
        tumin  = 1.0 / xke
        j2     =   0.00108262998905
        j3     =  -0.00000253215306
        j4     =  -0.00000161098761
        j3oj2  =  j3 / j2
        
        return tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2

    else
        @error ("Unknown gravity option $whichconst")
    end
end

"""
/*-----------------------------------------------------------------------------
*
*                           procedure initl
*
*  this procedure initializes the spg4 propagator. all the initialization is
*    consolidated here instead of having multiple loops inside other routines.
*
*  author        : david vallado                  719-573-2600   28 jun 2005
*
*  inputs        :
*    ecco        - eccentricity                           0.0 - 1.0
*    epoch       - epoch time in days from jan 0, 1950. 0 hr
*    inclo       - inclination of satellite
*    no          - mean motion of satellite
*    satn        - satellite number
*
*  outputs       :
*    ainv        - 1.0 / a
*    ao          - semi major axis
*    con41       -
*    con42       - 1.0 - 5.0 cos(i)
*    cosio       - cosine of inclination
*    cosio2      - cosio squared
*    eccsq       - eccentricity squared
*    method      - flag for deep space                    'd', 'n'
*    omeosq      - 1.0 - ecco * ecco
*    posq        - semi-parameter squared
*    rp          - radius of perigee
*    rteosq      - square root of (1.0 - ecco*ecco)
*    sinio       - sine of inclination
*    gsto        - gst at time of observation               rad
*    no          - mean motion of satellite
*
*  locals        :
*    ak          -
*    d1          -
*    del         -
*    adel        -
*    po          -
*
*  coupling      :
*    getgravconst
*    gstime      - find greenwich sidereal time from the julian date
*
*  references    :
*    hoots, roehrich, norad spacetrack report #3 1980
*    hoots, norad spacetrack report #6 1986
*    hoots, schumacher and glover 2004
*    vallado, crawford, hujsak, kelso  2006
  ----------------------------------------------------------------------------*/
"""
function initl(satn::Int, whichconst::gravconsttype, ecco::Real, epoch::Real,
            inclo::Real, no::Real, method::Char, afspc_mode::Bool)

    # sgp4fix use old way of finding gst

    #  ----------------------- earth constants ----------------------
    #  sgp4fix identify constants and allow alternate values
    tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2 = getgravconst(whichconst)
    x2o3   = 2.0 / 3.0

    #  ------------- calculate auxillary epoch quantities ----------
    eccsq  = ecco * ecco
    omeosq = 1.0 - eccsq
    rteosq = sqrt(omeosq)
    cosio  = cos(inclo)
    cosio2 = cosio * cosio

    #  ------------------ un-kozai the mean motion -----------------
    ak    = (xke / no) ^ x2o3
    d1    = 0.75 * j2 * (3.0 * cosio2 - 1.0) / (rteosq * omeosq)
    del_  = d1 / (ak * ak)
    adel  = ak * (1.0 - del_ * del_ - del_ *
            (1.0 / 3.0 + 134.0 * del_ * del_ / 81.0))
    del_  = d1/(adel * adel)
    no    = no / (1.0 + del_)

    ao    = (xke / no) ^ x2o3
    sinio = sin(inclo)
    po    = ao * omeosq
    con42 = 1.0 - 5.0 * cosio2
    con41 = -con42-cosio2-cosio2
    ainv  = 1.0 / ao
    posq  = po * po
    rp    = ao * (1.0 - ecco)
    method = 'n'

    #  sgp4fix modern approach to finding sidereal time
    if afspc_mode
        #  sgp4fix use old way of finding gst
        #  count integer number of days from 0 jan 1970
        ts70  = epoch - 7305.0
        ds70 = (ts70 + 1.0e-8) // 1.0
        tfrac = ts70 - ds70

        #  find greenwich location at epoch
        c1    = 1.72027916940703639e-2
        thgr70= 1.7321343856509374
        fk5r  = 5.07551419432269442e-15
        c1p2p = c1 + twopi
        gsto  = (thgr70 + c1*ds70 + c1p2p*tfrac + ts70*ts70*fk5r) % twopi
        
        if gsto < 0.0
            gsto = gsto + twopi
        end

    else
        gsto = gstime(epoch + 2433281.5)
    end

    return (
        no,
        method,
        ainv,  ao,    con41,  con42, cosio,
        cosio2,eccsq, omeosq, posq,
        rp,    rteosq,sinio , gsto,
    )
end

"""
/*-----------------------------------------------------------------------------
*
*                           procedure dscom
*
*  this procedure provides deep space common items used by both the secular
*    and periodics subroutines.  input is provided as shown. this routine
*    used to be called dpper, but the functions inside weren't well organized.
*
*  author        : david vallado                  719-573-2600   28 jun 2005
*
*  inputs        :
*    epoch       -
*    ep          - eccentricity
*    argpp       - argument of perigee
*    tc          -
*    inclp       - inclination
*    nodep       - right ascension of ascending node
*    np          - mean motion
*
*  outputs       :
*    sinim  , cosim  , sinomm , cosomm , snodm  , cnodm
*    day         -
*    e3          -
*    ee2         -
*    em          - eccentricity
*    emsq        - eccentricity squared
*    gam         -
*    peo         -
*    pgho        -
*    pho         -
*    pinco       -
*    plo         -
*    rtemsq      -
*    se2, se3         -
*    sgh2, sgh3, sgh4        -
*    sh2, sh3, si2, si3, sl2, sl3, sl4         -
*    s1, s2, s3, s4, s5, s6, s7          -
*    ss1, ss2, ss3, ss4, ss5, ss6, ss7, sz1, sz2, sz3         -
*    sz11, sz12, sz13, sz21, sz22, sz23, sz31, sz32, sz33        -
*    xgh2, xgh3, xgh4, xh2, xh3, xi2, xi3, xl2, xl3, xl4         -
*    nm          - mean motion
*    z1, z2, z3, z11, z12, z13, z21, z22, z23, z31, z32, z33         -
*    zmol        -
*    zmos        -
*
*  locals        :
*    a1, a2, a3, a4, a5, a6, a7, a8, a9, a10         -
*    betasq      -
*    cc          -
*    ctem, stem        -
*    x1, x2, x3, x4, x5, x6, x7, x8          -
*    xnodce      -
*    xnoi        -
*    zcosg  , zsing  , zcosgl , zsingl , zcosh  , zsinh  , zcoshl , zsinhl ,
*    zcosi  , zsini  , zcosil , zsinil ,
*    zx          -
*    zy          -
*
*  coupling      :
*    none.
*
*  references    :
*    hoots, roehrich, norad spacetrack report #3 1980
*    hoots, norad spacetrack report #6 1986
*    hoots, schumacher and glover 2004
*    vallado, crawford, hujsak, kelso  2006
  ----------------------------------------------------------------------------*/
"""
function dscom(epoch::Real, ep::Real, argpp::Real, tc::Real, inclp::Real,
       nodep::Real, np::Real, e3::Real, ee2::Real, peo::Real, pgho::Real,
       pho::Real, pinco::Real, plo::Real, se2::Real, se3::Real,
       sgh2::Real, sgh3::Real, sgh4::Real, sh2::Real, sh3::Real,
       si2::Real, si3::Real, sl2::Real, sl3::Real, sl4::Real,
       xgh2::Real, xgh3::Real, xgh4::Real, xh2::Real,
       xh3::Real, xi2::Real, xi3::Real, xl2::Real, xl3::Real,
       xl4::Real, zmol::Real, zmos::Real,
     )

    #  ------------------ initialize variables ----------------------
    s1 = 0.0;    s2 = 0.0;    s3 = 0.0;     s4 = 0.0;    s5 = 0.0;
    s6 = 0.0;    s7 = 0.0;    ss1 = 0.0;    ss2 = 0.0;   ss3 = 0.0;
    ss4 = 0.0;   ss5 = 0.0;   ss6 = 0.0;    ss7 = 0.0;   sz1 = 0.0;
    sz2 = 0.0;   sz3 = 0.0;   sz11 = 0.0;   sz12 = 0.0;  sz13 = 0.0;
    sz21 = 0.0;  sz22 = 0.0;  sz23 = 0.0;   sz31 = 0.0;  sz32 = 0.0;
    sz33 = 0.0;  z1 = 0.0;    z2 = 0.0;     z3 = 0.0;
    z11 = 0.0;   z12 = 0.0;   z13 = 0.0;    z21 = 0.0;   z22 = 0.0;
    z23 = 0.0;   z31 = 0.0;   z32 = 0.0;    z33 = 0.0;

    #  -------------------------- constants -------------------------
    zes     =  0.01675
    zel     =  0.05490
    c1ss    =  2.9864797e-6
    c1l     =  4.7968065e-7
    zsinis  =  0.39785416
    zcosis  =  0.91744867
    zcosgs  =  0.1945905
    zsings  = -0.98088458

    #  --------------------- local variables ------------------------
    nm     = np
    em     = ep
    snodm  = sin(nodep)
    cnodm  = cos(nodep)
    sinomm = sin(argpp)
    cosomm = cos(argpp)
    sinim  = sin(inclp)
    cosim  = cos(inclp)
    emsq   = em * em
    betasq = 1.0 - emsq
    rtemsq = sqrt(betasq)

    #  ----------------- initialize lunar solar terms ---------------
    peo    = 0.0
    pinco  = 0.0
    plo    = 0.0
    pgho   = 0.0
    pho    = 0.0
    day    = epoch + 18261.5 + tc / 1440.0
    xnodce = (4.5236020 - 9.2422029e-4 * day) % twopi
    stem   = sin(xnodce)
    ctem   = cos(xnodce)
    zcosil = 0.91375164 - 0.03568096 * ctem
    zsinil = sqrt(1.0 - zcosil * zcosil)
    zsinhl = 0.089683511 * stem / zsinil
    zcoshl = sqrt(1.0 - zsinhl * zsinhl)
    gam    = 5.8351514 + 0.0019443680 * day
    zx     = 0.39785416 * stem / zsinil
    zy     = zcoshl * ctem + 0.91744867 * zsinhl * stem
    zx     = atan(zx, zy)
    zx     = gam + zx - xnodce
    zcosgl = cos(zx)
    zsingl = sin(zx)

    #  ------------------------- do solar terms ---------------------
    zcosg = zcosgs
    zsing = zsings
    zcosi = zcosis
    zsini = zsinis
    zcosh = cnodm
    zsinh = snodm
    cc    = c1ss
    xnoi  = 1.0 / nm

    for lsflg in 1:2

        a1  =   zcosg * zcosh + zsing * zcosi * zsinh
        a3  =  -zsing * zcosh + zcosg * zcosi * zsinh
        a7  =  -zcosg * zsinh + zsing * zcosi * zcosh
        a8  =   zsing * zsini
        a9  =   zsing * zsinh + zcosg * zcosi * zcosh
        a10 =   zcosg * zsini
        a2  =   cosim * a7 + sinim * a8
        a4  =   cosim * a9 + sinim * a10
        a5  =  -sinim * a7 + cosim * a8
        a6  =  -sinim * a9 + cosim * a10

        x1  =  a1 * cosomm + a2 * sinomm
        x2  =  a3 * cosomm + a4 * sinomm
        x3  = -a1 * sinomm + a2 * cosomm
        x4  = -a3 * sinomm + a4 * cosomm
        x5  =  a5 * sinomm
        x6  =  a6 * sinomm
        x7  =  a5 * cosomm
        x8  =  a6 * cosomm

        z31 = 12.0 * x1 * x1 - 3.0 * x3 * x3
        z32 = 24.0 * x1 * x2 - 6.0 * x3 * x4
        z33 = 12.0 * x2 * x2 - 3.0 * x4 * x4
        z1  =  3.0 *  (a1 * a1 + a2 * a2) + z31 * emsq
        z2  =  6.0 *  (a1 * a3 + a2 * a4) + z32 * emsq
        z3  =  3.0 *  (a3 * a3 + a4 * a4) + z33 * emsq
        z11 = -6.0 * a1 * a5 + emsq *  (-24.0 * x1 * x7-6.0 * x3 * x5)
        z12 = -6.0 *  (a1 * a6 + a3 * a5) + emsq *
            (-24.0 * (x2 * x7 + x1 * x8) - 6.0 * (x3 * x6 + x4 * x5))
        z13 = -6.0 * a3 * a6 + emsq * (-24.0 * x2 * x8 - 6.0 * x4 * x6)
        z21 =  6.0 * a2 * a5 + emsq * (24.0 * x1 * x5 - 6.0 * x3 * x7)
        z22 =  6.0 *  (a4 * a5 + a2 * a6) + emsq *
            (24.0 * (x2 * x5 + x1 * x6) - 6.0 * (x4 * x7 + x3 * x8))
        z23 =  6.0 * a4 * a6 + emsq * (24.0 * x2 * x6 - 6.0 * x4 * x8)
        z1  = z1 + z1 + betasq * z31
        z2  = z2 + z2 + betasq * z32
        z3  = z3 + z3 + betasq * z33
        s3  = cc * xnoi
        s2  = -0.5 * s3 / rtemsq
        s4  = s3 * rtemsq
        s1  = -15.0 * em * s4
        s5  = x1 * x3 + x2 * x4
        s6  = x2 * x3 + x1 * x4
        s7  = x2 * x4 - x1 * x3

        #  ----------------------- do lunar terms -------------------
        if lsflg == 1

            ss1   = s1
            ss2   = s2
            ss3   = s3
            ss4   = s4
            ss5   = s5
            ss6   = s6
            ss7   = s7
            sz1   = z1
            sz2   = z2
            sz3   = z3
            sz11  = z11
            sz12  = z12
            sz13  = z13
            sz21  = z21
            sz22  = z22
            sz23  = z23
            sz31  = z31
            sz32  = z32
            sz33  = z33
            zcosg = zcosgl
            zsing = zsingl
            zcosi = zcosil
            zsini = zsinil
            zcosh = zcoshl * cnodm + zsinhl * snodm
            zsinh = snodm * zcoshl - cnodm * zsinhl
            cc    = c1l
        end
    end

    zmol = (4.7199672 + 0.22997150  * day - gam) % twopi
    zmos = (6.2565837 + 0.017201977 * day) % twopi

    #  ------------------------ do solar terms ----------------------
    se2  =   2.0 * ss1 * ss6
    se3  =   2.0 * ss1 * ss7
    si2  =   2.0 * ss2 * sz12
    si3  =   2.0 * ss2 * (sz13 - sz11)
    sl2  =  -2.0 * ss3 * sz2
    sl3  =  -2.0 * ss3 * (sz3 - sz1)
    sl4  =  -2.0 * ss3 * (-21.0 - 9.0 * emsq) * zes
    sgh2 =   2.0 * ss4 * sz32
    sgh3 =   2.0 * ss4 * (sz33 - sz31)
    sgh4 = -18.0 * ss4 * zes
    sh2  =  -2.0 * ss2 * sz22
    sh3  =  -2.0 * ss2 * (sz23 - sz21)

    #  ------------------------ do lunar terms ----------------------
    ee2  =   2.0 * s1 * s6
    e3   =   2.0 * s1 * s7
    xi2  =   2.0 * s2 * z12
    xi3  =   2.0 * s2 * (z13 - z11)
    xl2  =  -2.0 * s3 * z2
    xl3  =  -2.0 * s3 * (z3 - z1)
    xl4  =  -2.0 * s3 * (-21.0 - 9.0 * emsq) * zel
    xgh2 =   2.0 * s4 * z32
    xgh3 =   2.0 * s4 * (z33 - z31)
    xgh4 = -18.0 * s4 * zel
    xh2  =  -2.0 * s2 * z22
    xh3  =  -2.0 * s2 * (z23 - z21)

    return (
        snodm, cnodm, sinim,  cosim, sinomm,
        cosomm,day,   e3,     ee2,   em,
        emsq,  gam,   peo,    pgho,  pho,
        pinco, plo,   rtemsq, se2,   se3,
        sgh2,  sgh3,  sgh4,   sh2,   sh3,
        si2,   si3,   sl2,    sl3,   sl4,
        s1,    s2,    s3,     s4,    s5,
        s6,    s7,    ss1,    ss2,   ss3,
        ss4,   ss5,   ss6,    ss7,   sz1,
        sz2,   sz3,   sz11,   sz12,  sz13,
        sz21,  sz22,  sz23,   sz31,  sz32,
        sz33,  xgh2,  xgh3,   xgh4,  xh2,
        xh3,   xi2,   xi3,    xl2,   xl3,
        xl4,   nm,    z1,     z2,    z3,
        z11,   z12,   z13,    z21,   z22,
        z23,   z31,   z32,    z33,   zmol,
        zmos
        )
    end

"""
/* -----------------------------------------------------------------------------
*
*                           procedure dpper
*
*  this procedure provides deep space long period periodic contributions
*    to the mean elements.  by design, these periodics are zero at epoch.
*    this used to be dscom which included initialization, but it's really a
*    recurring function.
*
*  author        : david vallado                  719-573-2600   28 jun 2005
*
*  inputs        :
*    e3          -
*    ee2         -
*    peo         -
*    pgho        -
*    pho         -
*    pinco       -
*    plo         -
*    se2 , se3 , sgh2, sgh3, sgh4, sh2, sh3, si2, si3, sl2, sl3, sl4 -
*    t           -
*    xh2, xh3, xi2, xi3, xl2, xl3, xl4 -
*    zmol        -
*    zmos        -
*    ep          - eccentricity                           0.0 - 1.0
*    inclo       - inclination - needed for lyddane modification
*    nodep       - right ascension of ascending node
*    argpp       - argument of perigee
*    mp          - mean anomaly
*
*  outputs       :
*    ep          - eccentricity                           0.0 - 1.0
*    inclp       - inclination
*    nodep        - right ascension of ascending node
*    argpp       - argument of perigee
*    mp          - mean anomaly
*
*  locals        :
*    alfdp       -
*    betdp       -
*    cosip  , sinip  , cosop  , sinop  ,
*    dalf        -
*    dbet        -
*    dls         -
*    f2, f3      -
*    pe          -
*    pgh         -
*    ph          -
*    pinc        -
*    pl          -
*    sel   , ses   , sghl  , sghs  , shl   , shs   , sil   , sinzf , sis   ,
*    sll   , sls
*    xls         -
*    xnoh        -
*    zf          -
*    zm          -
*
*  coupling      :
*    none.
*
*  references    :
*    hoots, roehrich, norad spacetrack report #3 1980
*    hoots, norad spacetrack report #6 1986
*    hoots, schumacher and glover 2004
*    vallado, crawford, hujsak, kelso  2006
    ----------------------------------------------------------------------------*/
"""
function dpper(tle::TLE, inclo::Real, init::Char, 
        ep::Real, inclp::Real, nodep::Real, argpp::Real, mp::Real, afspc_mode::Bool)

    # Copy satellite attributes into local variables for convenience
    # and symmetry in writing formulae.

    e3 = tle.e3
    ee2 = tle.ee2
    peo = tle.peo
    pgho = tle.pgho
    pho = tle.pho
    pinco = tle.pinco
    plo = tle.plo
    se2 = tle.se2
    se3 = tle.se3
    sgh2 = tle.sgh2
    sgh3 = tle.sgh3
    sgh4 = tle.sgh4
    sh2 = tle.sh2
    sh3 = tle.sh3
    si2 = tle.si2
    si3 = tle.si3
    sl2 = tle.sl2
    sl3 = tle.sl3
    sl4 = tle.sl4
    t = tle.t
    xgh2 = tle.xgh2
    xgh3 = tle.xgh3
    xgh4 = tle.xgh4
    xh2 = tle.xh2
    xh3 = tle.xh3
    xi2 = tle.xi2
    xi3 = tle.xi3
    xl2 = tle.xl2
    xl3 = tle.xl3
    xl4 = tle.xl4
    zmol = tle.zmol
    zmos = tle.zmos

    #  ---------------------- constants -----------------------------
    zns   = 1.19459e-5
    zes   = 0.01675
    znl   = 1.5835218e-4
    zel   = 0.05490

    #  --------------- calculate time varying periodics -----------
    zm    = zmos + zns * t

    # be sure that the initial call has time set to zero
    if init == 'y'
        zm = zmos
    end

    zf    = zm + 2.0 * zes * sin(zm)
    sinzf = sin(zf)
    f2    =  0.5 * sinzf * sinzf - 0.25
    f3    = -0.5 * sinzf * cos(zf)
    ses   = se2* f2 + se3 * f3
    sis   = si2 * f2 + si3 * f3
    sls   = sl2 * f2 + sl3 * f3 + sl4 * sinzf
    sghs  = sgh2 * f2 + sgh3 * f3 + sgh4 * sinzf
    shs   = sh2 * f2 + sh3 * f3
    zm    = zmol + znl * t

    if init == 'y'
        zm = zmol
    end

    zf    = zm + 2.0 * zel * sin(zm)
    sinzf = sin(zf)
    f2    =  0.5 * sinzf * sinzf - 0.25
    f3    = -0.5 * sinzf * cos(zf)
    sel   = ee2 * f2 + e3 * f3
    sil   = xi2 * f2 + xi3 * f3
    sll   = xl2 * f2 + xl3 * f3 + xl4 * sinzf
    sghl  = xgh2 * f2 + xgh3 * f3 + xgh4 * sinzf
    shll  = xh2 * f2 + xh3 * f3
    pe    = ses + sel
    pinc  = sis + sil
    pl    = sls + sll
    pgh   = sghs + sghl
    ph    = shs + shll

    if init == 'n'

        pe    = pe - peo;
        pinc  = pinc - pinco;
        pl    = pl - plo;
        pgh   = pgh - pgho;
        ph    = ph - pho;
        inclp = inclp + pinc;
        ep    = ep + pe;
        sinip = sin(inclp);
        cosip = cos(inclp);

        """
        /* ----------------- apply periodics directly ------------ */
        //  sgp4fix for lyddane choice
        //  strn3 used original inclination - this is technically feasible
        //  gsfc used perturbed inclination - also technically feasible
        //  probably best to readjust the 0.2 limit value and limit discontinuity
        //  0.2 rad = 11.45916 deg
        //  use next line for original strn3 approach and original inclination
        //  if (inclo >= 0.2)
        //  use next line for gsfc version and perturbed inclination
        """

        if inclp >= 0.2

            ph /= sinip
            pgh -= cosip * ph
            argpp += pgh
            nodep += ph
            mp += pl

        else

            #  ---- apply periodics with lyddane modification ----
            sinop  = sin(nodep);
            cosop  = cos(nodep);
            alfdp  = sinip * sinop;
            betdp  = sinip * cosop;
            dalf   =  ph * cosop + pinc * cosip * sinop;
            dbet   = -ph * sinop + pinc * cosip * cosop;
            alfdp  = alfdp + dalf;
            betdp  = betdp + dbet;
            nodep = nodep >= 0.0 ? nodep % twopi : -(-nodep % twopi)

            #   sgp4fix for afspc written intrinsic functions
            #  nodep used without a trigonometric function ahead
            if nodep < 0.0 && afspc_mode
                nodep = nodep + twopi
            end
            xls = mp + argpp + pl + pgh + (cosip - pinc * sinip) * nodep
            xnoh   = nodep;
            nodep  = atan(alfdp, betdp);
            #   sgp4fix for afspc written intrinsic functions
            #  nodep used without a trigonometric function ahead
            if nodep < 0.0 && afspc_mode
                nodep = nodep + twopi;
            end
            if abs(xnoh - nodep) > pi
                if nodep < xnoh
                    nodep = nodep + twopi;
                else
                    nodep = nodep - twopi;
                end
            mp += pl
            argpp = xls - mp - cosip * nodep;
            end
        end
    end

    return ep, inclp, nodep, argpp, mp
end

"""
/*-----------------------------------------------------------------------------
*
*                           procedure dsinit
*
*  this procedure provides deep space contributions to mean motion dot due
*    to geopotential resonance with half day and one day orbits.
*
*  author        : david vallado                  719-573-2600   28 jun 2005
*
*  inputs        :
*    cosim, sinim-
*    emsq        - eccentricity squared
*    argpo       - argument of perigee
*    s1, s2, s3, s4, s5      -
*    ss1, ss2, ss3, ss4, ss5 -
*    sz1, sz3, sz11, sz13, sz21, sz23, sz31, sz33 -
*    t           - time
*    tc          -
*    gsto        - greenwich sidereal time                   rad
*    mo          - mean anomaly
*    mdot        - mean anomaly dot (rate)
*    no          - mean motion
*    nodeo       - right ascension of ascending node
*    nodedot     - right ascension of ascending node dot (rate)
*    xpidot      -
*    z1, z3, z11, z13, z21, z23, z31, z33 -
*    eccm        - eccentricity
*    argpm       - argument of perigee
*    inclm       - inclination
*    mm          - mean anomaly
*    xn          - mean motion
*    nodem       - right ascension of ascending node
*
*  outputs       :
*    em          - eccentricity
*    argpm       - argument of perigee
*    inclm       - inclination
*    mm          - mean anomaly
*    nm          - mean motion
*    nodem       - right ascension of ascending node
*    irez        - flag for resonance           0-none, 1-one day, 2-half day
*    atime       -
*    d2201, d2211, d3210, d3222, d4410, d4422, d5220, d5232, d5421, d5433    -
*    dedt        -
*    didt        -
*    dmdt        -
*    dndt        -
*    dnodt       -
*    domdt       -
*    del1, del2, del3        -
*    ses  , sghl , sghs , sgs  , shl  , shs  , sis  , sls
*    theta       -
*    xfact       -
*    xlamo       -
*    xli         -
*    xni
*
*  locals        :
*    ainv2       -
*    aonv        -
*    cosisq      -
*    eoc         -
*    f220, f221, f311, f321, f322, f330, f441, f442, f522, f523, f542, f543  -
*    g200, g201, g211, g300, g310, g322, g410, g422, g520, g521, g532, g533  -
*    sini2       -
*    temp        -
*    temp1       -
*    theta       -
*    xno2        -
*
*  coupling      :
*    getgravconst
*
*  references    :
*    hoots, roehrich, norad spacetrack report #3 1980
*    hoots, norad spacetrack report #6 1986
*    hoots, schumacher and glover 2004
*    vallado, crawford, hujsak, kelso  2006
  ----------------------------------------------------------------------------*/
"""
function dsinit(whichconst::gravconsttype,
       cosim::Real,  emsq::Real,   argpo::Real,   s1::Real,     s2::Real,
       s3::Real,     s4::Real,     s5::Real,      sinim::Real,  ss1::Real,
       ss2::Real,    ss3::Real,    ss4::Real,     ss5::Real,    sz1::Real,
       sz3::Real,    sz11::Real,   sz13::Real,    sz21::Real,   sz23::Real,
       sz31::Real,   sz33::Real,   t::Real,       tc::Real,     gsto::Real,
       mo::Real,     mdot::Real,   no::Real,      nodeo::Real,  nodedot::Real,
       xpidot::Real, z1::Real,     z3::Real,      z11::Real,    z13::Real,
       z21::Real,    z23::Real,    z31::Real,     z33::Real,    ecco::Real,
       eccsq::Real,  em::Real,    argpm::Real,  inclm::Real, mm::Real,
       nm::Real,    nodem::Real,
       irez::Int,
       atime::Real, d2201::Real, d2211::Real,  d3210::Real, d3222::Real,
       d4410::Real, d4422::Real, d5220::Real,  d5232::Real, d5421::Real,
       d5433::Real, dedt::Real,  didt::Real,   dmdt::Real,
       dnodt::Real, domdt::Real, del1::Real,   del2::Real,  del3::Real,
       xfact::Real, xlamo::Real, xli::Real,    xni::Real,)

    q22    = 1.7891679e-6
    q31    = 2.1460748e-6
    q33    = 2.2123015e-7
    root22 = 1.7891679e-6
    root44 = 7.3636953e-9
    root54 = 2.1765803e-9
    rptim  = 4.37526908801129966e-3 # equates to 7.29211514668855e-5 rad/sec
    root32 = 3.7393792e-7
    root52 = 1.1428639e-7
    x2o3   = 2.0 / 3.0
    znl    = 1.5835218e-4
    zns    = 1.19459e-5

    #  sgp4fix identify constants and allow alternate values
    tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2 = getgravconst(whichconst)
    # xke = whichconst.xke

    #  -------------------- deep space initialization ------------
    irez = 0
    if 0.0034906585 < nm < 0.0052359877
        irez = 1
    end

    if (8.26e-3 <= nm <= 9.24e-3) && em >= 0.5
        irez = 2
    end

    #  ------------------------ do solar terms -------------------
    ses  =  ss1 * zns * ss5
    sis  =  ss2 * zns * (sz11 + sz13)
    sls  = -zns * ss3 * (sz1 + sz3 - 14.0 - 6.0 * emsq)
    sghs =  ss4 * zns * (sz31 + sz33 - 6.0)
    shs  = -zns * ss2 * (sz21 + sz23)
    #  sgp4fix for 180 deg incl
    if inclm < 5.2359877e-2 || inclm > pi - 5.2359877e-2
        shs = 0.0
    end

    if sinim != 0.0
        shs = shs / sinim
        sgs  = sghs - cosim * shs
    end

    #  ------------------------- do lunar terms ------------------
    dedt = ses + s1 * znl * s5
    didt = sis + s2 * znl * (z11 + z13)
    dmdt = sls - znl * s3 * (z1 + z3 - 14.0 - 6.0 * emsq)
    sghl = s4 * znl * (z31 + z33 - 6.0)
    shll = -znl * s2 * (z21 + z23)

    #  sgp4fix for 180 deg incl
    if (inclm < 5.2359877e-2) || (inclm > pi - 5.2359877e-2)
        shll = 0.0
    end

    domdt = sgs + sghl
    dnodt = shs
    if sinim != 0.0
        domdt = domdt - cosim / sinim * shll
        dnodt = dnodt + shll / sinim
    end


    #  ----------- calculate deep space resonance effects --------
    dndt   = 0.0
    theta  = (gsto + tc * rptim) % twopi
    em     = em + dedt * t
    inclm  = inclm + didt * t
    argpm  = argpm + domdt * t
    nodem  = nodem + dnodt * t
    mm     = mm + dmdt * t
    """
    //   sgp4fix for negative inclinations
    //   the following if statement should be commented out
    //if (inclm < 0.0)
    //  {
    //    inclm  = -inclm;
    //    argpm  = argpm - pi;
    //    nodem = nodem + pi;
    //  }
    """

    #  -------------- initialize the resonance terms -------------
    if irez != 0

        aonv = (nm / xke)  ^ x2o3

        #  ---------- geopotential resonance for 12 hour orbits ------
        if irez == 2

            cosisq = cosim * cosim
            emo    = em
            em     = ecco
            emsqo  = emsq
            emsq   = eccsq
            eoc    = em * emsq
            g201   = -0.306 - (em - 0.64) * 0.440

            if em <= 0.65

                g211 =    3.616  -  13.2470 * em +  16.2900 * emsq
                g310 =  -19.302  + 117.3900 * em - 228.4190 * emsq +  156.5910 * eoc
                g322 =  -18.9068 + 109.7927 * em - 214.6334 * emsq +  146.5816 * eoc
                g410 =  -41.122  + 242.6940 * em - 471.0940 * emsq +  313.9530 * eoc
                g422 = -146.407  + 841.8800 * em - 1629.014 * emsq + 1083.4350 * eoc
                g520 = -532.114  + 3017.977 * em - 5740.032 * emsq + 3708.2760 * eoc

            else

                g211 =   -72.099 +   331.819 * em -   508.738 * emsq +   266.724 * eoc
                g310 =  -346.844 +  1582.851 * em -  2415.925 * emsq +  1246.113 * eoc
                g322 =  -342.585 +  1554.908 * em -  2366.899 * emsq +  1215.972 * eoc
                g410 = -1052.797 +  4758.686 * em -  7193.992 * emsq +  3651.957 * eoc
                g422 = -3581.690 + 16178.110 * em - 24462.770 * emsq + 12422.520 * eoc
                
                if em > 0.715
                    g520 =-5149.66 + 29936.92 * em - 54087.36 * emsq + 31324.56 * eoc
                else
                    g520 = 1464.74 -  4664.75 * em +  3763.64 * emsq
                end
            end

            if em < 0.7
                g533 = -919.22770 + 4988.6100 * em - 9064.7700 * emsq + 5542.21  * eoc
                g521 = -822.71072 + 4568.6173 * em - 8491.4146 * emsq + 5337.524 * eoc
                g532 = -853.66600 + 4690.2500 * em - 8624.7700 * emsq + 5341.4  * eoc

            else
                g533 =-37995.780 + 161616.52 * em - 229838.20 * emsq + 109377.94 * eoc
                g521 =-51752.104 + 218913.95 * em - 309468.16 * emsq + 146349.42 * eoc
                g532 =-40023.880 + 170470.89 * em - 242699.48 * emsq + 115605.82 * eoc
            end

            sini2=  sinim * sinim
            f220 =  0.75 * (1.0 + 2.0 * cosim+cosisq)
            f221 =  1.5 * sini2
            f321 =  1.875 * sinim  *  (1.0 - 2.0 * cosim - 3.0 * cosisq)
            f322 = -1.875 * sinim  *  (1.0 + 2.0 * cosim - 3.0 * cosisq)
            f441 = 35.0 * sini2 * f220
            f442 = 39.3750 * sini2 * sini2
            f522 =  9.84375 * sinim * (sini2 * (1.0 - 2.0 * cosim- 5.0 * cosisq) +
                    0.33333333 * (-2.0 + 4.0 * cosim + 6.0 * cosisq) )
            f523 = sinim * (4.92187512 * sini2 * (-2.0 - 4.0 * cosim +
                10.0 * cosisq) + 6.56250012 * (1.0+2.0 * cosim - 3.0 * cosisq))
            f542 = 29.53125 * sinim * (2.0 - 8.0 * cosim+cosisq *
                (-12.0 + 8.0 * cosim + 10.0 * cosisq))
            f543 = 29.53125 * sinim * (-2.0 - 8.0 * cosim+cosisq *
                (12.0 + 8.0 * cosim - 10.0 * cosisq))
            xno2  =  nm * nm
            ainv2 =  aonv * aonv
            temp1 =  3.0 * xno2 * ainv2
            temp  =  temp1 * root22
            d2201 =  temp * f220 * g201
            d2211 =  temp * f221 * g211
            temp1 =  temp1 * aonv
            temp  =  temp1 * root32
            d3210 =  temp * f321 * g310
            d3222 =  temp * f322 * g322
            temp1 =  temp1 * aonv
            temp  =  2.0 * temp1 * root44
            d4410 =  temp * f441 * g410
            d4422 =  temp * f442 * g422
            temp1 =  temp1 * aonv
            temp  =  temp1 * root52
            d5220 =  temp * f522 * g520
            d5232 =  temp * f523 * g532
            temp  =  2.0 * temp1 * root54
            d5421 =  temp * f542 * g521
            d5433 =  temp * f543 * g533
            xlamo =  (mo + nodeo + nodeo-theta - theta) % twopi
            xfact =  mdot + dmdt + 2.0 * (nodedot + dnodt - rptim) - no
            em    = emo
            emsq  = emsqo
        end

        #  ---------------- synchronous resonance terms --------------
        if irez == 1

            g200  = 1.0 + emsq * (-2.5 + 0.8125 * emsq)
            g310  = 1.0 + 2.0 * emsq
            g300  = 1.0 + emsq * (-6.0 + 6.60937 * emsq)
            f220  = 0.75 * (1.0 + cosim) * (1.0 + cosim)
            f311  = 0.9375 * sinim * sinim * (1.0 + 3.0 * cosim) - 0.75 * (1.0 + cosim)
            f330  = 1.0 + cosim
            f330  = 1.875 * f330 * f330 * f330
            del1  = 3.0 * nm * nm * aonv * aonv
            del2  = 2.0 * del1 * f220 * g200 * q22
            del3  = 3.0 * del1 * f330 * g300 * q33 * aonv
            del1  = del1 * f311 * g310 * q31 * aonv
            xlamo = (mo + nodeo + argpo - theta) % twopi
            xfact = mdot + xpidot - rptim + dmdt + domdt + dnodt - no
        end

        #  ------------ for sgp4, initialize the integrator ----------
        xli   = xlamo
        xni   = no
        atime = 0.0
        nm    = no + dndt
    end

    return (
        em,    argpm,  inclm, mm,
        nm,    nodem,
        irez, atime,
        d2201, d2211,  d3210, d3222,
        d4410, d4422, d5220,  d5232,
        d5421, d5433, dedt,  didt,
        dmdt,  dndt, dnodt, domdt,
        del1,   del2,  del3, xfact,
        xlamo, xli,    xni,
    )
end

"""
/*-----------------------------------------------------------------------------
*
*                           procedure dspace
*
*  this procedure provides deep space contributions to mean elements for
*    perturbing third body.  these effects have been averaged over one
*    revolution of the sun and moon.  for earth resonance effects, the
*    effects have been averaged over no revolutions of the satellite.
*    (mean motion)
*
*  author        : david vallado                  719-573-2600   28 jun 2005
*
*  inputs        :
*    d2201, d2211, d3210, d3222, d4410, d4422, d5220, d5232, d5421, d5433 -
*    dedt        -
*    del1, del2, del3  -
*    didt        -
*    dmdt        -
*    dnodt       -
*    domdt       -
*    irez        - flag for resonance           0-none, 1-one day, 2-half day
*    argpo       - argument of perigee
*    argpdot     - argument of perigee dot (rate)
*    t           - time
*    tc          -
*    gsto        - gst
*    xfact       -
*    xlamo       -
*    no          - mean motion
*    atime       -
*    em          - eccentricity
*    ft          -
*    argpm       - argument of perigee
*    inclm       - inclination
*    xli         -
*    mm          - mean anomaly
*    xni         - mean motion
*    nodem       - right ascension of ascending node
*
*  outputs       :
*    atime       -
*    em          - eccentricity
*    argpm       - argument of perigee
*    inclm       - inclination
*    xli         -
*    mm          - mean anomaly
*    xni         -
*    nodem       - right ascension of ascending node
*    dndt        -
*    nm          - mean motion
*
*  locals        :
*    delt        -
*    ft          -
*    theta       -
*    x2li        -
*    x2omi       -
*    xl          -
*    xldot       -
*    xnddt       -
*    xndt        -
*    xomi        -
*
*  coupling      :
*    none        -
*
*  references    :
*    hoots, roehrich, norad spacetrack report #3 1980
*    hoots, norad spacetrack report #6 1986
*    hoots, schumacher and glover 2004
*    vallado, crawford, hujsak, kelso  2006
  ----------------------------------------------------------------------------*/
"""
function dspace(irez::Int,
       d2201::Real,  d2211::Real,  d3210::Real,   d3222::Real,  d4410::Real,
       d4422::Real,  d5220::Real,  d5232::Real,   d5421::Real,  d5433::Real,
       dedt::Real,   del1::Real,   del2::Real,    del3::Real,   didt::Real,
       dmdt::Real,   dnodt::Real,  domdt::Real,   argpo::Real,  argpdot::Real,
       t::Real,      tc::Real,     gsto::Real,    xfact::Real,  xlamo::Real,
       no::Real,
       atime::Real, em::Real,    argpm::Real,  inclm::Real, xli::Real,
       mm::Real,    xni::Real,   nodem::Real,  nm::Real,)

    fasx2 = 0.13130908
    fasx4 = 2.8843198
    fasx6 = 0.37448087
    g22   = 5.7686396
    g32   = 0.95240898
    g44   = 1.8014998
    g52   = 1.0508330
    g54   = 4.4108898
    rptim = 4.37526908801129966e-3 # equates to 7.29211514668855e-5 rad/sec
    stepp =    720.0
    stepn =   -720.0
    step2 = 259200.0

    #  ----------- calculate deep space resonance effects -----------
    dndt   = 0.0
    theta  = (gsto + tc * rptim) % twopi
    em     = em + dedt * t

    inclm  = inclm + didt * t
    argpm  = argpm + domdt * t
    nodem  = nodem + dnodt * t
    mm     = mm + dmdt * t

    """
    //   sgp4fix for negative inclinations
    //   the following if statement should be commented out
    //  if (inclm < 0.0)
    // {
    //    inclm = -inclm
    //    argpm = argpm - pi
    //    nodem = nodem + pi
    //  }

    /* - update resonances : numerical (euler-maclaurin) integration - */
    /* ------------------------- epoch restart ----------------------  */
    //   sgp4fix for propagator problems
    //   the following integration works for negative time steps and periods
    //   the specific changes are unknown because the original code was so convoluted

    // sgp4fix take out atime = 0.0 and fix for faster operation
    """

    ft    = 0.0
    xndt  = 0.0
    xnddt = 0.0
    xldot = 0.0
    if irez != 0

        #  sgp4fix streamline check
        if atime == 0.0 || t * atime <= 0.0 || abs(t) < abs(atime)

            atime  = 0.0
            xni    = no
            xli    = xlamo
        end

        # sgp4fix move check outside loop
        if t > 0.0
            delt = stepp
        else
            delt = stepn
        end

        iretn = 381 # added for do loop
        iret  =   0 # added for loop
        while iretn == 381

            #  ------------------- dot terms calculated -------------
            #  ----------- near - synchronous resonance terms -------
            if irez != 2

                xndt  = del1 * sin(xli - fasx2) + del2 * sin(2.0 * (xli - fasx4)) +
                        del3 * sin(3.0 * (xli - fasx6))
                xldot = xni + xfact
                xnddt = del1 * cos(xli - fasx2) +
                        2.0 * del2 * cos(2.0 * (xli - fasx4)) +
                        3.0 * del3 * cos(3.0 * (xli - fasx6))
                xnddt = xnddt * xldot

            else

                # --------- near - half-day resonance terms --------
                xomi  = argpo + argpdot * atime
                x2omi = xomi + xomi
                x2li  = xli + xli
                xndt  = (d2201 * sin(x2omi + xli - g22) + d2211 * sin(xli - g22) +
                    d3210 * sin(xomi + xli - g32)  + d3222 * sin(-xomi + xli - g32)+
                    d4410 * sin(x2omi + x2li - g44)+ d4422 * sin(x2li - g44) +
                    d5220 * sin(xomi + xli - g52)  + d5232 * sin(-xomi + xli - g52)+
                    d5421 * sin(xomi + x2li - g54) + d5433 * sin(-xomi + x2li - g54))
                xldot = xni + xfact
                xnddt = (d2201 * cos(x2omi + xli - g22) + d2211 * cos(xli - g22) +
                    d3210 * cos(xomi + xli - g32) + d3222 * cos(-xomi + xli - g32) +
                    d5220 * cos(xomi + xli - g52) + d5232 * cos(-xomi + xli - g52) +
                    2.0 * (d4410 * cos(x2omi + x2li - g44) +
                    d4422 * cos(x2li - g44) + d5421 * cos(xomi + x2li - g54) +
                    d5433 * cos(-xomi + x2li - g54)))
                xnddt = xnddt * xldot
            end

            #  ----------------------- integrator -------------------
            #  sgp4fix move end checks to end of routine
            if abs(t - atime) >= stepp
                iret  = 0
                iretn = 381

            else
                ft    = t - atime
                iretn = 0
            end

            if iretn == 381

                xli   = xli + xldot * delt + xndt * step2
                xni   = xni + xndt * delt + xnddt * step2
                atime = atime + delt
            end
        end

        nm = xni + xndt * ft + xnddt * ft * ft * 0.5
        xl = xli + xldot * ft + xndt * ft * ft * 0.5
        if irez != 1
            mm   = xl - 2.0 * nodem + 2.0 * theta
            dndt = nm - no

        else
            mm   = xl - nodem - argpm + theta
            dndt = nm - no
        end

        nm = no + dndt
    end

    return (
        atime, em,    argpm,  inclm, xli,
        mm,    xni,   nodem,  dndt,  nm,
    )
end

"""
/*-----------------------------------------------------------------------------
*
*                             procedure sgp4init
*
*  this procedure initializes variables for sgp4.
*
*  author        : david vallado                  719-573-2600   28 jun 2005
*
*  inputs        :
*    afspc_mode  - use afspc or improved mode of operation
*    whichconst  - which set of constants to use  72, 84
*    satn        - satellite number
*    bstar       - sgp4 type drag coefficient              kg/m2er
*    ecco        - eccentricity
*    epoch       - epoch time in days from jan 0, 1950. 0 hr
*    argpo       - argument of perigee (output if ds)
*    inclo       - inclination
*    mo          - mean anomaly (output if ds)
*    no          - mean motion
*    nodeo       - right ascension of ascending node
*
*  outputs       :
*    tle      - common values for subsequent calls
*    return code - non-zero on error.
*                   1 - mean elements, ecc >= 1.0 or ecc < -0.001 or a < 0.95 er
*                   2 - mean motion less than 0.0
*                   3 - pert elements, ecc < 0.0  or  ecc > 1.0
*                   4 - semi-latus rectum < 0.0
*                   5 - epoch elements are sub-orbital
*                   6 - satellite has decayed
*
*  locals        :
*    cnodm  , snodm  , cosim  , sinim  , cosomm , sinomm
*    cc1sq  , cc2    , cc3
*    coef   , coef1
*    cosio4      -
*    day         -
*    dndt        -
*    em          - eccentricity
*    emsq        - eccentricity squared
*    eeta        -
*    etasq       -
*    gam         -
*    argpm       - argument of perigee
*    nodem       -
*    inclm       - inclination
*    mm          - mean anomaly
*    nm          - mean motion
*    perige      - perigee
*    pinvsq      -
*    psisq       -
*    qzms24      -
*    rtemsq      -
*    s1, s2, s3, s4, s5, s6, s7          -
*    sfour       -
*    ss1, ss2, ss3, ss4, ss5, ss6, ss7         -
*    sz1, sz2, sz3
*    sz11, sz12, sz13, sz21, sz22, sz23, sz31, sz32, sz33        -
*    tc          -
*    temp        -
*    temp1, temp2, temp3       -
*    tsi         -
*    xpidot      -
*    xhdot1      -
*    z1, z2, z3          -
*    z11, z12, z13, z21, z22, z23, z31, z32, z33         -
*
*  coupling      :
*    getgravconst-
*    initl       -
*    dscom       -
*    dpper       -
*    dsinit      -
*    sgp4        -
*
*  references    :
*    hoots, roehrich, norad spacetrack report #3 1980
*    hoots, norad spacetrack report #6 1986
*    hoots, schumacher and glover 2004
*    vallado, crawford, hujsak, kelso  2006
  ----------------------------------------------------------------------------*/
"""
function sgp4init(whichconst::gravconsttype, afspc_mode::Bool, satn::Int, epoch::Float64,
       xbstar::Real, xecco::Real, xargpo::Real,
       xinclo::Real, xmo::Real, xno::Real, xnodeo::Real,  tle::TLE)

    """
    /* ------------------------ initialization --------------------- */
    // sgp4fix divisor for divide by zero check on inclination
    // the old check used 1.0 + cos(pi-1.0e-9), but then compared it to
    // 1.5 e-12, so the threshold was changed to 1.5e-12 for consistency
    """
    temp4    =   1.5e-12;

    #  ----------- set all near earth variables to zero ------------
    tle.isimp   = 0;   tle.method = 'n'; tle.aycof    = 0.0;
    tle.con41   = 0.0; tle.cc1    = 0.0; tle.cc4      = 0.0;
    tle.cc5     = 0.0; tle.d2     = 0.0; tle.d3       = 0.0;
    tle.d4      = 0.0; tle.delmo  = 0.0; tle.eta      = 0.0;
    tle.argpdot = 0.0; tle.omgcof = 0.0; tle.sinmao   = 0.0;
    tle.t       = 0.0; tle.t2cof  = 0.0; tle.t3cof    = 0.0;
    tle.t4cof   = 0.0; tle.t5cof  = 0.0; tle.x1mth2   = 0.0;
    tle.x7thm1  = 0.0; tle.mdot   = 0.0; tle.nodedot  = 0.0;
    tle.xlcof   = 0.0; tle.xmcof  = 0.0; tle.nodecf   = 0.0;

    #  ----------- set all deep space variables to zero ------------
    tle.irez  = 0;   tle.d2201 = 0.0; tle.d2211 = 0.0;
    tle.d3210 = 0.0; tle.d3222 = 0.0; tle.d4410 = 0.0;
    tle.d4422 = 0.0; tle.d5220 = 0.0; tle.d5232 = 0.0;
    tle.d5421 = 0.0; tle.d5433 = 0.0; tle.dedt  = 0.0;
    tle.del1  = 0.0; tle.del2  = 0.0; tle.del3  = 0.0;
    tle.didt  = 0.0; tle.dmdt  = 0.0; tle.dnodt = 0.0;
    tle.domdt = 0.0; tle.e3    = 0.0; tle.ee2   = 0.0;
    tle.peo   = 0.0; tle.pgho  = 0.0; tle.pho   = 0.0;
    tle.pinco = 0.0; tle.plo   = 0.0; tle.se2   = 0.0;
    tle.se3   = 0.0; tle.sgh2  = 0.0; tle.sgh3  = 0.0;
    tle.sgh4  = 0.0; tle.sh2   = 0.0; tle.sh3   = 0.0;
    tle.si2   = 0.0; tle.si3   = 0.0; tle.sl2   = 0.0;
    tle.sl3   = 0.0; tle.sl4   = 0.0; tle.gsto  = 0.0;
    tle.xfact = 0.0; tle.xgh2  = 0.0; tle.xgh3  = 0.0;
    tle.xgh4  = 0.0; tle.xh2   = 0.0; tle.xh3   = 0.0;
    tle.xi2   = 0.0; tle.xi3   = 0.0; tle.xl2   = 0.0;
    tle.xl3   = 0.0; tle.xl4   = 0.0; tle.xlamo = 0.0;
    tle.zmol  = 0.0; tle.zmos  = 0.0; tle.atime = 0.0;
    tle.xli   = 0.0; tle.xni   = 0.0;

    """
    // sgp4fix - note the following variables are also passed directly via tle.
    // it is possible to streamline the sgp4init call by deleting the "x"
    // variables, but the user would need to set the tle.* values first. we
    // include the additional assignments in case twoline2rv is not used.
    """
    tle.bstar   = xbstar
    tle.ecco    = xecco
    tle.argpo   = xargpo
    tle.inclo   = xinclo
    tle.mo	    = xmo
    tle.no	    = xno
    tle.nodeo   = xnodeo

    #  sgp4fix add opsmode
    tle.afspc_mode = afspc_mode

    #  ------------------------ earth constants -----------------------
    #  sgp4fix identify constants and allow alternate values
    tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2 = getgravconst(whichconst)
    ss     = 78.0 / radiusearthkm + 1.0

    #  sgp4fix use multiply for speed instead of pow
    qzms2ttemp = (120.0 - 78.0) / radiusearthkm
    qzms2t = qzms2ttemp * qzms2ttemp * qzms2ttemp * qzms2ttemp
    x2o3   =  2.0 / 3.0

    tle.init = 'y'
    tle.t	 = 0.0

    (
        tle.no, method, ainv,  ao,    tle.con41,  con42, cosio, cosio2,eccsq, 
        omeosq, posq, rp, rteosq,sinio, tle.gsto,
    ) = initl(
        satn, whichconst, tle.ecco, epoch, tle.inclo, tle.no, tle.method,
        tle.afspc_mode
    )
    tle.error = 0

    """
    // sgp4fix remove this check as it is unnecessary
    // the mrt check in sgp4 handles decaying satellite cases even if the starting
    // condition is below the surface of te earth
    //     if (rp < 1.0)
    //       {
    //         printf("# *** satn%d epoch elts sub-orbital ***\n", satn)
    //         tle.error = 5
    //       }
    """

    if (omeosq >= 0.0 || tle.no >= 0.0)

        tle.isimp = 0
        if (rp < 220.0 / radiusearthkm + 1.0)
            tle.isimp = 1
        end

        sfour  = ss
        qzms24 = qzms2t
        perige = (rp - 1.0) * radiusearthkm

        #  - for perigees below 156 km, s and qoms2t are altered -
        if perige < 156.0

            sfour = perige - 78.0
            if perige < 98.0
                sfour = 20.0
            end

            #  sgp4fix use multiply for speed instead of pow
            qzms24temp =  (120.0 - sfour) / radiusearthkm
            qzms24 = qzms24temp * qzms24temp * qzms24temp * qzms24temp
            sfour  = sfour / radiusearthkm + 1.0

        end

        pinvsq = 1.0 / posq

        tsi  = 1.0 / (ao - sfour)
        tle.eta  = ao * tle.ecco * tsi
        etasq = tle.eta * tle.eta
        eeta  = tle.ecco * tle.eta
        psisq = abs(1.0 - etasq)
        coef  = qzms24 * (tsi ^ 4.0)
        coef1 = coef / (psisq ^ 3.5)
        cc2   = coef1 * tle.no * (ao * (1.0 + 1.5 * etasq + eeta *
                    (4.0 + etasq)) + 0.375 * j2 * tsi / psisq * tle.con41 *
                    (8.0 + 3.0 * etasq * (8.0 + etasq)))
        tle.cc1   = tle.bstar * cc2
        cc3   = 0.0
        if (tle.ecco > 1.0e-4)
            cc3 = -2.0 * coef * tsi * j3oj2 * tle.no * sinio / tle.ecco
        end

        tle.x1mth2 = 1.0 - cosio2
        tle.cc4    = 2.0 * tle.no * coef1 * ao * omeosq *
                        (tle.eta * (2.0 + 0.5 * etasq) + tle.ecco *
                        (0.5 + 2.0 * etasq) - j2 * tsi / (ao * psisq) *
                        (-3.0 * tle.con41 * (1.0 - 2.0 * eeta + etasq *
                        (1.5 - 0.5 * eeta)) + 0.75 * tle.x1mth2 *
                        (2.0 * etasq - eeta * (1.0 + etasq)) * cos(2.0 * tle.argpo)))
        tle.cc5 = 2.0 * coef1 * ao * omeosq * (1.0 + 2.75 *
                    (etasq + eeta) + eeta * etasq)
        cosio4 = cosio2 * cosio2
        temp1  = 1.5 * j2 * pinvsq * tle.no
        temp2  = 0.5 * temp1 * j2 * pinvsq
        temp3  = -0.46875 * j4 * pinvsq * pinvsq * tle.no
        tle.mdot     = tle.no + 0.5 * temp1 * rteosq * tle.con41 + 0.0625 *
                        temp2 * rteosq * (13.0 - 78.0 * cosio2 + 137.0 * cosio4)
        tle.argpdot  = (-0.5 * temp1 * con42 + 0.0625 * temp2 *
                            (7.0 - 114.0 * cosio2 + 395.0 * cosio4) +
                            temp3 * (3.0 - 36.0 * cosio2 + 49.0 * cosio4))
        xhdot1            = -temp1 * cosio
        tle.nodedot = xhdot1 + (0.5 * temp2 * (4.0 - 19.0 * cosio2) +
                            2.0 * temp3 * (3.0 - 7.0 * cosio2)) * cosio
        xpidot            =  tle.argpdot+ tle.nodedot
        tle.omgcof   = tle.bstar * cc3 * cos(tle.argpo)
        tle.xmcof    = 0.0
        if tle.ecco > 1.0e-4
            tle.xmcof = -x2o3 * coef * tle.bstar / eeta
        end
        tle.nodecf = 3.5 * omeosq * xhdot1 * tle.cc1
        tle.t2cof   = 1.5 * tle.cc1
        #  sgp4fix for divide by zero with xinco = 180 deg
        if abs(cosio+1.0) > 1.5e-12
            tle.xlcof = -0.25 * j3oj2 * sinio * (3.0 + 5.0 * cosio) / (1.0 + cosio)
        else
            tle.xlcof = -0.25 * j3oj2 * sinio * (3.0 + 5.0 * cosio) / temp4
        end
        tle.aycof   = -0.5 * j3oj2 * sinio
        #  sgp4fix use multiply for speed instead of pow
        delmotemp = 1.0 + tle.eta * cos(tle.mo)
        tle.delmo   = delmotemp * delmotemp * delmotemp
        tle.sinmao  = sin(tle.mo)
        tle.x7thm1  = 7.0 * cosio2 - 1.0

        #  --------------- deep space initialization -------------
        if 2*pi / tle.no >= 225.0

            tle.method = 'd'
            tle.isimp  = 1
            tc    =  0.0
            inclm = tle.inclo

            (
                snodm, cnodm, sinim,  cosim, sinomm,
                cosomm,day,   tle.e3,     tle.ee2,   em,
                emsq,  gam,   tle.peo,    tle.pgho,  tle.pho,
                tle.pinco, tle.plo,   rtemsq, tle.se2,   tle.se3,
                tle.sgh2, tle.sgh3, tle.sgh4, tle.sh2, tle.sh3,
                tle.si2,  tle.si3,  tle.sl2, tle.sl3, tle.sl4,
                s1,    s2,    s3,     s4,    s5,
                s6,    s7,    ss1,    ss2,   ss3,
                ss4,   ss5,   ss6,    ss7,   sz1,
                sz2,   sz3,   sz11,   sz12,  sz13,
                sz21,  sz22,  sz23,   sz31,  sz32,
                sz33,  tle.xgh2, tle.xgh3, tle.xgh4, tle.xh2,
                tle.xh3, tle.xi2, tle.xi3, tle.xl2, tle.xl3,
                tle.xl4,   nm,    z1,     z2,    z3,
                z11,   z12,   z13,    z21,   z22,
                z23,   z31,   z32,    z33,   tle.zmol,
                tle.zmos
            ) = dscom(
                epoch, tle.ecco, tle.argpo, tc, tle.inclo, tle.nodeo,
                tle.no,
                tle.e3, tle.ee2,
                tle.peo,  tle.pgho,   tle.pho, tle.pinco,
                tle.plo,        tle.se2, tle.se3,
                tle.sgh2, tle.sgh3,   tle.sgh4,
                tle.sh2,  tle.sh3,    tle.si2, tle.si3,
                tle.sl2,  tle.sl3,    tle.sl4,
                tle.xgh2, tle.xgh3,   tle.xgh4, tle.xh2,
                tle.xh3,  tle.xi2,    tle.xi3,  tle.xl2,
                tle.xl3,  tle.xl4,
                tle.zmol, tle.zmos
            )

            (tle.ecco, tle.inclo, tle.nodeo, tle.argpo, tle.mo
            ) = dpper(
                tle, inclm, tle.init,
                tle.ecco, tle.inclo, tle.nodeo, tle.argpo, tle.mo,
                tle.afspc_mode
            )

            argpm  = 0.0
            nodem  = 0.0
            mm     = 0.0

            (
                em,    argpm,  inclm, mm,
                nm,    nodem,
                tle.irez,  tle.atime,
                tle.d2201, tle.d2211, tle.d3210, tle.d3222,
                tle.d4410, tle.d4422, tle.d5220, tle.d5232,
                tle.d5421, tle.d5433, tle.dedt,  tle.didt,
                tle.dmdt, dndt,  tle.dnodt, tle.domdt,
                tle.del1,  tle.del2,  tle.del3,  tle.xfact,
                tle.xlamo, tle.xli,   tle.xni
            ) = dsinit(
                whichconst,
                cosim, emsq, tle.argpo, s1, s2, s3, s4, s5, sinim, ss1, ss2, ss3, ss4,
                ss5, sz1, sz3, sz11, sz13, sz21, sz23, sz31, sz33, tle.t, tc,
                tle.gsto, tle.mo, tle.mdot, tle.no, tle.nodeo,
                tle.nodedot, xpidot, z1, z3, z11, z13, z21, z23, z31, z33,
                tle.ecco, eccsq, em, argpm, inclm, mm, nm, nodem,
                tle.irez,  tle.atime,
                tle.d2201, tle.d2211, tle.d3210, tle.d3222 ,
                tle.d4410, tle.d4422, tle.d5220, tle.d5232,
                tle.d5421, tle.d5433, tle.dedt,  tle.didt,
                tle.dmdt,  tle.dnodt, tle.domdt,
                tle.del1,  tle.del2,  tle.del3,  tle.xfact,
                tle.xlamo, tle.xli,   tle.xni
            )
        end

        #----------- set variables if not deep space -----------
        if tle.isimp != 1

            cc1sq          = tle.cc1 * tle.cc1
            tle.d2    = 4.0 * ao * tsi * cc1sq
            temp           = tle.d2 * tsi * tle.cc1 / 3.0
            tle.d3    = (17.0 * ao + sfour) * temp
            tle.d4    = 0.5 * temp * ao * tsi * (221.0 * ao + 31.0 * sfour) *
                            tle.cc1
            tle.t3cof = tle.d2 + 2.0 * cc1sq
            tle.t4cof = 0.25 * (3.0 * tle.d3 + tle.cc1 *
                            (12.0 * tle.d2 + 10.0 * cc1sq))
            tle.t5cof = 0.2 * (3.0 * tle.d4 +
                            12.0 * tle.cc1 * tle.d3 +
                            6.0 * tle.d2 * tle.d2 +
                            15.0 * cc1sq * (2.0 * tle.d2 + cc1sq))
        end
    end

    """
    /* finally propogate to zero epoch to initialize all others. */
    // sgp4fix take out check to let satellites process until they are actually below earth surface
    //       if(tle.error == 0)
    """
    # sgp4(tle, 0.0)

    tle.init = 'n'

    # sgp4fix return boolean. tle.error contains any error codes
    return true
end

export sgp4
"""
/*-----------------------------------------------------------------------------
*
*                             procedure sgp4
*
*  this procedure is the sgp4 prediction model from space command. this is an
*    updated and combined version of sgp4 and sdp4, which were originally
*    published separately in spacetrack report #3. this version follows the
*    methodology from the aiaa paper (2006) describing the history and
*    development of the code.
*
*  author        : david vallado                  719-573-2600   28 jun 2005
*
*  inputs        :
*    tle	 - initialised structure from sgp4init() call.
*    tsince	 - time eince epoch (minutes)
*
*  outputs       :
*    r           - position vector                     km
*    v           - velocity                            km/sec
*  return code - non-zero on error.
*                   1 - mean elements, ecc >= 1.0 or ecc < -0.001 or a < 0.95 er
*                   2 - mean motion less than 0.0
*                   3 - pert elements, ecc < 0.0  or  ecc > 1.0
*                   4 - semi-latus rectum < 0.0
*                   5 - epoch elements are sub-orbital
*                   6 - satellite has decayed
*
*  locals        :
*    am          -
*    axnl, aynl        -
*    betal       -
*    cosim   , sinim   , cosomm  , sinomm  , cnod    , snod    , cos2u   ,
*    sin2u   , coseo1  , sineo1  , cosi    , sini    , cosip   , sinip   ,
*    cosisq  , cossu   , sinsu   , cosu    , sinu
*    delm        -
*    delomg      -
*    dndt        -
*    eccm        -
*    emsq        -
*    ecose       -
*    el2         -
*    eo1         -
*    eccp        -
*    esine       -
*    argpm       -
*    argpp       -
*    omgadf      -
*    pl          -
*    r           -
*    rtemsq      -
*    rdotl       -
*    rl          -
*    rvdot       -
*    rvdotl      -
*    su          -
*    t2  , t3   , t4    , tc
*    tem5, temp , temp1 , temp2  , tempa  , tempe  , templ
*    u   , ux   , uy    , uz     , vx     , vy     , vz
*    inclm       - inclination
*    mm          - mean anomaly
*    nm          - mean motion
*    nodem       - right asc of ascending node
*    xinc        -
*    xincp       -
*    xl          -
*    xlm         -
*    mp          -
*    xmdf        -
*    xmx         -
*    xmy         -
*    nodedf      -
*    xnode       -
*    nodep       -
*    np          -
*
*  coupling      :
*    getgravconst-
*    dpper
*    dpspace
*
*  references    :
*    hoots, roehrich, norad spacetrack report #3 1980
*    hoots, norad spacetrack report #6 1986
*    hoots, schumacher and glover 2004
*    vallado, crawford, hujsak, kelso  2006
  ----------------------------------------------------------------------------*/
"""
function sgp4(tle::TLE, tsince::Real)

    # Initialize Values
    mrt = 0.0
    whichconst = tle.whichconst

    rv = zeros(Float64, 6)

    """
    /* ------------------ set mathematical constants --------------- */
    // sgp4fix divisor for divide by zero check on inclination
    // the old check used 1.0 + cos(pi-1.0e-9), but then compared it to
    // 1.5 e-12, so the threshold was changed to 1.5e-12 for consistency
    """
    temp4 =   1.5e-12
    twopi = 2.0 * pi
    x2o3  = 2.0 / 3.0
    #  sgp4fix identify constants and allow alternate values
    tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2 = getgravconst(whichconst)
    vkmpersec     = radiusearthkm * xke/60.0

    #  --------------------- clear sgp4 error flag -----------------
    tle.t     = tsince
    tle.error = 0

    #  ------- update for secular gravity and atmospheric drag -----
    xmdf    = tle.mo + tle.mdot * tle.t
    argpdf  = tle.argpo + tle.argpdot * tle.t
    nodedf  = tle.nodeo + tle.nodedot * tle.t
    argpm   = argpdf
    mm      = xmdf
    t2      = tle.t * tle.t
    nodem   = nodedf + tle.nodecf * t2
    tempa   = 1.0 - tle.cc1 * tle.t
    tempe   = tle.bstar * tle.cc4 * tle.t
    templ   = tle.t2cof * t2

    if tle.isimp != 1
        delomg = tle.omgcof * tle.t
        #  sgp4fix use mutliply for speed instead of pow
        delmtemp =  1.0 + tle.eta * cos(xmdf)
        delm   = tle.xmcof *
                (delmtemp * delmtemp * delmtemp -
                tle.delmo)
        temp   = delomg + delm
        mm     = xmdf + temp
        argpm  = argpdf - temp
        t3     = t2 * tle.t
        t4     = t3 * tle.t
        tempa  = tempa - tle.d2 * t2 - tle.d3 * t3 -
                        tle.d4 * t4
        tempe  = tempe + tle.bstar * tle.cc5 * (sin(mm) -
                        tle.sinmao)
        templ  = templ + tle.t3cof * t3 + t4 * (tle.t4cof +
                        tle.t * tle.t5cof)
    end

    nm    = tle.no
    em    = tle.ecco
    inclm = tle.inclo
    if tle.method == 'd'

        tc = tle.t
        (
            atime, em,    argpm,  inclm, xli,
            mm,    xni,   nodem,  dndt,  nm,
        ) = dspace(
            tle.irez,
            tle.d2201, tle.d2211, tle.d3210,
            tle.d3222, tle.d4410, tle.d4422,
            tle.d5220, tle.d5232, tle.d5421,
            tle.d5433, tle.dedt,  tle.del1,
            tle.del2,  tle.del3,  tle.didt,
            tle.dmdt,  tle.dnodt, tle.domdt,
            tle.argpo, tle.argpdot, tle.t, tc,
            tle.gsto, tle.xfact, tle.xlamo,
            tle.no, tle.atime,
            em, argpm, inclm, tle.xli, mm, tle.xni,
            nodem, nm
        )
    end

    if nm <= 0.0

        @error "mean motion $nm is less than zero"
        tle.error = 2

        #  sgp4fix add return
        return false
    end

    am = (xke / nm) ^ x2o3 * tempa * tempa;
    nm = (xke / am) ^ 1.5
    em = em - tempe

    #  fix tolerance for error recognition
    #  sgp4fix am is fixed from the previous nm check
    if em >= 1.0 || em < -0.001  # || (am < 0.95)

        @error "mean eccentricity $em not within range 0.0 <= e < 1.0"
        tle.error = 1

        #  sgp4fix to return if there is an error in eccentricity
        return false
    end

    #  sgp4fix fix tolerance to avoid a divide by zero
    if em < 1.0e-6
        em  = 1.0e-6
    end

    mm     = mm + tle.no * templ
    xlm    = mm + argpm + nodem
    emsq   = em * em
    temp   = 1.0 - emsq

    nodem  = nodem >= 0.0 ? (nodem % twopi) : -(-nodem % twopi)
    argpm  = argpm % twopi
    xlm    = xlm % twopi
    mm     = (xlm - argpm - nodem) % twopi

    #  ----------------- compute extra mean quantities -------------
    sinim = sin(inclm)
    cosim = cos(inclm)

    #  -------------------- add lunar-solar periodics --------------
    ep     = em
    xincp  = inclm
    argpp  = argpm
    nodep  = nodem
    mp     = mm
    sinip  = sinim
    cosip  = cosim
    if tle.method == 'd'

        ep, xincp, nodep, argpp, mp = dpper(
            tle, tle.inclo,
            'n', ep, xincp, nodep, argpp, mp, tle.afspc_mode
            )

        if xincp < 0.0

            xincp  = -xincp
            nodep = nodep + pi
            argpp  = argpp - pi
        end

        if ep < 0.0 || ep > 1.0

            @error "perturbed eccentricity $ep not within range 0.0 <= e <= 1.0"
            tle.error = 3

            return false
        end
    end

    #  -------------------- long period periodics ------------------
    if tle.method == 'd'
        sinip =  sin(xincp)
        cosip =  cos(xincp)
        tle.aycof = -0.5*j3oj2*sinip
        #  sgp4fix for divide by zero for xincp = 180 deg
        if abs(cosip+1.0) > 1.5e-12
            tle.xlcof = -0.25 * j3oj2 * sinip * (3.0 + 5.0 * cosip) / (1.0 + cosip)
        else
            tle.xlcof = -0.25 * j3oj2 * sinip * (3.0 + 5.0 * cosip) / temp4
        end
    end

    axnl = ep * cos(argpp)
    temp = 1.0 / (am * (1.0 - ep * ep))
    aynl = ep* sin(argpp) + temp * tle.aycof
    xl   = mp + argpp + nodep + temp * tle.xlcof * axnl

    #  --------------------- solve kepler's equation ---------------
    u    = (xl - nodep) % twopi
    eo1  = u
    tem5 = 9999.9
    ktr = 1
    sineo1 = 0.0
    coseo1 = 0.0
    #    sgp4fix for kepler iteration
    #    the following iteration needs better limits on corrections
    while abs(tem5) >= 1.0e-12 && ktr <= 10

        sineo1 = sin(eo1)
        coseo1 = cos(eo1)
        tem5   = 1.0 - coseo1 * axnl - sineo1 * aynl
        tem5   = (u - aynl * coseo1 + axnl * sineo1 - eo1) / tem5

        if abs(tem5) >= 0.95
            tem5 = tem5 > 0.0 ? 0.95 : -0.95
        end

        eo1 = eo1 + tem5
        ktr = ktr + 1
    end

    #  ------------- short period preliminary quantities -----------
    ecose = axnl*coseo1 + aynl*sineo1
    esine = axnl*sineo1 - aynl*coseo1
    el2   = axnl*axnl + aynl*aynl
    pl    = am*(1.0-el2)

    if pl < 0.0
        @error "semilatus rectum $pl is less than zero"
        tle.error = 4

        return false
    else

        rl     = am * (1.0 - ecose)
        rdotl  = sqrt(am) * esine/rl
        rvdotl = sqrt(pl) / rl
        betal  = sqrt(1.0 - el2)
        temp   = esine / (1.0 + betal)
        sinu   = am / rl * (sineo1 - aynl - axnl * temp)
        cosu   = am / rl * (coseo1 - axnl + aynl * temp)
        su     = atan(sinu, cosu)
        sin2u  = (cosu + cosu) * sinu
        cos2u  = 1.0 - 2.0 * sinu * sinu
        temp   = 1.0 / pl
        temp1  = 0.5 * j2 * temp
        temp2  = temp1 * temp

        #  -------------- update for short period periodics ------------
        if tle.method == 'd'

            cosisq     = cosip * cosip
            tle.con41  = 3.0*cosisq - 1.0
            tle.x1mth2 = 1.0 - cosisq
            tle.x7thm1 = 7.0*cosisq - 1.0
        end

        mrt   = rl * (1.0 - 1.5 * temp2 * betal * tle.con41) +
                0.5 * temp1 * tle.x1mth2 * cos2u
        su    = su - 0.25 * temp2 * tle.x7thm1 * sin2u
        xnode = nodep + 1.5 * temp2 * cosip * sin2u
        xinc  = xincp + 1.5 * temp2 * cosip * sinip * cos2u
        mvt   = rdotl - nm * temp1 * tle.x1mth2 * sin2u / xke
        rvdot = rvdotl + nm * temp1 * (tle.x1mth2 * cos2u +
                1.5 * tle.con41) / xke

        #  --------------------- orientation vectors -------------------
        sinsu =  sin(su)
        cossu =  cos(su)
        snod  =  sin(xnode)
        cnod  =  cos(xnode)
        sini  =  sin(xinc)
        cosi  =  cos(xinc)
        xmx   = -snod * cosi
        xmy   =  cnod * cosi
        ux    =  xmx * sinsu + cnod * cossu
        uy    =  xmy * sinsu + snod * cossu
        uz    =  sini * sinsu
        vx    =  xmx * cossu - cnod * sinsu
        vy    =  xmy * cossu - snod * sinsu
        vz    =  sini * cossu

        #  --------- position and velocity (in km and km/sec) ----------
        _mr = mrt * radiusearthkm
        rv[1] = _mr * ux
        rv[2] = _mr * uy
        rv[3] = _mr * uz
        rv[4] = (mvt * ux + rvdot * vx) * vkmpersec
        rv[5] = (mvt * uy + rvdot * vy) * vkmpersec
        rv[6] = (mvt * uz + rvdot * vz) * vkmpersec
    end


    #  sgp4fix for decaying satellites
    if mrt < 1.0
        @error "mrt $mrt is less than 1.0 indicating the satellite has decayed"
        tle.error = 6
        return false
    end

    return rv
end

function sgp4(tle::TLE, epc::Epoch)
    return sgp4(tle, (epc - tle.epoch)/60.0)
end

# Constructor to initialize from just strings
function TLE(line1::String, line2::String, whichconst::gravconsttype=wgs84, afspc_mode::Bool=false)
    # Validate input lines
    if length(line1) != 69
        throw(ArgumentError("Invalid line 1 length. Expecting 69 found $(length(line1))"))
    end

    if length(line2) != 69
        throw(ArgumentError("Invalid line 2 length. Expecting 69 found $(length(line2))"))
    end

    csum = tle_checksum(line1[1:end-1])
    if csum != parse(Int, line1[end])
        throw(ArgumentError("Invalid checksum for line 1. Expecting $csum found $(line1[end])"))
    end

    csum = tle_checksum(line2[1:end-1])
    if csum != parse(Int, line2[end])
        throw(ArgumentError("Invalid checksum for line 2. Expecting $csum found $(line2[end])"))
    end

    # --- Initialize Variable Values and Extract TLE Data --- #

    # Computation values
    init  = 'n'
    error = 0

    # Initialize Generic Values
    satnum         = 0
    classification = 'U'
    intldesg       = ""
    epochyr        = 0
    epoch          = Epoch(2000, 1, 1, tsys="UTC")
    elnum          = 0
    revnum         = 0
    method         = 'n'

    # Initialize Common values
    a = 0.0
    altp = 0.0
    alta = 0.0
    epochdays = 0.0
    jdsatepoch = 0.0
    nddot = 0.0
    ndot = 0.0
    bstar = 0.0
    rcse = 0.0
    inclo = 0.0
    nodeo = 0.0
    ecco = 0.0
    argpo = 0.0
    mo = 0.0
    no = 0.0

    # Near Earth Values
    isimp = 0
    aycof, con41, cc1, cc4, cc5 = 0.0, 0.0, 0.0, 0.0, 0.0
    d2, d3, d4, delmo = 0.0, 0.0, 0.0, 0.0
    eta, argpdot, omgcof, sinmao = 0.0, 0.0, 0.0, 0.0
    t, t2cof, t3cof, t4cof, t5cof = 0.0, 0.0, 0.0, 0.0, 0.0
    x1mth2, x7thm1, mdot, nodedot = 0.0, 0.0, 0.0, 0.0
    xlcof, xmcof, nodecf = 0.0, 0.0, 0.0

    # Deep space values 
    irez = 0
    d2201, d2211, d3210, d3222, d4410 = 0.0, 0.0, 0.0, 0.0, 0.0
    d4422, d5220, d5232, d5421, d5433 = 0.0, 0.0, 0.0, 0.0, 0.0
    dedt, del1, del2, del3, = 0.0, 0.0, 0.0, 0.0, 0.0
    didt, dmdt, dnodt, domdt = 0.0, 0.0, 0.0, 0.0
    e3, ee2, peo, pgho, pho, pinco, plo = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    se2, se3, sgh2, sgh3, sgh4, sh2, sh3 = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    si2, si3, sl2, sl3, sl4, gsto, xfact = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    xgh2, xgh3, xgh4, xh2 = 0.0, 0.0, 0.0, 0.0
    xh3, xi2, xi3, xl2, xl3, xl4 = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    xlamo, zmol, zmos, atime, xli, xni = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0

    # --- Create TLE Object --- #
    tle = TLE(line1, line2, whichconst, afspc_mode, init,
        satnum, classification, intldesg, epochyr, epoch, elnum, revnum, method,
        a, altp, alta, epochdays, jdsatepoch, nddot, ndot, bstar, rcse, inclo, 
        nodeo, ecco, argpo, mo, no, 
        isimp, aycof, con41, cc1, cc4, cc5, d2, d3, d4, delmo, eta, argpdot, 
        omgcof, sinmao, t, t2cof, t3cof, t4cof, t5cof, x1mth2, x7thm1, mdot, 
        nodedot, xlcof, xmcof, nodecf, 
        irez, d2201, d2211, d3210, d3222, d4410, d4422, d5220, d5232, d5421, 
        d5433, dedt, del1, del2, del3, didt, dmdt, dnodt, domdt, e3, ee2, peo, 
        pgho, pho, pinco, plo, se2, se3, sgh2, sgh3, sgh4, sh2, sh3, si2, si3, 
        sl2, sl3, sl4, gsto, xfact, xgh2, xgh3, xgh4, xh2, xh3, xi2, xi3, xl2, 
        xl3, xl4, xlamo, zmol, zmos, atime, xli, xni, error)

    # Extract information from line 1
    tle.satnum         = parse(Int, line1[3:7])
    tle.classification = line1[8]
    tle.intldesg       = line1[10:17]
    tle.epochyr        = parse(Int, line1[19:20])
    tle.epochdays      = parse(Float64, line1[21:32])
    tle.ndot           = parse(Float64, line1[34] * "0" * lstrip(line1[35:43]))
    tle.nddot          = parse(Float64, line1[45] * "0." * lstrip(line1[46:50]))*10.0^parse(Int, line1[51:52])
    tle.bstar          = parse(Float64, line1[54] * "0." * lstrip(line1[55:59]))*10.0^parse(Int, line1[60:61])
    tle.elnum          = parse(Int, line1[65:68])
 
 
    # Extract values from line 2
    tle.inclo  = parse(Float64, line2[9:16])
    tle.nodeo  = parse(Float64, line2[18:25])
    tle.ecco   = parse(Float64, "0." * line2[27:33])
    tle.argpo  = parse(Float64, line2[35:42])
    tle.mo     = parse(Float64, line2[44:51])
    tle.no     = parse(Float64, line2[53:63])
    tle.revnum = parse(Int, line2[64:68])
 
    #  ---- find no, ndot, nddot ----
    tle.no   = tle.no / xpdotp #   rad/min
    # tle.nddot= tle.nddot * 10.0 ^ nexp  # Already done during parse
    # tle.bstar= tle.bstar * 10.0 ^ ibexp # Already done during parse

    #  ---- convert to sgp4 units ----
    tumin, _ = getgravconst(whichconst)
    tle.a    = tle.no*tumin ^ (-2.0/3.0)
    tle.ndot = tle.ndot  / (xpdotp*1440.0)  #   ? * minperday
    tle.nddot= tle.nddot / (xpdotp*1440.0*1440)

    #  ---- find standard orbital elements ----
    tle.inclo = tle.inclo  * DEG2RAD
    tle.nodeo = tle.nodeo  * DEG2RAD
    tle.argpo = tle.argpo  * DEG2RAD
    tle.mo    = tle.mo     * DEG2RAD

    tle.alta = tle.a*(1.0 + tle.ecco) - 1.0
    tle.altp = tle.a*(1.0 - tle.ecco) - 1.0

    """
    // ----------------------------------------------------------------
    // find sgp4epoch time of element set
    // remember that sgp4 uses units of days from 0 jan 1950 (sgp4epoch)
    // and minutes from the epoch (time)
    // ----------------------------------------------------------------

    // ---------------- temp fix for years from 1957-2056 -------------------
    // --------- correct fix will occur when year is 4-digit in tle ---------
    """
    year = 0
    if tle.epochyr < 57
        year = tle.epochyr + 2000
    else
        year = tle.epochyr + 1900
    end

    mon,day,hr,minute,sec = days2mdhms(year, tle.epochdays)

    tle.epochyr = year
    tle.jdsatepoch = jday(year,mon,day,hr,minute,sec)
    tle.epoch = Epoch(year, mon, day, hr, minute, sec, tsys="UTC")
    #  ---------------- initialize the orbit at sgp4epoch -------------------
    sgp4init(whichconst, afspc_mode, tle.satnum, tle.jdsatepoch-2433281.5, tle.bstar,
             tle.ecco, tle.argpo, tle.inclo, tle.mo, tle.no,
             tle.nodeo, tle)

    return tle
end

################
# SGP4 Wrapper #
################

# Wrapper functions to make calling into TLE propagator easy and compatible with
# Rest of package reference frames. 

export state
"""
Compute the satellite state at the given epoch as output from the SGP4 propagator.

Arguments:
- `tle::TLE` Two-Line element object 
- `epc::Epoch` Epoch for SGP4 state
- `opscode::Char` Operatitonal mode of propagators

Returns:
- `rv::Array{Float64, 1}` Position and velocity as output by TLE propagator. 
    Units [m m/s]
"""
function state(tle::TLE, epc::Epoch)
    rv = sgp4(tle, epc)

    for i in 1:length(rv)
        rv[i] = rv[i]*1000
    end

    return rv
end

export ecef
```
Compute the satellite state in the Earth-fixed frame.
```
function ecef(tle::TLE, epc::Epoch)
    s    = state(tle, epc)
    ecef = zeros(Float64, 6)
    gst  = gstime(epc - tle.epoch + 2433281.5)
    omega_vec = [0, 0, OMEGA_EARTH]
    R = Rz(gst)
    ecef[1:3] = R * s[1:3]
    ecef[4:6] = R * s[4:6] - cross(omega_vec, R * s[1:3])

    return ecef
end

export eci
```
Compute the satellite state in the Earth-centered inertial.
```
function eci(tle::TLE, epc::Epoch)
    s = ecef(tle, epc)

    return sECEFtoECI(epc, s)
end