var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "#SatelliteDynamics.jl-Documentation-1",
    "page": "Home",
    "title": "SatelliteDynamics.jl Documentation",
    "category": "section",
    "text": "Welcome to the SatelliteDynamics.jl documentation! The package is meant to address the needs of the satellite operator, academic researcher, and public enthusiast  communities. In particular, it is designed and built to fill the following needs:Open-source, MIT-lencesed software to remove barriers to entry of doing high-fidelity satellite modeling.\nHigh-quality, validated, tested library for orbit and attitude dynamics modeling.\nEasily acceible API design to make implementation of simulation and analysis intuitive."
},

{
    "location": "#Getting-Started:-Installation-And-First-Steps-1",
    "page": "Home",
    "title": "Getting Started: Installation And First Steps",
    "category": "section",
    "text": "To install the package, use the following command inside the Julia REPL:Pkg.add(\"SatelliteDynamics\")To load the package, use the command:using SatelliteDynamics"
},

{
    "location": "#Package-Structure-1",
    "page": "Home",
    "title": "Package Structure",
    "category": "section",
    "text": "The package is divided into a number of submodules each designed to provide a single, well-defined set of functions. The details on Pages = [\n    \"modules/constants.md\",\n    \"modules/universe.md\",\n]\nDepth = 2"
},

{
    "location": "#Examples-1",
    "page": "Home",
    "title": "Examples",
    "category": "section",
    "text": "The best way to learn how to use any software is to Pages = [\n    \"tutorials/universe_example.md\",\n    \"tutorials/time_example.md\",\n    \"tutorials/epoch_example.md\",\n]\nDepth = 2"
},

{
    "location": "#Citing-1",
    "page": "Home",
    "title": "Citing",
    "category": "section",
    "text": "The software in this package was developed as part of academic research. If you would like to help support it, please star the repository as such metrics on usage help guage interest and secure funding. If you use the software as part of your research, teaching, or other activities, we would be grateful if you could cite our work."
},

{
    "location": "modules/constants/#",
    "page": "Constants",
    "title": "Constants",
    "category": "page",
    "text": ""
},

{
    "location": "modules/constants/#SatelliteDynamics.Constants.AS2RAD",
    "page": "Constants",
    "title": "SatelliteDynamics.Constants.AS2RAD",
    "category": "constant",
    "text": "Constant to convert arcseconds to radians. Equal to 2pi/(360*3600). [rad/as]\n\n\n\n\n\n"
},

{
    "location": "modules/constants/#SatelliteDynamics.Constants.RAD2AS",
    "page": "Constants",
    "title": "SatelliteDynamics.Constants.RAD2AS",
    "category": "constant",
    "text": "Constant to convert radians to arcseconds. Equal to 2pi/(360*3600) [as/ras]\n\n\n\n\n\n"
},

{
    "location": "modules/constants/#SatelliteDynamics.Constants.C_LIGHT",
    "page": "Constants",
    "title": "SatelliteDynamics.Constants.C_LIGHT",
    "category": "constant",
    "text": "Speed of light in vacuum. [m/s]\n\nD. Vallado, Fundamentals of Astrodynamics and Applications (4th Ed.), 2010\n\n\n\n\n\n"
},

{
    "location": "modules/constants/#SatelliteDynamics.Constants.AU",
    "page": "Constants",
    "title": "SatelliteDynamics.Constants.AU",
    "category": "constant",
    "text": "Astronomical Unit. Equal to the mean distance of the Earth from the sun. TDB-compatible value. [m]\n\nP. Gérard and B. Luzum, IERS Technical Note 36, 2010\n\n\n\n\n\n"
},

{
    "location": "modules/constants/#SatelliteDynamics.Constants.MJD_ZERO",
    "page": "Constants",
    "title": "SatelliteDynamics.Constants.MJD_ZERO",
    "category": "constant",
    "text": "Offset of Modified Julian Days representation with respect to Julian Days. For  a time, t, MJD_ZERO is equal to:\n\nMJD_ZERO = t_jd - t_mjd\n\nWhere tjd is the epoch represented in Julian Days, and tmjd is the epoch in Modified Julian Days.\n\nO. Montenbruck, and E. Gill, Satellite Orbits: Models, Methods and  Applications, 2012.\n\n\n\n\n\n"
},

{
    "location": "modules/constants/#SatelliteDynamics.Constants.MJD2000",
    "page": "Constants",
    "title": "SatelliteDynamics.Constants.MJD2000",
    "category": "constant",
    "text": "Modified Julian Date of January 1, 2000 00:00:00. Value is independent of time scale.\n\nO. Montenbruck, and E. Gill, Satellite Orbits: Models, Methods and  Applications, 2012.\n\n\n\n\n\n"
},

{
    "location": "modules/constants/#SatelliteDynamics.Constants.GPS_TAI",
    "page": "Constants",
    "title": "SatelliteDynamics.Constants.GPS_TAI",
    "category": "constant",
    "text": "Offset of GPS time system with respect to TAI time system.\n\nO. Montenbruck, and E. Gill, Satellite Orbits: Models, Methods and  Applications, 2012.\n\n\n\n\n\n"
},

{
    "location": "modules/constants/#SatelliteDynamics.Constants.TAI_GPS",
    "page": "Constants",
    "title": "SatelliteDynamics.Constants.TAI_GPS",
    "category": "constant",
    "text": "Offset of TAI time system with respect to GPS time system.\n\nO. Montenbruck, and E. Gill, Satellite Orbits: Models, Methods and  Applications, 2012.\n\n\n\n\n\n"
},

{
    "location": "modules/constants/#SatelliteDynamics.Constants.TT_TAI",
    "page": "Constants",
    "title": "SatelliteDynamics.Constants.TT_TAI",
    "category": "constant",
    "text": "Offset of TT time system with respect to TAI time system.\n\nO. Montenbruck, and E. Gill, Satellite Orbits: Models, Methods and  Applications, 2012.\n\n\n\n\n\n"
},

{
    "location": "modules/constants/#SatelliteDynamics.Constants.TAI_TT",
    "page": "Constants",
    "title": "SatelliteDynamics.Constants.TAI_TT",
    "category": "constant",
    "text": "Offset of TAI time system with respect to TT time system.\n\nO. Montenbruck, and E. Gill, Satellite Orbits: Models, Methods and  Applications, 2012.\n\n\n\n\n\n"
},

{
    "location": "modules/constants/#SatelliteDynamics.Constants.GPS_TT",
    "page": "Constants",
    "title": "SatelliteDynamics.Constants.GPS_TT",
    "category": "constant",
    "text": "Offset of GPS time system with respect to TT time system.\n\nO. Montenbruck, and E. Gill, Satellite Orbits: Models, Methods and  Applications, 2012.\n\n\n\n\n\n"
},

{
    "location": "modules/constants/#SatelliteDynamics.Constants.TT_GPS",
    "page": "Constants",
    "title": "SatelliteDynamics.Constants.TT_GPS",
    "category": "constant",
    "text": "Offset of TT time system with respect to GPS time system.\n\nO. Montenbruck, and E. Gill, Satellite Orbits: Models, Methods and  Applications, 2012.\n\n\n\n\n\n"
},

{
    "location": "modules/constants/#SatelliteDynamics.Constants.GPS_ZERO",
    "page": "Constants",
    "title": "SatelliteDynamics.Constants.GPS_ZERO",
    "category": "constant",
    "text": "Modified Julian Date of the start of the GPS time system in the GPS time system. This date was January 6, 1980 0H as reckond in the UTC time system.\n\nO. Montenbruck, and E. Gill, Satellite Orbits: Models, Methods and  Applications, 2012.\n\n\n\n\n\n"
},

{
    "location": "modules/constants/#SatelliteDynamics.Constants.R_EARTH",
    "page": "Constants",
    "title": "SatelliteDynamics.Constants.R_EARTH",
    "category": "constant",
    "text": "Earth\'s equatorial radius. [m]\n\nGGM05s Gravity Model\n\n\n\n\n\n"
},

{
    "location": "modules/constants/#SatelliteDynamics.Constants.WGS84_a",
    "page": "Constants",
    "title": "SatelliteDynamics.Constants.WGS84_a",
    "category": "constant",
    "text": "Earth\'s semi-major axis as defined by the WGS84 geodetic system. [m]\n\nNIMA Technical Report TR8350.2\n\n\n\n\n\n"
},

{
    "location": "modules/constants/#SatelliteDynamics.Constants.WGS84_f",
    "page": "Constants",
    "title": "SatelliteDynamics.Constants.WGS84_f",
    "category": "constant",
    "text": "Earth\'s ellipsoidal flattening.  WGS84 Value.\n\nNIMA Technical Report TR8350.2\n\n\n\n\n\n"
},

{
    "location": "modules/constants/#SatelliteDynamics.Constants.GM_EARTH",
    "page": "Constants",
    "title": "SatelliteDynamics.Constants.GM_EARTH",
    "category": "constant",
    "text": "Earth\'s Gravitational constant [m^3/s^2]\n\nO. Montenbruck, and E. Gill, Satellite Orbits: Models, Methods and  Applications, 2012.\n\n\n\n\n\n"
},

{
    "location": "modules/constants/#SatelliteDynamics.Constants.e_EARTH",
    "page": "Constants",
    "title": "SatelliteDynamics.Constants.e_EARTH",
    "category": "constant",
    "text": "Earth\'s first eccentricity. WGS84 Value. [dimensionless]\n\nNIMA Technical Report TR8350.2\n\n\n\n\n\n"
},

{
    "location": "modules/constants/#SatelliteDynamics.Constants.J2_EARTH",
    "page": "Constants",
    "title": "SatelliteDynamics.Constants.J2_EARTH",
    "category": "constant",
    "text": "Earth\'s first zonal harmonic. [dimensionless]\n\nGGM05s Gravity Model.\n\n\n\n\n\n"
},

{
    "location": "modules/constants/#SatelliteDynamics.Constants.OMEGA_EARTH",
    "page": "Constants",
    "title": "SatelliteDynamics.Constants.OMEGA_EARTH",
    "category": "constant",
    "text": "Earth axial rotation rate. [rad/s]\n\nD. Vallado, Fundamentals of Astrodynamics and Applications (4th Ed.), p. 222, 2010\n\n\n\n\n\n"
},

{
    "location": "modules/constants/#SatelliteDynamics.Constants.GM_SUN",
    "page": "Constants",
    "title": "SatelliteDynamics.Constants.GM_SUN",
    "category": "constant",
    "text": "Gravitational constant of the Sun. [m^3/s^2]\n\nO. Montenbruck, and E. Gill, Satellite Orbits: Models, Methods and  Applications, 2012.\n\n\n\n\n\n"
},

{
    "location": "modules/constants/#SatelliteDynamics.Constants.R_SUN",
    "page": "Constants",
    "title": "SatelliteDynamics.Constants.R_SUN",
    "category": "constant",
    "text": "Nominal solar photospheric radius. [m]\n\nO. Montenbruck, and E. Gill, Satellite Orbits: Models, Methods and  Applications, 2012.\n\n\n\n\n\n"
},

{
    "location": "modules/constants/#SatelliteDynamics.Constants.P_SUN",
    "page": "Constants",
    "title": "SatelliteDynamics.Constants.P_SUN",
    "category": "constant",
    "text": "Nominal solar radiation pressure at 1 AU. [N/m^2]\n\nO. Montenbruck, and E. Gill, Satellite Orbits: Models, Methods and  Applications, 2012.\n\n\n\n\n\n"
},

{
    "location": "modules/constants/#SatelliteDynamics.Constants.GM_MOON",
    "page": "Constants",
    "title": "SatelliteDynamics.Constants.GM_MOON",
    "category": "constant",
    "text": "Gravitational constant of the Moon. [m^3/s^2]\n\nO. Montenbruck, and E. Gill, Satellite Orbits: Models, Methods and  Applications, 2012.\n\n\n\n\n\n"
},

{
    "location": "modules/constants/#SatelliteDynamics.Constants.GM_MERCURY",
    "page": "Constants",
    "title": "SatelliteDynamics.Constants.GM_MERCURY",
    "category": "constant",
    "text": "Gravitational constant of the Mercury. [m^3/s^2]\n\nO. Montenbruck, and E. Gill, Satellite Orbits: Models, Methods and  Applications, 2012.\n\n\n\n\n\n"
},

{
    "location": "modules/constants/#SatelliteDynamics.Constants.GM_VENUS",
    "page": "Constants",
    "title": "SatelliteDynamics.Constants.GM_VENUS",
    "category": "constant",
    "text": "Gravitational constant of the Venus. [m^3/s^2]\n\nO. Montenbruck, and E. Gill, Satellite Orbits: Models, Methods and  Applications, 2012.\n\n\n\n\n\n"
},

{
    "location": "modules/constants/#SatelliteDynamics.Constants.GM_MARS",
    "page": "Constants",
    "title": "SatelliteDynamics.Constants.GM_MARS",
    "category": "constant",
    "text": "Gravitational constant of the Mars. [m^3/s^2]\n\nO. Montenbruck, and E. Gill, Satellite Orbits: Models, Methods and  Applications, 2012.\n\n\n\n\n\n"
},

{
    "location": "modules/constants/#SatelliteDynamics.Constants.GM_JUPITER",
    "page": "Constants",
    "title": "SatelliteDynamics.Constants.GM_JUPITER",
    "category": "constant",
    "text": "Gravitational constant of the Jupiter. [m^3/s^2]\n\nO. Montenbruck, and E. Gill, Satellite Orbits: Models, Methods and  Applications, 2012.\n\n\n\n\n\n"
},

{
    "location": "modules/constants/#SatelliteDynamics.Constants.GM_SATURN",
    "page": "Constants",
    "title": "SatelliteDynamics.Constants.GM_SATURN",
    "category": "constant",
    "text": "Gravitational constant of the Saturn. [m^3/s^2]\n\nO. Montenbruck, and E. Gill, Satellite Orbits: Models, Methods and  Applications, 2012.\n\n\n\n\n\n"
},

{
    "location": "modules/constants/#SatelliteDynamics.Constants.GM_URANUS",
    "page": "Constants",
    "title": "SatelliteDynamics.Constants.GM_URANUS",
    "category": "constant",
    "text": "Gravitational constant of the Uranus. [m^3/s^2]\n\nO. Montenbruck, and E. Gill, Satellite Orbits: Models, Methods and  Applications, 2012.\n\n\n\n\n\n"
},

{
    "location": "modules/constants/#SatelliteDynamics.Constants.GM_NEPTUNE",
    "page": "Constants",
    "title": "SatelliteDynamics.Constants.GM_NEPTUNE",
    "category": "constant",
    "text": "Gravitational constant of the Neptune. [m^3/s^2]\n\nO. Montenbruck, and E. Gill, Satellite Orbits: Models, Methods and  Applications, 2012.\n\n\n\n\n\n"
},

{
    "location": "modules/constants/#SatelliteDynamics.Constants.GM_PLUTO",
    "page": "Constants",
    "title": "SatelliteDynamics.Constants.GM_PLUTO",
    "category": "constant",
    "text": "Gravitational constant of the Pluto. [m^3/s^2]\n\nO. Montenbruck, and E. Gill, Satellite Orbits: Models, Methods and  Applications, 2012.\n\n\n\n\n\n"
},

{
    "location": "modules/constants/#Constants-1",
    "page": "Constants",
    "title": "Constants",
    "category": "section",
    "text": "The constants submodule constains common mathematical, physical, astronomical constants.AS2RAD\nRAD2AS\nC_LIGHT\nAU\nMJD_ZERO\nMJD2000\nGPS_TAI\nTAI_GPS\nTT_TAI\nTAI_TT\nGPS_TT\nTT_GPS\nGPS_ZERO\nR_EARTH\nWGS84_a\nWGS84_f\nGM_EARTH\ne_EARTH\nJ2_EARTH\nOMEGA_EARTH\nGM_SUN\nR_SUN\nP_SUN\nGM_MOON\nGM_MERCURY\nGM_VENUS\nGM_MARS\nGM_JUPITER\nGM_SATURN\nGM_URANUS\nGM_NEPTUNE\nGM_PLUTO"
},

{
    "location": "modules/universe/#",
    "page": "Universe",
    "title": "Universe",
    "category": "page",
    "text": ""
},

{
    "location": "modules/universe/#SatelliteDynamics.Universe.EarthOrientationData",
    "page": "Universe",
    "title": "SatelliteDynamics.Universe.EarthOrientationData",
    "category": "type",
    "text": "The EarthOrientationData constains a single data member of type  Dict{Int32, Tuple{Float64, Float64, Float64}} that stores the Earth Orientation parameters UT1-UTC, xp, and yp whose units are meters,  radians, and radians, respectively. xp and yp are the x- and  y-components of Earth\'s polar motion. The dictionary key is the Epoch the  parameters are for as a Modified Julian Day at 0h UTC.\n\nArguments:\n\nproduct::Symbol The IERS product type can be :C04_14, :C04_80, or :FINALS_2000\n\n\n\n\n\n"
},

{
    "location": "modules/universe/#SatelliteDynamics.Universe.EOP",
    "page": "Universe",
    "title": "SatelliteDynamics.Universe.EOP",
    "category": "constant",
    "text": "Module-wide global EarthOrientationData object. This data object is used as the default source of Earth Orientation Data by reference system transformations if no explicit EarthOrientationData file is provided to those transformations.\n\nThis value can be overridden in your own code as follows:\n\nSatelliteDynamics.EOP = EarthOrientationData(:EOP_PRODUCT_CHOICE)\n\nThis global variable defaults to use the module\'s internal version of :FINALS_2000  if it is not otherwise set/provided.\n\n\n\n\n\n"
},

{
    "location": "modules/universe/#SatelliteDynamics.Universe.UT1_UTC",
    "page": "Universe",
    "title": "SatelliteDynamics.Universe.UT1_UTC",
    "category": "function",
    "text": "Compute the offset between the UT1 and UTC time systems in seconds. If the EarthOrientationData argument is ommitted the function will use the default module-global value.\n\nArguments:\n\neop::EarthOrientationData EarthOrientationData object to use to compute the offset\nmjd::Real Modified Julian Date in UTC of the Epoch for which the UT1-UTC offset is desired.\ninterp::Bool Whether to linearly interpolate the parameter data to the input MJD.\n\nReturns:\n\nut1_utc::Float UT1 - UTC offset. [s] \n\n\n\n\n\n"
},

{
    "location": "modules/universe/#SatelliteDynamics.Universe.POLE_LOCATOR",
    "page": "Universe",
    "title": "SatelliteDynamics.Universe.POLE_LOCATOR",
    "category": "function",
    "text": "Compute the location of the pole. Returns x- and y- components as a tuple with the units of [radians].  If the EarthOrientationData argument is ommitted the function will use the default module-global value.\n\nArguments:\n\neop::EarthOrientationData EarthOrientationData object to use to compute the offset\nmjd::Real Modified Julian Date in UTC of the Epoch for which the pole locator is desired.\ninterp::Bool Whether to linearly interpolate the parameter data to the input MJD.\n\nReturns:\n\npole_locator::Tuple{ -Float, Float} (x, y) pole location in radians.\n\n\n\n\n\n"
},

{
    "location": "modules/universe/#SatelliteDynamics.Universe.XP",
    "page": "Universe",
    "title": "SatelliteDynamics.Universe.XP",
    "category": "function",
    "text": "Compute the x-component of the pole locator in [radians]. If the first EarthOrientationData argument is ommitted the function will use the default module-global value.\n\nArguments:\n\neop::EarthOrientationData EarthOrientationData object to use to compute the offset\nmjd::Real Modified Julian Date in UTC of the Epoch for which the xp value is desired.\ninterp::Bool Whether to linearly interpolate the parameter data to the input MJD.\n\nReturns:\n\nxp::Float x-component of pole locator in radians.\n\n\n\n\n\n"
},

{
    "location": "modules/universe/#SatelliteDynamics.Universe.YP",
    "page": "Universe",
    "title": "SatelliteDynamics.Universe.YP",
    "category": "function",
    "text": "Compute the y-component of the pole locator in [radians]. If the first EarthOrientationData argument is ommitted the function will use the default module-global value.\n\nArguments:\n\neop::EarthOrientationData EarthOrientationData object to use to compute the offset\nmjd::Real Modified Julian Date in UTC of the Epoch for which the yp value is desired.\ninterp::Bool Whether to linearly interpolate the parameter data to the input MJD.\n\nReturns:\n\nyp::Float y-component of pole locator in radians.\n\n\n\n\n\n"
},

{
    "location": "modules/universe/#SatelliteDynamics.Universe.set_eop",
    "page": "Universe",
    "title": "SatelliteDynamics.Universe.set_eop",
    "category": "function",
    "text": "Set Earth orientation data values for a specific date in the module global EarthOrientationData object.\n\nArguments:\n\nmjd::Real Modified Julian Date in UTC of the Epoch for which the Earth orientation data is aligned to.\nut1_utc::Real Offset between UT1 and UTC in seconds.\nxp::Real x-component of the pole locator in radians.\nyp::Real y-component of the pole locator in radians.\n\n\n\n\n\n"
},

{
    "location": "modules/universe/#SatelliteDynamics.Universe.load_eop",
    "page": "Universe",
    "title": "SatelliteDynamics.Universe.load_eop",
    "category": "function",
    "text": "Load new Earth orientation data into the module global EarthOrientationData object. The product can be one of the symbols: :C04_14, :C04_80, or :FINALS_2000.\n\nArguments:\n\nproduct::Symbol Loads a different set of EarthOrientationData values into the module-wide global EarthOrientationData parameters.\n\n\n\n\n\n"
},

{
    "location": "modules/universe/#SatelliteDynamics.Universe.update_eop",
    "page": "Universe",
    "title": "SatelliteDynamics.Universe.update_eop",
    "category": "function",
    "text": "Download updated Earth orientation datafiles for included products IERS products.\n\nArguments:\n\nproduct::Symbol The IERS product type can be :C04_14, :C04_80, or :FINALS_2000\n\n\n\n\n\n"
},

{
    "location": "modules/universe/#SatelliteDynamics.Universe.GravModel",
    "page": "Universe",
    "title": "SatelliteDynamics.Universe.GravModel",
    "category": "type",
    "text": "GravModel stores a spherical harmonic gravity field in memory. Can store normalized or denomalized coefficients. Package contains EGM2008, GGM01S, and GGM0S gravity models, as well as the default gravity model of EGM2008 truncated to degree and order 90.\n\nAdditional gravity field models can be downloaded from: http://icgem.gfz-potsdam.de/home\n\nArguments:\n\nfilepath::string Path to spherical harmonic gravity model file.\n\n\n\n\n\n"
},

{
    "location": "modules/universe/#SatelliteDynamics.Universe.GRAVITY_MODEL",
    "page": "Universe",
    "title": "SatelliteDynamics.Universe.GRAVITY_MODEL",
    "category": "constant",
    "text": "Module-wide global GravityModel object. This data object is used as the default spherical harmonic gravity field unless one is otherwise provided.\n\nThis value can be overridden in your own code as follows:\n\nSatelliteDynamics.GravityModel = GravityModel(PATH_TO_YOUR_GRAVITY_MODEL)\n\nThis global variable defaults to use the module\'s internal version of the EGM2008 model truncated to order and degree 90, if it is not otherwise set.\n\n\n\n\n\n"
},

{
    "location": "modules/universe/#SatelliteDynamics.Universe.load_gravity_model",
    "page": "Universe",
    "title": "SatelliteDynamics.Universe.load_gravity_model",
    "category": "function",
    "text": "Load new gravity model into module global EarthOrientationData object. The product can be one of the symbols: :EGM2008_20, :EGM2008_90, :GGM01S, :GGM05S, or the filepath to a text-encoded gravity model file.\n\nArguments:\n\ngfc_file::String File path of gravity field model\nproduct_name::Symbol OR a symbol of a known gravity field product. Valid ones are: :EGM2008_20, :EGM2008_90, :GGM01S, :GGM05S\n\n\n\n\n\n"
},

{
    "location": "modules/universe/#Universe-1",
    "page": "Universe",
    "title": "Universe",
    "category": "section",
    "text": "The Universe submodule defines simulation-specific data files which are constants  of most simulations. In particular it provides data structures for storing and accessing Earth orientation parameters and spherical harmonic gravity field  models.The module defines the global variables EOP and GRAVITY_MODEL which are loaded at runtime and .EOP defaults to use the rapid Earth orientation data file finals.all (IAU 2000) distributed by the IERS. The module also supports IERS C04 product files.GRAVITY_MODEL defaults to use the EGM2008 spherical harmonic gravity model,  truncated to order and degree 90.EarthOrientationData\nEOP\nUT1_UTC\nPOLE_LOCATOR\nXP\nYP\nset_eop\nload_eop\nupdate_eop\nGravModel\nGRAVITY_MODEL\nload_gravity_model"
},

{
    "location": "modules/time/#",
    "page": "Time",
    "title": "Time",
    "category": "page",
    "text": ""
},

{
    "location": "modules/time/#SatelliteDynamics.Time.caldate_to_mjd",
    "page": "Time",
    "title": "SatelliteDynamics.Time.caldate_to_mjd",
    "category": "function",
    "text": "Convert a Gregorian calendar date to the equivalent Modified Julian Date representation of that time instant.\n\nAguments:\n\nyear::Int Year\nyear::Int Month\nyear::Int Day\nhour::Int Hour\nminute::Int Minute \nsecond::Real Seconds\nnanoseconds::Real Microseconds\n\nReturns:\n\nmjd::Float64 Modified Julian Date of Epoch\n\n\n\n\n\n"
},

{
    "location": "modules/time/#SatelliteDynamics.Time.mjd_to_caldate",
    "page": "Time",
    "title": "SatelliteDynamics.Time.mjd_to_caldate",
    "category": "function",
    "text": "Convert a Modified Julian Date to the equivalent Gregorian calendar date representation of the same instant in time.\n\nAguments:\n\nmjd::Real: Modified Julian Date of Epoch\n\nReturns:\n\nyear::Int32: Year\nyear::Int32: Month\nyear::Int32: Day\nhour::Int32: Hour\nminute::Int32: Minute \nsecond::Float64: Seconds\nnanoseconds::Float64: nanosecondss\n\n\n\n\n\n"
},

{
    "location": "modules/time/#SatelliteDynamics.Time.caldate_to_jd",
    "page": "Time",
    "title": "SatelliteDynamics.Time.caldate_to_jd",
    "category": "function",
    "text": "Convert a Gregorian calendar date to the equivalent Julian Date representation of that time instant.\n\nAguments:\n\nyear::Int: Year\nyear::Int: Month\nyear::Int: Day\nhour::Int: Hour\nminute::Int: Minute \nsecond::Real: Seconds\nnanoseconds::Real: nanosecondss\n\nReturns:\n\nmjd::Float64: Julian Date of Epoch\n\n\n\n\n\n"
},

{
    "location": "modules/time/#SatelliteDynamics.Time.jd_to_caldate",
    "page": "Time",
    "title": "SatelliteDynamics.Time.jd_to_caldate",
    "category": "function",
    "text": "Convert a Julian Date to the equivalent Gregorian calendar date representation of the same instant in time.\n\nAguments:\n\nmjd::Real: Julian Date of Epoch\n\nReturns:\n\nyear::Int32: Year\nyear::Int32: Month\nyear::Int32: Day\nhour::Int32: Hour\nminute::Int32: Minute \nsecond::Float64: Seconds\nmicrosecond::Float64: Microseconds\n\n\n\n\n\n"
},

{
    "location": "modules/time/#SatelliteDynamics.Time.elapsed_from_epoch",
    "page": "Time",
    "title": "SatelliteDynamics.Time.elapsed_from_epoch",
    "category": "function",
    "text": "Compute the number of elapsed seconds since a given Epoch from the day number. Can be used to compute the elapsed time since a given Julian or Modified Julian Date.\n\nArguments:\n\nday_number::Real: Day number, can contain fractional days. Asummes all days are a uniform 86400.0 seconds in length.\nday_epoch::Real: Modified Julian Date of Epoch\n\nReturns:\n\nt::Float: Number of elapsed seconds between the Provided Modified   Julian date and the epoch.\n\n\n\n\n\n"
},

{
    "location": "modules/time/#SatelliteDynamics.Time.days_from_elapsed",
    "page": "Time",
    "title": "SatelliteDynamics.Time.days_from_elapsed",
    "category": "function",
    "text": "Computes the day number in a given time scale given the elapsed time since epoch and the epoch itself.\n\nAssumes all days are counted using a uniform 86400.0 seconds over the time span.\n\nArguments:\n\nt::Real: Elapsed seconds since the day_epoch.\nday_epoch::Real: Day number of the epoch. Common values are SatelliteDynamics.Constants.MJD_ZERO (to get the Julian Day number) or SatelliteDynamics.Constants.MJD2000 (to get Modified Julian Days if reckoning time from January 1, 2000 0H)\n\nReturns:\n\ndays::Float: Number of elapsed days in the time scale.\n\n\n\n\n\n"
},

{
    "location": "modules/time/#Time-1",
    "page": "Time",
    "title": "Time",
    "category": "section",
    "text": "The Time submodule contains common time transformations such as converting between different date representations or converting a specific instant in time between different time systems.The module also defines the Epoch class which pMost of the transformations are make backend calls to the SOFA C-library functions provide the package SOFA.jlcaldate_to_mjd\nmjd_to_caldate\ncaldate_to_jd\njd_to_caldate\nelapsed_from_epoch\ndays_from_elapsed"
},

{
    "location": "modules/function_index/#",
    "page": "Function Index",
    "title": "Function Index",
    "category": "page",
    "text": ""
},

{
    "location": "modules/function_index/#Function-Library-1",
    "page": "Function Index",
    "title": "Function Library",
    "category": "section",
    "text": ""
},

{
    "location": "modules/function_index/#Constants-1",
    "page": "Function Index",
    "title": "Constants",
    "category": "section",
    "text": "Modules = [SatelliteDynamics.Constants]"
},

{
    "location": "modules/function_index/#Universe-1",
    "page": "Function Index",
    "title": "Universe",
    "category": "section",
    "text": "Modules = [SatelliteDynamics.Universe]"
},

{
    "location": "modules/function_index/#Time-1",
    "page": "Function Index",
    "title": "Time",
    "category": "section",
    "text": "Modules = [SatelliteDynamics.Time]"
},

]}
