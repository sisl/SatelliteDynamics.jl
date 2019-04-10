let 
    tol = 1.0e-7

    # Test known position conversions
    geoc1 = [0.0, 0.0, 0.0]
    ecef1 = sGEOCtoECEF(geoc1)

    @test isapprox(ecef1[1], WGS84_a, atol=tol)
    @test isapprox(ecef1[2], 0, atol=tol)
    @test isapprox(ecef1[3], 0, atol=tol)

    geoc2 = [90.0, 0.0, 0.0]
    ecef2 = sGEOCtoECEF(geoc2, use_degrees=true)

    @test isapprox(ecef2[1], 0, atol=tol)
    @test isapprox(ecef2[2], WGS84_a, atol=tol)
    @test isapprox(ecef2[3], 0, atol=tol)

    geoc3 = [0.0, 90.0, 0.0]
    ecef3 = sGEOCtoECEF(geoc3, use_degrees=true)

    @test isapprox(ecef3[1], 0, atol=tol)
    @test isapprox(ecef3[2], 0, atol=tol)
    @test isapprox(ecef3[3], WGS84_a, atol=tol)

    # Test two-input format 
    geoc = [0.0, 0.0]
    ecef = sGEOCtoECEF(geoc)

    @test isapprox(ecef[1], WGS84_a, atol=tol)
    @test isapprox(ecef[2], 0, atol=tol)
    @test isapprox(ecef[3], 0, atol=tol)

    geoc = [90.0, 0.0]
    ecef = sGEOCtoECEF(geoc, use_degrees=true)

    @test isapprox(ecef[1], 0, atol=tol)
    @test isapprox(ecef[2], WGS84_a, atol=tol)
    @test isapprox(ecef[3], 0, atol=tol)

    geoc = [0.0, 90.0]
    ecef = sGEOCtoECEF(geoc, use_degrees=true)

    @test isapprox(ecef[1], 0, atol=tol)
    @test isapprox(ecef[2], 0, atol=tol)
    @test isapprox(ecef[3], WGS84_a, atol=tol)

    # Test circularity
    geoc4 = sECEFtoGEOC(ecef1, use_degrees=true)
    geoc5 = sECEFtoGEOC(ecef2, use_degrees=true)
    geoc6 = sECEFtoGEOC(ecef3, use_degrees=true)

    @test isapprox(geoc4[1], geoc1[1], atol=tol)
    @test isapprox(geoc4[2], geoc1[2], atol=tol)
    @test isapprox(geoc4[3], geoc1[3], atol=tol)

    @test isapprox(geoc5[1], geoc2[1], atol=tol)
    @test isapprox(geoc5[2], geoc2[2], atol=tol)
    @test isapprox(geoc5[3], geoc2[3], atol=tol)

    @test isapprox(geoc6[1], geoc3[1], atol=tol)
    @test isapprox(geoc6[2], geoc3[2], atol=tol)
    @test isapprox(geoc6[3], geoc3[3], atol=tol)

    # Random point circularity
    geoc  = [77.875000,    20.975200,     0.000000]
    ecef  = sGEOCtoECEF(geoc, use_degrees=true)
    geocc = sECEFtoGEOC(ecef, use_degrees=true)
    @test isapprox(geoc[1], geocc[1], atol=tol)
    @test isapprox(geoc[2], geocc[2], atol=tol)
    @test isapprox(geoc[3], geocc[3], atol=tol)

    # Test Error Condition
    @test_throws ArgumentError sGEOCtoECEF([0.0,  90.1], use_degrees=true)
    @test_throws ArgumentError sGEOCtoECEF([0.0, -90.1], use_degrees=true)    
end

let 
    tol = 1.0e-7

    # Test known position conversions
    geod1 = [0, 0, 0]
    ecef1 = sGEODtoECEF(geod1)

    @test isapprox(ecef1[1], WGS84_a, atol=tol)
    @test isapprox(ecef1[2], 0, atol=tol)
    @test isapprox(ecef1[3], 0, atol=tol)

    geod2 = [90.0, 0.0, 0.0]
    ecef2 = sGEODtoECEF(geod2, use_degrees=true)

    @test isapprox(ecef2[1], 0, atol=tol)
    @test isapprox(ecef2[2], WGS84_a, atol=tol)
    @test isapprox(ecef2[3], 0, atol=tol)

    geod3 = [0, 90.0, 0]
    ecef3 = sGEODtoECEF(geod3, use_degrees=true)

    @test isapprox(ecef3[1], 0, atol=tol)
    @test isapprox(ecef3[2], 0, atol=tol)
    @test isapprox(ecef3[3], WGS84_a*(1.0-WGS84_f), atol=tol)

    # Test two input format
    geod = [0.0, 0.0]
    ecef = sGEODtoECEF(geod)

    @test isapprox(ecef[1], WGS84_a, atol=tol)
    @test isapprox(ecef[2], 0, atol=tol)
    @test isapprox(ecef[3], 0, atol=tol)

    geod = [90.0, 0.0]
    ecef = sGEODtoECEF(geod, use_degrees=true)

    @test isapprox(ecef[1], 0, atol=tol)
    @test isapprox(ecef[2], WGS84_a, atol=tol)
    @test isapprox(ecef[3], 0, atol=tol)

    geod = [0.0, 90.0]
    ecef = sGEODtoECEF(geod, use_degrees=true)

    @test isapprox(ecef[1], 0, atol=tol)
    @test isapprox(ecef[2], 0, atol=tol)
    @test isapprox(ecef[3], WGS84_a*(1.0-WGS84_f), atol=tol)

    # Test circularity
    geod4 = sECEFtoGEOD(ecef1, use_degrees=true)
    geod5 = sECEFtoGEOD(ecef2, use_degrees=true)
    geod6 = sECEFtoGEOD(ecef3, use_degrees=true)

    @test isapprox(geod4[1], geod1[1], atol=tol)
    @test isapprox(geod4[2], geod1[2], atol=tol)
    @test isapprox(geod4[3], geod1[3], atol=tol)

    @test isapprox(geod5[1], geod2[1], atol=tol)
    @test isapprox(geod5[2], geod2[2], atol=tol)
    @test isapprox(geod5[3], geod2[3], atol=tol)

    @test isapprox(geod6[1], geod3[1], atol=tol)
    @test isapprox(geod6[2], geod3[2], atol=tol)
    @test isapprox(geod6[3], geod3[3], atol=tol)

    geod  = [77.875000,    20.975200,     0.000000]
    ecef  = sGEODtoECEF(geod, use_degrees=true)
    geodc = sECEFtoGEOD(ecef, use_degrees=true)
    @test isapprox(geod[1], geodc[1], atol=tol)
    @test isapprox(geod[2], geodc[2], atol=tol)
    @test isapprox(geod[3], geodc[3], atol=tol)

    # Test Error Condition
    @test_throws ArgumentError sGEODtoECEF([0.0,  90.1], use_degrees=true)
    @test_throws ArgumentError sGEODtoECEF([0.0, -90.1], use_degrees=true)    
end

let
    tol = 1.0e-8

    station_ecef = [0, R_EARTH, 0]

    R_ecef_enz = rECEFtoENZ(station_ecef, conversion="geocentric")
    R_enz_ecef = rENZtoECEF(station_ecef, conversion="geocentric")

    @test R_ecef_enz == R_enz_ecef'

    # State conversion
    epc  = Epoch(2018,1,1,12,0,0)
    oe   = [R_EARTH + 500e3, 1e-3, 97.8, 75, 25, 45]
    ecef = sECItoECEF(epc, sOSCtoCART(oe, use_degrees=true))

    station_ecef = sGEODtoECEF([-122.4056, 37.7716, 0.0], use_degrees=true)

    enz   = sECEFtoENZ(station_ecef, ecef)
    ecef2 = sENZtoECEF(station_ecef, enz)

    @test isapprox(ecef[1], ecef2[1], atol=tol)
    @test isapprox(ecef[2], ecef2[2], atol=tol)
    @test isapprox(ecef[3], ecef2[3], atol=tol)
    @test isapprox(ecef[4], ecef2[4], atol=tol)
    @test isapprox(ecef[5], ecef2[5], atol=tol)
    @test isapprox(ecef[6], ecef2[6], atol=tol)

    ecef         = sGEODtoECEF([-122.4, 37.78, 200.0],    use_degrees=true)
    station_ecef = sGEODtoECEF([-122.4056, 37.7716, 0.0], use_degrees=true)

    enz   = sECEFtoENZ(station_ecef, ecef, conversion="geocentric")
    ecef2 = sENZtoECEF(station_ecef, enz, conversion="geocentric")

    @test isapprox(ecef[1], ecef2[1], atol=tol)
    @test isapprox(ecef[2], ecef2[2], atol=tol)
    @test isapprox(ecef[3], ecef2[3], atol=tol)

    # Test ENZ Error Conditions
    @test_throws ArgumentError rECEFtoENZ([R_EARTH, 0.0])
    @test_throws ArgumentError rECEFtoENZ([R_EARTH, 0.0, 0.0], conversion="unknown")
    @test_throws ArgumentError rENZtoECEF([R_EARTH, 0.0])
    @test_throws ArgumentError sECEFtoENZ([R_EARTH, 0.0], [R_EARTH + 100.0, 0.0, 0.0])
    @test_throws ArgumentError sECEFtoENZ([R_EARTH, 0.0, 0.0], [R_EARTH + 100.0, 0.0])
    @test_throws ArgumentError sENZtoECEF([R_EARTH, 0.0], [0.0, 0.0, 0.0])
    @test_throws ArgumentError sENZtoECEF([R_EARTH, 0.0, 0.0], [0.0, 0.0])

    # Test length of return is 3
    enz = sECEFtoENZ(station_ecef, ecef[1:3], conversion="geocentric")
    @test length(enz) == 3
end

let
    tol = 1.0e-8
    station_ecef = [0, R_EARTH, 0]

    R_ecef_sez = rECEFtoSEZ(station_ecef, conversion="geocentric")
    R_sez_ecef = rSEZtoECEF(station_ecef, conversion="geocentric")

    @test isapprox(R_ecef_sez[1, 1], R_sez_ecef[1, 1], atol=tol)
    @test isapprox(R_ecef_sez[2, 1], R_sez_ecef[1, 2], atol=tol)
    @test isapprox(R_ecef_sez[3, 1], R_sez_ecef[1, 3], atol=tol)

    @test isapprox(R_ecef_sez[1, 2], R_sez_ecef[2, 1], atol=tol)
    @test isapprox(R_ecef_sez[2, 2], R_sez_ecef[2, 2], atol=tol)
    @test isapprox(R_ecef_sez[3, 2], R_sez_ecef[2, 3], atol=tol)

    @test isapprox(R_ecef_sez[1, 3], R_sez_ecef[3, 1], atol=tol)
    @test isapprox(R_ecef_sez[2, 3], R_sez_ecef[3, 2], atol=tol)
    @test isapprox(R_ecef_sez[3, 3], R_sez_ecef[3, 3], atol=tol)

    # State conversion
    epc  = Epoch(2018,1,1,12,0,0)
    oe   = [R_EARTH + 500e3, 1e-3, 97.8, 75, 25, 45]
    ecef = sECItoECEF(epc, sOSCtoCART(oe, use_degrees=true))

    station_ecef = sGEODtoECEF([-122.4056, 37.7716, 0.0], use_degrees=true)

    sez   = sECEFtoSEZ(station_ecef, ecef)
    ecef2 = sSEZtoECEF(station_ecef, sez)

    @test isapprox(ecef[1], ecef2[1], atol=tol)
    @test isapprox(ecef[2], ecef2[2], atol=tol)
    @test isapprox(ecef[3], ecef2[3], atol=tol)
    @test isapprox(ecef[4], ecef2[4], atol=tol)
    @test isapprox(ecef[5], ecef2[5], atol=tol)
    @test isapprox(ecef[6], ecef2[6], atol=tol)

    ecef         = sGEODtoECEF([-122.4, 37.78, 200.0],    use_degrees=true)
    station_ecef = sGEODtoECEF([-122.4056, 37.7716, 0.0], use_degrees=true)

    sez   = sECEFtoSEZ(station_ecef, ecef, conversion="geocentric")
    ecef2 = sSEZtoECEF(station_ecef, sez, conversion="geocentric")

    @test isapprox(ecef[1], ecef2[1], atol=tol)
    @test isapprox(ecef[2], ecef2[2], atol=tol)
    @test isapprox(ecef[3], ecef2[3], atol=tol)

    # Test SEZ Error Conditions
    @test_throws ArgumentError rECEFtoSEZ([R_EARTH, 0.0])
    @test_throws ArgumentError rECEFtoSEZ([R_EARTH, 0.0, 0.0], conversion="unknown")
    @test_throws ArgumentError rSEZtoECEF([R_EARTH, 0.0])
    @test_throws ArgumentError sECEFtoSEZ([R_EARTH, 0.0], [R_EARTH + 100.0, 0.0, 0.0])
    @test_throws ArgumentError sECEFtoSEZ([R_EARTH, 0.0, 0.0], [R_EARTH + 100.0, 0.0])
    @test_throws ArgumentError sSEZtoECEF([R_EARTH, 0.0], [0.0, 0.0, 0.0])
    @test_throws ArgumentError sSEZtoECEF([R_EARTH, 0.0, 0.0], [0.0, 0.0])

    # Test length of return is 3
    sez = sECEFtoSEZ(station_ecef, ecef[1:3], conversion="geocentric")
    @test length(sez) == 3
end

let
    # Test taken from Montenbruck and Gill Exercise 2.4
    # It mixes geodetic and geocentric coordinations in a strange way, but the
    # mixing is retained here for consistenty with the source material test
    epc = Epoch(1997, 1, 1, 0, 0, 0, tsys=:UTC)
    oe  = [6378.137e3 + 960e3, 0, 97, 130.7, 0, 0]
    dt  = 15*60

    # Get Satellite position at 15 minutes
    n = mean_motion(oe[1], use_degrees=true)
    oe[6] += n*dt

    sat_eci = sOSCtoCART(oe, use_degrees=true)

    # Low precision ECEF transform
    d = (dt/86400.0 + mjd(epc) - 51544.5)
    # O = 280.4606 + 360.9856473*d
    O = 1.82289510683
    sat_ecef = Rz(0, use_degrees=false) * sat_eci[1:3]

    # Station coordinates
    station_ecef = sGEODtoECEF([48.0, 11.0, 0.0], use_degrees=true)


    # Compute enz and sez state
    enz   = sECEFtoENZ(station_ecef, sat_ecef, conversion="geocentric")
    sez   = sECEFtoSEZ(station_ecef, sat_ecef, conversion="geocentric")


    # Compute azimuth and elevation from topocentric coordinates
    azel_enz = sENZtoAZEL(enz, use_degrees=true)
    azel_sez = sSEZtoAZEL(sez, use_degrees=true)

    @test azel_enz[1] == azel_sez[1]
    @test azel_enz[2] == azel_sez[2]
    @test azel_enz[3] == azel_sez[3]
end

let
    # State conversion
    epc  = Epoch(2018,1,1,12,0,0)
    oe   = [R_EARTH + 500e3, 1e-3, 97.8, 75, 25, 45]
    ecef = sECItoECEF(epc, sOSCtoCART(oe, use_degrees=true))

    station_ecef = sGEODtoECEF([-122.4056, 37.7716, 0.0], use_degrees=true)

    enz = sECEFtoENZ(station_ecef, ecef)
    sez = sECEFtoSEZ(station_ecef, ecef)

    # Compute azimuth and elevation from topocentric coordinates
    azel_enz = sENZtoAZEL(enz, use_degrees=true)
    azel_sez = sSEZtoAZEL(sez, use_degrees=true)

    @test azel_enz[1] == azel_sez[1]
    @test azel_enz[2] == azel_sez[2]
    @test azel_enz[3] == azel_sez[3]
    @test azel_enz[4] == azel_sez[4]
    @test azel_enz[5] == azel_sez[5]
    @test azel_enz[6] == azel_sez[6]
end

let 
    # Test Error Conditions
    enz = [0.0, 0.0, 100, 90.0, 0.0, 0.0]

    # Non-standard input length
    @test_throws ArgumentError sENZtoAZEL(enz[1:2])

    # Cant resolve azimuth without range information
    azel = sENZtoAZEL(enz[1:3])
    @test azel[1] == 0.0

    # Test ability to resolve azimuth ambiguity
    azel = sENZtoAZEL(enz)
    @test azel[1] != 0
    @test azel[2] != 0
    @test azel[3] != 0
end

let 
    # Test Error Conditions
    sez = [0.0, 0.0, 100, 90.0, 0.0, 0.0]

    # Non-standard input length
    @test_throws ArgumentError sSEZtoAZEL(sez[1:2])

    # Cant resolve azimuth without range information
    azel = sSEZtoAZEL(sez[1:3])
    @test azel[1] == 0.0

    # Test ability to resolve azimuth ambiguity
    azel = sSEZtoAZEL(sez)
    @test azel[1] != 0
    @test azel[2] != 0
    @test azel[3] != 0
end