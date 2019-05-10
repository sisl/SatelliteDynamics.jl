let
    using SatelliteDynamics.EarthEnvironment.NRLMSISE00: NRLMSISE_Input, NRLMSISE_Output, NRLMSISE_Flags, gtd7!
    # Test core NRLMSISE - Test cases match Daniel Brodo's C implementation
    input  = NRLMSISE_Input[NRLMSISE_Input() for i in 1:17]
    output = NRLMSISE_Output[NRLMSISE_Output() for i in 1:17]
    flags  = NRLMSISE_Flags()

    # Set input values
    aph = ones(Float64, 7).*100.0

    flags.switches[1] = 0
    for i in 2:24
        flags.switches[i] = 1
    end

    for i in 1:17
        input[i].doy   = 172
        input[i].year  = 0
        input[i].sec   = 29000
        input[i].alt   = 400
        input[i].g_lat = 60
        input[i].g_lon = -70
        input[i].lst   = 16
        input[i].f107A = 150
        input[i].f107  = 150
        input[i].ap    = 4
    end

    input[2].doy       = 81
    input[3].sec       = 75000
    input[3].alt       = 1000
    input[4].alt       = 100
    input[11].alt      = 0
    input[12].alt      = 10
    input[13].alt      = 30
    input[14].alt      = 50
    input[15].alt      = 70
    input[17].alt      = 100
    input[5].g_lat     = 0
    input[6].g_lon     = 0
    input[7].lst       = 4
    input[8].f107A     = 70
    input[9].f107      = 180
    input[10].ap       = 40
    input[16].ap_array = aph
    input[17].ap_array = aph

    # Run model
    gtd7!(input[1], flags, output[1])

    # Evalute 1 to 15
    for i in 1:15
        gtd7!(input[i], flags, output[i])
    end

    # Evaluate 16 and 17
    for i in 16:17
        gtd7!(input[i], flags, output[i])
    end

    # Test Case 1
    @test isapprox(output[1].t[1], 1250.54, atol=1.0e-2)
    @test isapprox(output[1].t[2], 1241.42, atol=1.0e-2)
    @test isapprox(output[1].d[1], 6.6651769e+05, rtol=1.0e-6)
    @test isapprox(output[1].d[2], 1.1388056e+08, rtol=1.0e-6)
    @test isapprox(output[1].d[3], 1.9982109e+07, rtol=1.0e-6)
    @test isapprox(output[1].d[4], 4.0227636e+05, rtol=1.0e-6)
    @test isapprox(output[1].d[5], 3.5574650e+03, rtol=1.0e-6)
    @test isapprox(output[1].d[7], 3.4753124e+04, rtol=1.0e-6)
    @test isapprox(output[1].d[8], 4.0959133e+06, rtol=1.0e-6)
    @test isapprox(output[1].d[9], 2.6672732e+04, rtol=1.0e-6)
    @test isapprox(output[1].d[6], 4.0747135e-15, rtol=1.0e-6)

    # Test Case 2
    @test isapprox(output[2].t[1], 1166.75, atol=1.0e-2)
    @test isapprox(output[2].t[2], 1161.71, atol=1.0e-2)
    @test isapprox(output[2].d[1], 3.4072932e+06, rtol=1.0e-6)
    @test isapprox(output[2].d[2], 1.5863334e+08, rtol=1.0e-6)
    @test isapprox(output[2].d[3], 1.3911174e+07, rtol=1.0e-6)
    @test isapprox(output[2].d[4], 3.2625595e+05, rtol=1.0e-6)
    @test isapprox(output[2].d[5], 1.5596182e+03, rtol=1.0e-6)
    @test isapprox(output[2].d[7], 4.8542085e+04, rtol=1.0e-6)
    @test isapprox(output[2].d[8], 4.3809667e+06, rtol=1.0e-6)
    @test isapprox(output[2].d[9], 6.9566820e+03, rtol=1.0e-6)
    @test isapprox(output[2].d[6], 5.0018457e-15, rtol=1.0e-6)

    # Test Case 3
    @test isapprox(output[3].t[1], 1239.89, atol=1.0e-2)
    @test isapprox(output[3].t[2], 1239.89, atol=1.0e-2)
    @test isapprox(output[3].d[1], 1.1237672e+05, rtol=1.0e-6)
    @test isapprox(output[3].d[2], 6.9341301e+04, rtol=1.0e-6)
    @test isapprox(output[3].d[3], 4.2471052e+01, rtol=1.0e-6)
    @test isapprox(output[3].d[4], 1.3227501e-01, rtol=1.0e-6)
    @test isapprox(output[3].d[5], 2.6188484e-05, rtol=1.0e-6)
    @test isapprox(output[3].d[7], 2.0167499e+04, rtol=1.0e-6)
    @test isapprox(output[3].d[8], 5.7412559e+03, rtol=1.0e-6)
    @test isapprox(output[3].d[9], 2.3743942e+04, rtol=1.0e-6)
    @test isapprox(output[3].d[6], 2.7567723e-18, rtol=1.0e-6)

    # Test Case 4
    @test isapprox(output[4].t[1], 1027.32, atol=1.0e-2)
    @test isapprox(output[4].t[2], 206.89, atol=1.0e-2)
    @test isapprox(output[4].d[1], 5.4115544e+07, rtol=1.0e-6)
    @test isapprox(output[4].d[2], 1.9188934e+11, rtol=1.0e-6)
    @test isapprox(output[4].d[3], 6.1158256e+12, rtol=1.0e-6)
    @test isapprox(output[4].d[4], 1.2252011e+12, rtol=1.0e-6)
    @test isapprox(output[4].d[5], 6.0232120e+10, rtol=1.0e-6)
    @test isapprox(output[4].d[7], 1.0598797e+07, rtol=1.0e-6)
    @test isapprox(output[4].d[8], 2.6157367e+05, rtol=1.0e-6)
    @test isapprox(output[4].d[9], 2.8198794e-42, rtol=1.0e-6)
    @test isapprox(output[4].d[6], 3.5844263e-10, rtol=1.0e-6)

    # Test Case 5
    @test isapprox(output[5].t[1], 1212.40, atol=1.0e-2)
    @test isapprox(output[5].t[2], 1208.14, atol=1.0e-2)
    @test isapprox(output[5].d[1], 1.8511225e+06, rtol=1.0e-6)
    @test isapprox(output[5].d[2], 1.4765548e+08, rtol=1.0e-6)
    @test isapprox(output[5].d[3], 1.5793562e+07, rtol=1.0e-6)
    @test isapprox(output[5].d[4], 2.6337950e+05, rtol=1.0e-6)
    @test isapprox(output[5].d[5], 1.5887814e+03, rtol=1.0e-6)
    @test isapprox(output[5].d[7], 5.8161668e+04, rtol=1.0e-6)
    @test isapprox(output[5].d[8], 5.4789845e+06, rtol=1.0e-6)
    @test isapprox(output[5].d[9], 1.2644459e+03, rtol=1.0e-6)
    @test isapprox(output[5].d[6], 4.8096302e-15, rtol=1.0e-6)

    # Test Case 6
    @test isapprox(output[6].t[1], 1220.15, atol=1.0e-2)
    @test isapprox(output[6].t[2], 1212.71, atol=1.0e-2)
    @test isapprox(output[6].d[1], 8.6730952e+05, rtol=1.0e-6)
    @test isapprox(output[6].d[2], 1.2788618e+08, rtol=1.0e-6)
    @test isapprox(output[6].d[3], 1.8225766e+07, rtol=1.0e-6)
    @test isapprox(output[6].d[4], 2.9222142e+05, rtol=1.0e-6)
    @test isapprox(output[6].d[5], 2.4029624e+03, rtol=1.0e-6)
    @test isapprox(output[6].d[7], 3.6863892e+04, rtol=1.0e-6)
    @test isapprox(output[6].d[8], 3.8972755e+06, rtol=1.0e-6)
    @test isapprox(output[6].d[9], 2.6672732e+04, rtol=1.0e-6)
    @test isapprox(output[6].d[6], 4.3558656e-15, rtol=1.0e-6)

    # Test Case 7
    @test isapprox(output[7].t[1], 1116.39, atol=1.0e-2)
    @test isapprox(output[7].t[2], 1113.00, atol=1.0e-2)
    @test isapprox(output[7].d[1], 5.7762512e+05, rtol=1.0e-6)
    @test isapprox(output[7].d[2], 6.9791387e+07, rtol=1.0e-6)
    @test isapprox(output[7].d[3], 1.2368136e+07, rtol=1.0e-6)
    @test isapprox(output[7].d[4], 2.4928677e+05, rtol=1.0e-6)
    @test isapprox(output[7].d[5], 1.4057387e+03, rtol=1.0e-6)
    @test isapprox(output[7].d[7], 5.2919856e+04, rtol=1.0e-6)
    @test isapprox(output[7].d[8], 1.0698141e+06, rtol=1.0e-6)
    @test isapprox(output[7].d[9], 2.6672732e+04, rtol=1.0e-6)
    @test isapprox(output[7].d[6], 2.4706514e-15, rtol=1.0e-6)

    # Test Case 8
    @test isapprox(output[8].t[1], 1031.25, atol=1.0e-2)
    @test isapprox(output[8].t[2], 1024.85, atol=1.0e-2)
    @test isapprox(output[8].d[1], 3.7403041e+05, rtol=1.0e-6)
    @test isapprox(output[8].d[2], 4.7827201e+07, rtol=1.0e-6)
    @test isapprox(output[8].d[3], 5.2403800e+06, rtol=1.0e-6)
    @test isapprox(output[8].d[4], 1.7598746e+05, rtol=1.0e-6)
    @test isapprox(output[8].d[5], 5.5016488e+02, rtol=1.0e-6)
    @test isapprox(output[8].d[7], 8.8967757e+04, rtol=1.0e-6)
    @test isapprox(output[8].d[8], 1.9797408e+06, rtol=1.0e-6)
    @test isapprox(output[8].d[9], 9.1218149e+03, rtol=1.0e-6)
    @test isapprox(output[8].d[6], 1.5718887e-15, rtol=1.0e-6)

    # Test Case 9
    @test isapprox(output[9].t[1], 1306.05, atol=1.0e-2)
    @test isapprox(output[9].t[2], 1293.37, atol=1.0e-2)
    @test isapprox(output[9].d[1], 6.7483388e+05, rtol=1.0e-6)
    @test isapprox(output[9].d[2], 1.2453153e+08, rtol=1.0e-6)
    @test isapprox(output[9].d[3], 2.3690095e+07, rtol=1.0e-6)
    @test isapprox(output[9].d[4], 4.9115832e+05, rtol=1.0e-6)
    @test isapprox(output[9].d[5], 4.5787811e+03, rtol=1.0e-6)
    @test isapprox(output[9].d[7], 3.2445948e+04, rtol=1.0e-6)
    @test isapprox(output[9].d[8], 5.3708331e+06, rtol=1.0e-6)
    @test isapprox(output[9].d[9], 2.6672732e+04, rtol=1.0e-6)
    @test isapprox(output[9].d[6], 4.5644202e-15, rtol=1.0e-6)

    # Test Case 10
    @test isapprox(output[10].t[1], 1361.87, atol=1.0e-2)
    @test isapprox(output[10].t[2], 1347.39, atol=1.0e-2)
    @test isapprox(output[10].d[1], 5.5286008e+05, rtol=1.0e-6)
    @test isapprox(output[10].d[2], 1.1980413e+08, rtol=1.0e-6)
    @test isapprox(output[10].d[3], 3.4957978e+07, rtol=1.0e-6)
    @test isapprox(output[10].d[4], 9.3396184e+05, rtol=1.0e-6)
    @test isapprox(output[10].d[5], 1.0962548e+04, rtol=1.0e-6)
    @test isapprox(output[10].d[7], 2.6864279e+04, rtol=1.0e-6)
    @test isapprox(output[10].d[8], 4.8899742e+06, rtol=1.0e-6)
    @test isapprox(output[10].d[9], 2.8054448e+04, rtol=1.0e-6)
    @test isapprox(output[10].d[6], 4.9745431e-15, rtol=1.0e-6)

    # Test Case 11
    @test isapprox(output[11].t[1], 1027.32, atol=1.0e-2)
    @test isapprox(output[11].t[2], 281.46, atol=1.0e-2)
    @test isapprox(output[11].d[1], 1.3754876e+14, rtol=1.0e-6)
    @test isapprox(output[11].d[2], 0.0000000e+00, rtol=1.0e-6)
    @test isapprox(output[11].d[3], 2.0496870e+19, rtol=1.0e-6)
    @test isapprox(output[11].d[4], 5.4986954e+18, rtol=1.0e-6)
    @test isapprox(output[11].d[5], 2.4517332e+17, rtol=1.0e-6)
    @test isapprox(output[11].d[7], 0.0000000e+00, rtol=1.0e-6)
    @test isapprox(output[11].d[8], 0.0000000e+00, rtol=1.0e-6)
    @test isapprox(output[11].d[9], 0.0000000e+00, rtol=1.0e-6)
    @test isapprox(output[11].d[6], 1.2610657e-03, rtol=1.0e-6)

    # Test Case 12
    @test isapprox(output[12].t[1], 1027.32, atol=1.0e-2)
    @test isapprox(output[12].t[2], 227.42, atol=1.0e-2)
    @test isapprox(output[12].d[1], 4.4274426e+13, rtol=1.0e-6)
    @test isapprox(output[12].d[2], 0.0000000e+00, rtol=1.0e-6)
    @test isapprox(output[12].d[3], 6.5975672e+18, rtol=1.0e-6)
    @test isapprox(output[12].d[4], 1.7699293e+18, rtol=1.0e-6)
    @test isapprox(output[12].d[5], 7.8916800e+16, rtol=1.0e-6)
    @test isapprox(output[12].d[7], 0.0000000e+00, rtol=1.0e-6)
    @test isapprox(output[12].d[8], 0.0000000e+00, rtol=1.0e-6)
    @test isapprox(output[12].d[9], 0.0000000e+00, rtol=1.0e-6)
    @test isapprox(output[12].d[6], 4.0591394e-04, rtol=1.0e-6)

    # Test Case 13
    @test isapprox(output[13].t[1], 1027.32, atol=1.0e-2)
    @test isapprox(output[13].t[2], 237.44, atol=1.0e-2)
    @test isapprox(output[13].d[1], 2.1278288e+12, rtol=1.0e-6)
    @test isapprox(output[13].d[2], 0.0000000e+00, rtol=1.0e-6)
    @test isapprox(output[13].d[3], 3.1707906e+17, rtol=1.0e-6)
    @test isapprox(output[13].d[4], 8.5062798e+16, rtol=1.0e-6)
    @test isapprox(output[13].d[5], 3.7927411e+15, rtol=1.0e-6)
    @test isapprox(output[13].d[7], 0.0000000e+00, rtol=1.0e-6)
    @test isapprox(output[13].d[8], 0.0000000e+00, rtol=1.0e-6)
    @test isapprox(output[13].d[9], 0.0000000e+00, rtol=1.0e-6)
    @test isapprox(output[13].d[6], 1.9508222e-05, rtol=1.0e-6)

    # Test Case 14
    @test isapprox(output[14].t[1], 1027.32, atol=1.0e-2)
    @test isapprox(output[14].t[2], 279.56, atol=1.0e-2)
    @test isapprox(output[14].d[1], 1.4121835e+11, rtol=1.0e-6)
    @test isapprox(output[14].d[2], 0.0000000e+00, rtol=1.0e-6)
    @test isapprox(output[14].d[3], 2.1043696e+16, rtol=1.0e-6)
    @test isapprox(output[14].d[4], 5.6453924e+15, rtol=1.0e-6)
    @test isapprox(output[14].d[5], 2.5171417e+14, rtol=1.0e-6)
    @test isapprox(output[14].d[7], 0.0000000e+00, rtol=1.0e-6)
    @test isapprox(output[14].d[8], 0.0000000e+00, rtol=1.0e-6)
    @test isapprox(output[14].d[9], 0.0000000e+00, rtol=1.0e-6)
    @test isapprox(output[14].d[6], 1.2947090e-06, rtol=1.0e-6)

    # Test Case 15
    @test isapprox(output[15].t[1], 1027.32, atol=1.0e-2)
    @test isapprox(output[15].t[2], 219.07, atol=1.0e-2)
    @test isapprox(output[15].d[1], 1.2548844e+10, rtol=1.0e-6)
    @test isapprox(output[15].d[2], 0.0000000e+00, rtol=1.0e-6)
    @test isapprox(output[15].d[3], 1.8745328e+15, rtol=1.0e-6)
    @test isapprox(output[15].d[4], 4.9230510e+14, rtol=1.0e-6)
    @test isapprox(output[15].d[5], 2.2396854e+13, rtol=1.0e-6)
    @test isapprox(output[15].d[7], 0.0000000e+00, rtol=1.0e-6)
    @test isapprox(output[15].d[8], 0.0000000e+00, rtol=1.0e-6)
    @test isapprox(output[15].d[9], 0.0000000e+00, rtol=1.0e-6)
    @test isapprox(output[15].d[6], 1.1476677e-07, rtol=1.0e-6)
end

let
    # Date and conversions
    epc = Epoch(2018, 1, 1, 12, 0, 0)
    oe = [R_EARTH + 500e3, 0.001, 97.8, 75, 25, 45]
    eci  = sOSCtoCART(oe, use_degrees=true)
    ecef = sECItoECEF(epc, eci)
    geod = sECEFtoGEOD(ecef, use_degrees=true)

    rho = density_nrlmsise00(epc, geod, use_degrees=true)

    @test 0.0 <= rho <= 1.0e-13
end