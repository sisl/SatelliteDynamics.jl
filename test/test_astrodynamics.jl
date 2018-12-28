let
    n = mean_motion(R_EARTH + 500e3, use_degrees=false)
    @test isapprox(n, 0.0011067836148773837, atol=1e-12)

    n = mean_motion(R_EARTH + 500e3, use_degrees=true)
    @test isapprox(n, 0.0634140299667068, atol=1e-12)
end

let
    a = semimajor_axis(0.0011067836148773837, use_degrees=false)
    @test isapprox(a, R_EARTH + 500e3, atol=1e-6)

    a = semimajor_axis(0.0634140299667068, use_degrees=true)
    @test isapprox(a, R_EARTH + 500e3, atol=1e-6)
end

let 
    T = orbit_period(R_EARTH + 500e3)
    @test isapprox(T, 5676.977164028288, atol=1e-9)
end

let
    iss = sun_sync_inclination(R_EARTH + 574e3, 0.0, use_degrees=false)
    @test isapprox(iss, 97.685*pi/180, atol=1.0e-3)

    iss = sun_sync_inclination(R_EARTH + 574e3, 0.0, use_degrees=true)
    @test isapprox(iss, 97.685, atol=1.0e-3)
end

let
    # 0 
    M = anomaly_eccentric_to_mean(0.0, 0.0, use_degrees=false)
    @test M == 0

    M = anomaly_eccentric_to_mean(0.0, 0.0, use_degrees=true)
    @test M == 0

    # 180
    M = anomaly_eccentric_to_mean(pi/2, 0.1, use_degrees=false)
    @test isapprox(M, 1.4707, atol=1e-3)

    M = anomaly_eccentric_to_mean(90.0, 0.1, use_degrees=true)
    @test isapprox(M, 84.270, atol=1e-3)

    # 180
    M = anomaly_eccentric_to_mean(pi, 0.0, use_degrees=false)
    @test isapprox(M, pi, atol=1e-12)

    M = anomaly_eccentric_to_mean(180.0, 0.0, use_degrees=true)
    @test M == 180.0
end

let
    # 0 
    E = anomaly_mean_to_eccentric(0.0, 0.0, use_degrees=false)
    @test E == 0

    E = anomaly_mean_to_eccentric(0.0, 0.0, use_degrees=true)
    @test E == 0

    # 180
    E = anomaly_mean_to_eccentric(1.4707, 0.1, use_degrees=false)
    @test isapprox(E, pi/2, atol=1e-3)

    E = anomaly_mean_to_eccentric(84.270, 0.1, use_degrees=true)
    @test isapprox(E, 90.0, atol=1e-3)

    # 180
    E = anomaly_mean_to_eccentric(pi, 0.0, use_degrees=false)
    @test isapprox(E, pi, atol=1e-12)

    E = anomaly_mean_to_eccentric(180.0, 0.0, use_degrees=true)
    @test E == 180.0

    # Large Eccentricities
    E = anomaly_mean_to_eccentric(180.0, 0.9, use_degrees=true)
    @test E == 180.0

    # Max Iteration exit on hyperbolic trajectory
    # @test_throws ErrorException anomaly_mean_to_eccentric(pi, 2.0)
end

let
    oe  = [R_EARTH + 500e3, 0, 90.0, 0, 0, 0]
    eci = sOSCtoCART(oe, use_degrees=true)

    tol = 1e-6
    @test isapprox(eci[1], R_EARTH + 500e3, atol=tol)
    @test isapprox(eci[2], 0.0, atol=tol)
    @test isapprox(eci[3], 0.0, atol=tol)
    @test isapprox(eci[4], 0.0, atol=tol)
    @test isapprox(eci[5], 0.0, atol=tol)
    @test isapprox(eci[6], sqrt(GM_EARTH/(R_EARTH + 500e3)), atol=tol)
end

let 
    # Using radians
    oe   = [R_EARTH + 500e3, 0, pi/2.0, 0, 0, 0]
    eci  = sOSCtoCART(oe, use_degrees=false)
    eci2 = sOSCtoCART(oe, use_degrees=false)

    tol = 1e-6
    @test isapprox(eci[1], eci2[1], atol=tol)
    @test isapprox(eci[2], eci2[2], atol=tol)
    @test isapprox(eci[3], eci2[3], atol=tol)
    @test isapprox(eci[4], eci2[4], atol=tol)
    @test isapprox(eci[5], eci2[5], atol=tol)
    @test isapprox(eci[6], eci2[6], atol=tol)  

    # Using degrees
    oe   = [R_EARTH + 500e3, 0, 90.0, 0, 0, 0]
    eci  = sOSCtoCART(oe, use_degrees=true)
    eci2 = sOSCtoCART(oe, use_degrees=true)

    tol = 1e-6
    @test isapprox(eci[1], eci2[1], atol=tol)
    @test isapprox(eci[2], eci2[2], atol=tol)
    @test isapprox(eci[3], eci2[3], atol=tol)
    @test isapprox(eci[4], eci2[4], atol=tol)
    @test isapprox(eci[5], eci2[5], atol=tol)
    @test isapprox(eci[6], eci2[6], atol=tol)
end

let 
    eci   = [R_EARTH + 500e3, 100e3, 575e3, 0, 0, 7300]
    eci2  = sOSCtoCART(sCARTtoOSC(eci))

    tol = 1e-6
    @test isapprox(eci[1], eci2[1], atol=tol)
    @test isapprox(eci[2], eci2[2], atol=tol)
    @test isapprox(eci[3], eci2[3], atol=tol)
    @test isapprox(eci[4], eci2[4], atol=tol)
    @test isapprox(eci[5], eci2[5], atol=tol)
    @test isapprox(eci[6], eci2[6], atol=tol)

    # Equatorial circulator
    a   = R_EARTH + 1000e3
    e   = 0.0
    eci = [a, 0, 0, 0, 0, sqrt(GM_EARTH/a)]
    oe  = sCARTtoOSC(eci, use_degrees=true)

    tol = 1e-6
    @test isapprox(oe[1], a,     atol=tol)
    @test isapprox(oe[2], e,     atol=tol)
    @test isapprox(oe[3], 90.0,  atol=tol)
end

let
    # Test near-circular conversions
    a   = R_EARTH + 500e3
    eci = [a, 0.0, 0.0, 0.0, 0.0, sqrt(GM_EARTH/a)]
    oe  = sCARTtoOSC(eci, use_degrees=true)

    tol = 1e-6
    @test isapprox(oe[1], a, atol=tol)
    @test isapprox(oe[2], 0.0, atol=tol)
    @test isapprox(oe[3], 90.0, atol=tol)
    @test isapprox(oe[4], 0.0, atol=tol)
    @test isapprox(oe[5], 0.0, atol=tol)
    @test isapprox(oe[6], 0.0, atol=tol)
end