let 
    # Test initialization of orbit
    epc0 = Epoch(2019, 1, 1, 12, 0, 0, 0.0) 
    oe0  = [R_EARTH + 500e3, 0.0, 90.0, 0, 0, 0]
    eci0 = sOSCtoCART(oe0, use_degrees=true)

    # Simulate orbit for one orbit
    epcf = epc0 + orbit_period(oe0[1])
    t, epc, eci = propagate_orbit(epc0, eci0, epcf, solver=Tsit5(), timestep=0.1, dtmax=1)

    # Test that initial state equals intial propagated state
    @test eci[1, 1] == eci0[1]
    @test eci[2, 1] == eci0[2]
    @test eci[3, 1] == eci0[3]
    @test eci[4, 1] == eci0[4]
    @test eci[5, 1] == eci0[5]
    @test eci[6, 1] == eci0[6]

    oe02 = sCARTtoOSC(eci[:, 1], use_degrees=true)
    @test isapprox(oe02[1], oe0[1], atol=1e-9)
    @test isapprox(oe02[2], oe0[2], atol=1e-9)
    @test isapprox(oe02[3], oe0[3], atol=1e-9)
    @test isapprox(oe02[4], oe0[4], atol=1e-9)
    @test isapprox(oe02[5], oe0[5], atol=1e-9)
    @test isapprox(oe02[6], oe0[6], atol=1e-9)

    # Test that initial state equals final state after one orbit
    # Exact equality isn't possible due to numerical integration errors
    tol=1e-4
    @test isapprox(eci[1, end], eci0[1], atol=tol)
    @test isapprox(eci[2, end], eci0[2], atol=tol)
    @test isapprox(eci[3, end], eci0[3], atol=tol)
    @test isapprox(eci[4, end], eci0[4], atol=tol)
    @test isapprox(eci[5, end], eci0[5], atol=tol)
    @test isapprox(eci[6, end], eci0[6], atol=tol)
end

# Test all options enabled
let
    # Test initialization of orbit
    epc0 = Epoch(2019, 1, 1, 12, 0, 0, 0.0) 
    oe0  = [R_EARTH + 500e3, 0.0, 90.0, 0, 0, 0]
    eci0 = sOSCtoCART(oe0, use_degrees=true)

    # Simulate orbit for one orbit
    epcf = epc0 + orbit_period(oe0[1])
    t, epc, eci = propagate_orbit(epc0, eci0, epcf, solver=Tsit5(), timestep=0.1, dtmax=1,
    mass=100.0, area_drag=1.0, coef_drag=2.3, arear_srp=1.0, coef_srp=1.8, n_grav=20, m_grav=20, drag=true, srp=true, moon=true, sun=true, relativity=true)
end