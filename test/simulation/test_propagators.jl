let
    # Initial Conditions
    epc = Epoch(2019, 1, 1)
    oe  = [R_EARTH + 500e3, 0.01, 90.0, 45.0, 75.0, 180.0]
    eci = sOSCtoCART(oe, use_degrees=true)

    # Initialize orbit
    orb = EarthInertialState(epc, eci, dt=1.0,
            mass=1.0, n_grav=0, m_grav=0,
            drag=false, srp=false,
            moon=false, sun=false,
            relativity=false
    )

    # Set propagation period
    epcf = epc + orbit_period(oe[1])

    # Run propagation
    stepto!(orb, epcf)

    # Final state should exactly equal the inittial state
    ecif = orb.x
    oef  = sCARTtoOSC(orb.x, use_degrees=true)
    for i in 1:6
        @test isapprox(ecif[i], eci[i], atol=1e-5)
        @test isapprox(oef[i], oe[i], atol=1e-5)
    end
end

# Test integration of variational equations
let
    # Initial Conditions
    epc = Epoch(2019, 1, 1)
    oe  = [R_EARTH + 500e3, 0.01, 90.0, 45.0, 75.0, 180.0]
    eci = sOSCtoCART(oe, use_degrees=true)

    # Initialize orbit
    orb = EarthInertialState(epc, eci, dt=1.0)
    reinit!(orb) # Initialize STM to enable propagation

    # Set propagation period
    epcf = epc + orbit_period(oe[1])
    
    # Run propagation
    stepto!(orb, epcf)

    # Test to confirm non-zero STM
    for i in 1:length(eci)
        for j in 1:length(eci)
            @test orb.phi[i,j] != 0.0
        end
    end
end

# Test agreegation of output information
let
    # Initial Conditions
    epc = Epoch(2019, 1, 1)
    oe  = [R_EARTH + 500e3, 0.01, 90.0, 45.0, 75.0, 180.0]
    eci = sOSCtoCART(oe, use_degrees=true)

    # Initialize orbit
    orb = EarthInertialState(epc, eci, dt=1.0)
    reinit!(orb) # Initialize STM to enable propagation

    # Set propagation period
    epcf = 5

    # Simulate output to agregate state information
    sim_t, sim_epc, sim_x, sim_phi = sim!(orb, epcf)

    for i in 1:(epcf+1)
        @test sim_t[i] == i-1
    end
end