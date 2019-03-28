let
    epc   = Epoch(2018, 3, 20, 16, 15, 0) # Test on Vernal equinox
    oe    = [R_EARTH + 500e3, 0, 90.0, 0, 0, 0]
    x     = sOSCtoCART(oe, use_degrees = true)
    R     = rECItoECEF(epc)

    a_grav_point = accel_point_mass(x[1:3])

    a_grav_00 = accel_gravity(x, R, 0, 0)

    tol=1e-9
    @test isapprox(a_grav_00[1], a_grav_point[1], atol=tol)
    @test isapprox(a_grav_00[2], a_grav_point[2], atol=tol)
    @test isapprox(a_grav_00[3], a_grav_point[3], atol=tol)

    a_grav_20 = accel_gravity(x, R, 20, 20)
    @test norm(a_grav_20) > norm(a_grav_00)
end

let 
    # Define Initial State
    epc   = Epoch(2018, 3, 20, 16, 15, 0) # Test on Vernal equinox
    oe    = [R_EARTH + 500e3, 0, 90.0, 0, 0, 0]
    x     = sOSCtoCART(oe, use_degrees = true)
    r_sun = sun_position(epc)

    # Call function
    a_sun = accel_thirdbody_sun(epc, x)

    @test a_sun[1] < 1e-5
    @test a_sun[2] < 1e-5
    @test a_sun[3] < 1e-5

    a_sun2 = accel_thirdbody_sun(x, r_sun)

    @test a_sun[1] == a_sun2[1]
    @test a_sun[2] == a_sun2[2]
    @test a_sun[3] == a_sun2[3]
end

let 
    # Define Initial State
    epc    = Epoch(2018, 3, 20, 16, 15, 0) # Test on Vernal equinox
    oe     = [R_EARTH + 500e3, 0, 90.0, 0, 0, 0]
    x      = sOSCtoCART(oe, use_degrees = true)
    r_moon = moon_position(epc)

    # Call function
    a_moon = accel_thirdbody_moon(epc, x)

    @test a_moon[1] < 1e-5
    @test a_moon[2] < 1e-5
    @test a_moon[3] < 1e-5

    a_moon2 = accel_thirdbody_moon(x, r_moon)

    @test a_moon[1] == a_moon2[1]
    @test a_moon[2] == a_moon2[2]
    @test a_moon[3] == a_moon2[3]
end


let 
    # Define Initial State
    epc   = Epoch(2018, 3, 20, 16, 15, 0) # Test on Vernal equinox
    oe    = [R_EARTH + 500e3, 0, 90.0, 0, 0, 0]
    x     = sOSCtoCART(oe, use_degrees = true)
    T     = rECItoECEF(epc)
    
    # Call function
    rho    = density_harris_priester(epc, x)
    a_drag = accel_drag(x, rho, 1, 1, 2.3, T)

    @test a_drag[1] < 1e-5
    @test a_drag[2] < 1e-5
    @test a_drag[3] < 1e-5
end

let 
    # Define Initial State
    epc   = Epoch(2018, 3, 20, 16, 15, 0) # Test on Vernal equinox
    oe    = [R_EARTH + 500e3, 0, 0, 0, 0, 0]
    x     = sOSCtoCART(oe, use_degrees = true)
    r_sun = sun_position(epc)
    
    # Call function
    nu = eclipse_cylindrical(x, r_sun)
    @test nu == 1.0
    @test nu == eclipse_cylindrical(epc, x)

    oe    = [R_EARTH + 500e3, 0, 0, 180.0, 0, 0]
    x     = sOSCtoCART(oe, use_degrees = true)
    nu = eclipse_cylindrical(x, r_sun)
    @test nu == 0.0
    @test nu == eclipse_cylindrical(epc, x)
end

let 
    # Define Initial State
    epc   = Epoch(2018, 3, 20, 16, 15, 0) # Test on Vernal equinox
    oe    = [R_EARTH + 500e3, 0, 0, 0, 0, 0]
    x     = sOSCtoCART(oe, use_degrees = true)
    r_sun = sun_position(epc)
    
    # Call function
    nu = eclipse_conical(x, r_sun)
    @test nu == 0.0
    @test nu == eclipse_conical(epc, x)

    oe    = [R_EARTH + 500e3, 0, 0, 180.0, 0, 0]
    x     = sOSCtoCART(oe, use_degrees = true)
    nu = eclipse_conical(x, r_sun)
    @test nu == 1.0
    @test nu == eclipse_conical(epc, x)
end

let 
    # Define Initial State
    epc   = Epoch(2018, 1, 1, 0, 0, 0)
    oe    = [R_EARTH + 500e3, 0, 90.0, 0, 0, 0]
    x     = sOSCtoCART(oe, use_degrees = true)
    r_sun = sun_position(epc)
    
    # Call function
    a_srp = accel_srp(x, r_sun, 1, 1)

    @test a_srp[1] < 1e-5
    @test a_srp[2] < 1e-5
    @test a_srp[3] < 1e-5
end

let 
    # Define Initial State
    oe = [R_EARTH + 500e3, 0, 90.0, 0, 0, 0]
    x  = sOSCtoCART(oe, use_degrees=true)
    
    # Call function
    a_rel = accel_relativity(x)

    @test a_rel[1] < 1e-7
    @test a_rel[2] < 1e-7
    @test a_rel[3] < 1e-7
end