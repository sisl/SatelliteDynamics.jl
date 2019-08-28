"""
Test function of a parabola
"""
function parabola(t, x; active::Bool=false)
    if active == true
        return [3*t^2]
    else
        return [0.0]
    end
end

"""
Simple test function for spherical earth dynamics
"""
function point_earth(epc::Epoch, x::Array{<:Real, 1})
    # Restrict inputs to position only. Considered in body frame
    r  = x[1:3]
    v  = x[4:6]

    # Acceleration
    a = -GM_EARTH * r/norm(r)^3

    return vcat(v, a)
end

# Test basic integration of parabola 
let
    rk4 = RK4(parabola, active=false)
    @test istep(rk4, 0.0, 1.0, [0.0])[1] == 0.0

    rk4 = RK4(parabola, active=true)
    @test istep(rk4, 0.0, 1.0, [0.0])[1] == 1.0

    t = 0.0
    for i in 1:10
        t = istep(rk4, 0, i, [0.0])[1]
    end

    @test t == 1000
end

# Test stability of integration of conservative motion (point-mass orbit)
let
    rk4 = RK4(point_earth)

    epc0 = Epoch(2019, 1, 1)
    oe0  = [R_EARTH + 500e3, 0.01, 90.0, 0, 0, 0]
    eci0 = sOSCtoCART(oe0, use_degrees=true)

    epcf = epc0 + orbit_period(oe0[1])
    
    # Manually propagate orbit
    epc = deepcopy(epc0)
    eci = copy(eci0)
    while epc < epcf
        dt  = min(1.0, epcf-epc) # Set step-size
        eci = istep(rk4, epc, dt, eci)
        epc += dt
    end

    # Test that initial state equals intial propagated state
    @test isapprox(eci[1, 1], eci0[1], atol=1e-5)
    @test isapprox(eci[2, 1], eci0[2], atol=1e-5)
    @test isapprox(eci[3, 1], eci0[3], atol=1e-5)
    @test isapprox(eci[4, 1], eci0[4], atol=1e-5)
    @test isapprox(eci[5, 1], eci0[5], atol=1e-5)
    @test isapprox(eci[6, 1], eci0[6], atol=1e-5)

    oe = sCARTtoOSC(eci[:, 1], use_degrees=true)
    @test isapprox(oe[1], oe0[1], atol=1e-5)
    @test isapprox(oe[2], oe0[2], atol=1e-5)
    @test isapprox(oe[3], oe0[3], atol=1e-5)
    @test isapprox(oe[4], oe0[4], atol=1e-5)
    @test isapprox(oe[5], oe0[5], atol=1e-5)
    @test isapprox(oe[6], oe0[6], atol=1e-5)
end

let
    # Initialize integrator
    rk4 = RK4(point_earth)

    epc = Epoch(2019, 1, 1)
    eci = sOSCtoCART([R_EARTH + 500e3, 0.01, 90.0, 0, 0, 0], use_degrees=true)
    phi = diagm(0 => ones(Float64, length(eci)))

    # No step variational equations STM should be identity
    xu, phiu = istep(rk4, epc, 0.0, eci, phi)

    for i in 1:length(eci)
        @test phiu[i,i] == 1.0
    end

    # Small step STM should be nearly identity but non-zero
    xu, phiu = istep(rk4, epc, 1.0, eci, phi)

    for i in 1:length(eci)
        @test phiu[i, i] != 0.0
        @test isapprox(phiu[i, i], 1.0, atol=1.0e-5)
    end

    # Actually performing a perturbed step should equal the STM output
    pert = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    # Perturbed position from integrating pertturbed initial conditions
    ecip   = istep(rk4, epc, 1.0, eci + pert)

    # Pertrubed position from variational equation approximation
    ecistm = istep(rk4, epc, 1.0, eci) + phiu*pert

    # Assert equal
    for i in 1:6
        @test isapprox(ecip[i], ecistm[i], atol=1.0e-9)
    end
end
