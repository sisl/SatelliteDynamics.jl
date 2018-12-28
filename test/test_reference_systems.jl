let
    epc = Epoch(2018,1,1,12,0,0)

    a   = R_EARTH + 500e3
    oe  = [a, 0, 0, 0, 0, sqrt(GM_EARTH/a)]
    eci = sOSCtoCART(oe, use_degrees=true)

    r_rtn   = rRTNtoECI(eci)
    r_rtn_t = rECItoRTN(eci)

    r = r_rtn * r_rtn_t
    
    tol = 1e-8
    @test isapprox(r[1, 1], 1.0, atol=tol)
    @test isapprox(r[1, 2], 0.0, atol=tol)
    @test isapprox(r[1, 3], 0.0, atol=tol)

    @test isapprox(r[2, 1], 0.0, atol=tol)
    @test isapprox(r[2, 2], 1.0, atol=tol)
    @test isapprox(r[2, 3], 0.0, atol=tol)

    @test isapprox(r[3, 1], 0.0, atol=tol)
    @test isapprox(r[3, 2], 0.0, atol=tol)
    @test isapprox(r[3, 3], 1.0, atol=tol)
end

let 
    epc = Epoch(2018,1,1,12,0,0)

    oe  = [R_EARTH + 500e3, 0, 0, 0, 0, 0]
    eci = sOSCtoCART(oe, use_degrees=true)

    xt = deepcopy(eci) + [100, 0, 0, 0, 0, 0]

    x_rtn = sECItoRTN(eci, xt)

    tol = 1e-8
    @test isapprox(x_rtn[1], 100.0, atol=tol)
    @test isapprox(x_rtn[2], 0.0, atol=tol)
    @test isapprox(x_rtn[3], 0.0, atol=tol)
    @test isapprox(x_rtn[4], 0.0, atol=tol)
    @test isapprox(x_rtn[5], 0.0, atol=0.5)
    @test isapprox(x_rtn[6], 0.0, atol=tol)

    xt2 = sRTNtoECI(eci, x_rtn)

    @test isapprox(xt[1], xt2[1], atol=tol)
    @test isapprox(xt[2], xt2[2], atol=tol)
    @test isapprox(xt[3], xt2[3], atol=tol)
    @test isapprox(xt[4], xt2[4], atol=tol)
    @test isapprox(xt[5], xt2[5], atol=tol)
    @test isapprox(xt[6], xt2[6], atol=tol)
end

let 
    epc = Epoch(2007, 4, 5, 12, 0, 0, tsys=:UTC)

    set_eop(54195.5, -0.072073685, 0.0349282, 0.4833163)

    rc2i = bias_precession_nutation(epc)

    tol = 1e-8
    @test isapprox(rc2i[1, 1], +0.999999746339445, atol=tol)
    @test isapprox(rc2i[1, 2], -0.000000005138822, atol=tol)
    @test isapprox(rc2i[1, 3], -0.000712264730072, atol=tol)

    @test isapprox(rc2i[2, 1], -0.000000026475227, atol=tol)
    @test isapprox(rc2i[2, 2], +0.999999999014975, atol=tol)
    @test isapprox(rc2i[2, 3], -0.000044385242827, atol=tol)

    @test isapprox(rc2i[3, 1], +0.000712264729599, atol=tol)
    @test isapprox(rc2i[3, 2], +0.000044385250426, atol=tol)
    @test isapprox(rc2i[3, 3], +0.999999745354420, atol=tol)
end

let 
    epc = Epoch(2007, 4, 5, 12, 0, 0, tsys=:UTC)

    set_eop(54195.5, -0.072073685, 0.0349282, 0.4833163)

    r = earth_rotation(epc) * bias_precession_nutation(epc)

    tol = 1e-8
    @test isapprox(r[1, 1], +0.973104317573127, atol=tol)
    @test isapprox(r[1, 2], +0.230363826247709, atol=tol)
    @test isapprox(r[1, 3], -0.000703332818845, atol=tol)

    @test isapprox(r[2, 1], -0.230363798804182, atol=tol)
    @test isapprox(r[2, 2], +0.973104570735574, atol=tol)
    @test isapprox(r[2, 3], +0.000120888549586, atol=tol)

    @test isapprox(r[3, 1], +0.000712264729599, atol=tol)
    @test isapprox(r[3, 2], +0.000044385250426, atol=tol)
    @test isapprox(r[3, 3], +0.999999745354420, atol=tol)
end

let 
    epc = Epoch(2007, 4, 5, 12, 0, 0, tsys=:UTC)

    set_eop(54195.5, -0.072073685, 0.0349282, 0.4833163)

    r = rECItoECEF(epc)

    tol = 1e-8
    @test isapprox(r[1, 1], +0.973104317697535, atol=tol)
    @test isapprox(r[1, 2], +0.230363826239128, atol=tol)
    @test isapprox(r[1, 3], -0.000703163482198, atol=tol)

    @test isapprox(r[2, 1], -0.230363800456037, atol=tol)
    @test isapprox(r[2, 2], +0.973104570632801, atol=tol)
    @test isapprox(r[2, 3], +0.000118545366625, atol=tol)

    @test isapprox(r[3, 1], +0.000711560162668, atol=tol)
    @test isapprox(r[3, 2], +0.000046626403995, atol=tol)
    @test isapprox(r[3, 3], +0.999999745754024, atol=tol)
end

let 
    epc = Epoch(2007, 4, 5, 12, 0, 0, tsys=:UTC)

    set_eop(54195.5, -0.072073685, 0.0349282, 0.4833163)

    r = rECEFtoECI(epc)

    tol = 1e-8
    @test isapprox(r[1, 1], +0.973104317697535, atol=tol)
    @test isapprox(r[1, 2], -0.230363800456037, atol=tol)
    @test isapprox(r[1, 3], +0.000711560162668, atol=tol)

    @test isapprox(r[2, 1], +0.230363826239128, atol=tol)
    @test isapprox(r[2, 2], +0.973104570632801, atol=tol)
    @test isapprox(r[2, 3], +0.000046626403995, atol=tol)

    @test isapprox(r[3, 1], -0.000703163482198, atol=tol)
    @test isapprox(r[3, 2], +0.000118545366625, atol=tol)
    @test isapprox(r[3, 3], +0.999999745754024, atol=tol)
end

let
    epc = Epoch(2018,1,1,12,0,0)

    oe  = [R_EARTH + 500e3, 1e-3, 97.8, 75, 25, 45]
    eci = sOSCtoCART(oe, use_degrees=true)

    # Perform circular transformations
    ecef  = sECItoECEF(epc, eci)
    eci2  = sECEFtoECI(epc, ecef)
    ecef2 = sECItoECEF(epc, eci2)

    tol=1e-6
    # Check equivalence of ECI transforms
    @test isapprox(eci2[1], eci[1],  atol=tol)
    @test isapprox(eci2[2], eci[2],  atol=tol)
    @test isapprox(eci2[3], eci[3],  atol=tol)
    @test isapprox(eci2[4], eci[4],  atol=tol)
    @test isapprox(eci2[5], eci[5],  atol=tol)
    @test isapprox(eci2[6], eci[6],  atol=tol)

    # Check equivalence of ECEF transforms
    @test isapprox(ecef2[1], ecef[1], atol=tol)
    @test isapprox(ecef2[2], ecef[2], atol=tol)
    @test isapprox(ecef2[3], ecef[3], atol=tol)
    @test isapprox(ecef2[4], ecef[4], atol=tol)
    @test isapprox(ecef2[5], ecef[5], atol=tol)
    @test isapprox(ecef2[6], ecef[6], atol=tol)
end