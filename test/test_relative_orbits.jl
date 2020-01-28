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