let 
    @test typeof(EOP.data) == Dict{Int32,Tuple{Float64,Float64,Float64}}

    @test UT1_UTC(58747)      == -0.1897929
    @test POLE_LOCATOR(58747) == (0.230292*AS2RAD, 0.332704*AS2RAD)
    @test XP(58747)           == 0.230292*AS2RAD
    @test YP(58747)           == 0.332704*AS2RAD
end

let
    @test typeof(GRAVITY_MODEL.data) == Array{Float64, 2}

    @test GRAVITY_MODEL.name == "EGM2008"
    @test GRAV_COEF(0, 0)   == 1.0
    @test GRAV_COEF(90, 90) == 0.733188520723327e-09
end