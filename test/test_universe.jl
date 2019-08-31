##########################
# Earth Orientation Data #
##########################

# Test internal data type
let 
    @test typeof(EOP.data) == Dict{Int,Tuple{Float64,Float64,Float64}}
end

# Test constructors
let
    finals_2000 = EarthOrientationData("FINALS_2000")
    eop_c04     = EarthOrientationData("C04_14")
end

# Test Data valutes
let
    @test isapprox(UT1_UTC(58483), -0.0351948, atol=1e-5)
    xp, yp = POLE_LOCATOR(58483)
    @test isapprox(xp, 0.088504*AS2RAD, atol=1e-5)
    @test isapprox(yp, 0.270753*AS2RAD, atol=1e-5)
    @test isapprox(XP(58483), 0.088504*AS2RAD, atol=1e-5)
    @test isapprox(YP(58483), 0.270753*AS2RAD, atol=1e-5)
end

# Test Set values
let
    set_eop(58747, -0.2, 0.225, 0.3)

    @test UT1_UTC(58747)      == -0.2
    @test POLE_LOCATOR(58747) == (0.225*AS2RAD, 0.3*AS2RAD)
    @test XP(58747)           == 0.225*AS2RAD
    @test YP(58747)           == 0.3*AS2RAD

    set_eop(58747, -0.1897929, 0.230292, 0.332704)
    @test UT1_UTC(58747)      == -0.1897929
    @test POLE_LOCATOR(58747) == (0.230292*AS2RAD, 0.332704*AS2RAD)
    @test XP(58747)           == 0.230292*AS2RAD
    @test YP(58747)           == 0.332704*AS2RAD
end

# Test interpolation values
let
    ut1_utc = (UT1_UTC(58748) + UT1_UTC(58747))/2.0
    @test UT1_UTC(58747.5, interp=true) == ut1_utc

    x1, y1 = POLE_LOCATOR(58747)
    x2, y2 = POLE_LOCATOR(58748)
    pole_locator = ((x2+x1)/2.0, (y2+y1)/2.0)
    @test POLE_LOCATOR(58747.5, interp=true)[1] == pole_locator[1]
    @test POLE_LOCATOR(58747.5, interp=true)[2] == pole_locator[2]

    xp = (XP(58748) + XP(58747))/2.0
    @test XP(58747.5, interp=true) == xp

    yp = (YP(58748) + YP(58747))/2.0
    @test YP(58747.5, interp=true) == yp
end

# Test loading new products into project global vairable
let 
    load_eop("C04_14")
    @test UT1_UTC(37665) == 0.0326338

    # Return to original state
    load_eop("FINALS_2000")
end

# Test Update of EOP
# let
#     update_eop("C04_14")
#     update_eop("C04_80")
#     update_eop("FINALS_2000")
# end

#################
# Gravity Model #
#################

let
    @test typeof(GRAVITY_MODEL.data) == Array{Float64, 2}

    @test GRAVITY_MODEL.name == "EGM2008"
    @test GRAV_COEF(0, 0)   == 1.0
    @test GRAV_COEF(90, 90) == 0.733188520723327e-09
end

let
    gmdl = GravModel(:EGM2008_90)
end

let 
    load_gravity_model(:GGM01S)
end