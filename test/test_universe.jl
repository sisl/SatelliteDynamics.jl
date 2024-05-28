##########################
# Earth Orientation Data #
##########################

# Test internal data type
let 
    @test typeof(EOP.data) == Dict{Int,Tuple{Float64,Float64,Float64}}
end

# Test initialization from bundled data
let
    finals_2000 = EarthOrientationData(:FINALS_2000)
    eop_c04     = EarthOrientationData(:C04_20)
end

# Load static test data
let 
    load_eop(TEST_FINALS_EOP_DATA, :FINALS_2000)
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
    set_eop(1, -0.2, 0.225, 0.3)

    @test UT1_UTC(1)      == -0.2
    @test POLE_LOCATOR(1) == (0.225, 0.3)
    @test XP(1)           == 0.225
    @test YP(1)           == 0.3

    set_eop(1, -0.1897929, 0.230292, 0.332704)
    @test UT1_UTC(1)      == -0.1897929
    @test POLE_LOCATOR(1) == (0.230292, 0.332704)
    @test XP(1)           == 0.230292
    @test YP(1)           == 0.332704

    # Reset global EOP values
    load_eop(TEST_FINALS_EOP_DATA, :FINALS_2000)
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

# Test usage of C04 test data
let
    # Test loading new products into project global vairable
    load_eop(TEST_C04_EOP_DATA, :C04_20)

    # Confirm output values are expected
    @test UT1_UTC(60149)      == -0.0195468
    @test POLE_LOCATOR(60149) == (0.242675*AS2RAD, 0.487119*AS2RAD)
    @test XP(60149)           == 0.242675*AS2RAD
    @test YP(60149)           == 0.487119*AS2RAD

    # Test dX and dY values
    @test DX(60149) ==  0.000375 * AS2RAD
    @test DY(60149) == -0.000186 * AS2RAD

    # Return to original state
    load_eop(TEST_FINALS_EOP_DATA, :FINALS_2000)
end

# Test paring and retrieval of dX, dY values from FINALS 2000
let
    load_eop(TEST_FINALS_EOP_DATA, :FINALS_2000)
    
    # Confirm output values are expected
    @test UT1_UTC(60149)      == -0.0027438
    @test POLE_LOCATOR(60149) == (0.236429*AS2RAD, 0.474971*AS2RAD)
    @test XP(60149)           == 0.236429*AS2RAD
    @test YP(60149)           == 0.474971*AS2RAD

    # Test dX and dY values
    @test isapprox(DX(60149),  0.371 * 1.0e-3 * AS2RAD, rtol=1e-15)
    @test isapprox(DY(60149), -0.049 * 1.0e-3 * AS2RAD, rtol=1e-15)
end

# Test test global EOP from default data
let 
    load_eop(:C04_20)

    load_eop(:FINALS_2000)

    # Return to original state for futur tests
    load_eop(TEST_FINALS_EOP_DATA, :FINALS_2000)
end

# Test out-of-bounds error
let 
    try 
        UT1_UTC(99999)
    catch e
        @test e isa ErrorException
        # Test message contains
        @test occursin("is beyond the last", e.msg)
    end
    try 
        POLE_LOCATOR(99999)
    catch e
        @test e isa ErrorException
        # Test message contains
        @test occursin("is beyond the last", e.msg)
    end
    try 
        XP(99999)
    catch e
        @test e isa ErrorException
        # Test message contains
        @test occursin("is beyond the last", e.msg)
    end
    try 
        YP(99999)
    catch e
        @test e isa ErrorException
        # Test message contains
        @test occursin("is beyond the last", e.msg)
    end
end

# #################
# # Gravity Model #
# #################

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