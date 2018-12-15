let
    @test C_LIGHT     == 299792458.0
    @test AU          == 1.49597870700e11
    @test MJD_ZERO    == 2400000.5
    @test MJD2000     == 51544.0
    @test GPS_TAI     == -19.0
    @test TAI_GPS     == -GPS_TAI
    @test TT_TAI      == 32.184
    @test TAI_TT      == -TT_TAI
    @test GPS_TT      == GPS_TAI + TAI_TT
    @test TT_GPS      == -GPS_TT
    @test GPS_ZERO    == 44244.0
    @test R_EARTH     == 6.378136300e6
    @test WGS84_a     == 6378137.0
    @test WGS84_f     == 1.0/298.257223563
    @test GM_EARTH    == 3.986004415e14
    @test e_EARTH     == 8.1819190842622e-2
    @test J2_EARTH    == 0.0010826358191967
    @test OMEGA_EARTH == 7.292115146706979e-5
    @test GM_SUN      == 132712440041.939400*1e9
    @test R_SUN       == 6.957*1e8
    @test P_SUN       == 4.560E-6
    @test GM_MOON     == 4902.800066*1e9
    @test GM_MERCURY  == 22031.780000*1e9
    @test GM_VENUS    == 324858.592000*1e9
    @test GM_MARS     == 42828.37521*1e9
    @test GM_JUPITER  == 126712764.8*1e9
    @test GM_SATURN   == 37940585.2*1e9
    @test GM_URANUS   == 5794548.6*1e9
    @test GM_NEPTUNE  == 6836527.100580*1e9
    @test GM_PLUTO    == 977.000000*1e9
end