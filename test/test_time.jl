let 
    @test caldate_to_mjd(2000, 1, 1, 12, 0, 0) == 51544.5
end

let
    year, month, day, hour, minute, second, microsecond = mjd_to_caldate(51544.5)
    
    @test year        == 2000
    @test month       == 1
    @test day         == 1
    @test hour        == 12
    @test minute      == 0
    @test second      == 0
    @test microsecond == 0
end

let
    @test caldate_to_jd(2000, 1, 1, 12, 0, 0) == 2451545.0
end

let
    year, month, day, hour, minute, second, microsecond = jd_to_caldate(2451545.0)
    @test year        == 2000
    @test month       == 1
    @test day         == 1
    @test hour        == 12
    @test minute      == 0
    @test second      == 0
    @test microsecond == 0
end

let
    # MJD_ZERO = SatelliteDynamics.MJD_ZERO
    t = elapsed_from_epoch(MJD_ZERO + 175.500011574074074074073e-5, MJD_ZERO)
    @test isapprox(t, 15163200.000011574074074074073e-5, rtol=1e-6)

    t = elapsed_from_epoch(MJD_ZERO + 365.25*20, MJD_ZERO)
    @test isapprox(t, 631152000, rtol=1e-6)
end

let
    mjd = days_from_elapsed(15163200.000011574074074074073e-5, 0)
    @test isapprox(mjd, 175.500011574074074074073e-5, rtol=1e-3)

    mjd = days_from_elapsed(631152000, 0)
    @test isapprox(mjd, 7305, rtol=1e-3)
end

# let
#     jd, fd = SatelliteDynamics.align_jd_fd(1, 1.1)

#     @test isapprox(jd, 2, rtol=1e-12)
#     @test isapprox(fd, 0.1, rtol=1e-12)
# end

# let
#     mjd = caldate_to_mjd(2018, 6, 1)

#     dutc = -37.0

#     # GPS
#     @test time_system_offset(mjd, "GPS", "GPS") == 0
#     @test time_system_offset(mjd, "GPS", "TT")  == TT_GPS
#     @test time_system_offset(mjd, "GPS", "UTC") == -18
#     @test time_system_offset(mjd, "GPS", "UT1") == -17.92267
#     @test time_system_offset(mjd, "GPS", "TAI") == TAI_GPS

#     # TT
#     @test time_system_offset(mjd, "TT", "GPS") == GPS_TT
#     @test time_system_offset(mjd, "TT", "TT")  == 0
#     @test time_system_offset(mjd, "TT", "UTC") == dutc + TAI_TT
#     @test time_system_offset(mjd, "TT", "UT1") == -69.10667
#     @test time_system_offset(mjd, "TT", "TAI") == TAI_TT

#     # UTC
#     @test time_system_offset(mjd, "UTC", "GPS") == 18
#     @test time_system_offset(mjd, "UTC", "TT")  == -dutc + TT_TAI
#     @test time_system_offset(mjd, "UTC", "UTC") == 0.0
#     @test time_system_offset(mjd, "UTC", "UT1") == 0.0769964
#     @test time_system_offset(mjd, "UTC", "TAI") == -dutc

#     # UT1
#     @test time_system_offset(mjd, "UT1", "GPS") == 17.9230036
#     @test time_system_offset(mjd, "UT1", "TT")  == 69.1070036
#     @test isapprox(time_system_offset(mjd, "UT1", "UTC"), -0.0769964, atol=1e-7)
#     @test time_system_offset(mjd, "UT1", "UT1") == 0
#     @test time_system_offset(mjd, "UT1", "TAI") == 36.9230036

#     # TAI
#     @test time_system_offset(mjd, "TAI", "GPS") == GPS_TAI
#     @test time_system_offset(mjd, "TAI", "TT")  == TT_TAI
#     @test time_system_offset(mjd, "TAI", "UTC") == dutc
#     @test time_system_offset(mjd, "TAI", "UT1") == -36.92267
#     @test time_system_offset(mjd, "TAI", "TAI") == 0
# end

# let
#     epc = Epoch(2000, 1, 1, 12, 0, 0, tsys="TAI")

#     @test epc.jd == 2451544.5
#     @test epc.fd == 0.5

#     year, month, day, hour, minute, second = cal(epc)
#     @test year   == 2000
#     @test month  == 1
#     @test day    == 1
#     @test hour   == 12
#     @test minute == 0
#     @test second == 0

#     @test mjd(epc) == 51544.5

#     @test floor(day_of_year(epc)) == 1

#     @test 0.0 <= gast(epc) <= 2*pi 
#     @test 0.0 <= gmst(epc) <= 2*pi
# end

# let
#     epc  = Epoch(2018, 10, 11, 8, 42, 23)
#     epc += 60
#     @test epc == Epoch(2018, 10, 11, 8, 43, 23)

#     epc = epc - 23


#     @test epc == Epoch(2018, 10, 11, 8, 43, 0)

#     epc = Epoch(2018, 12, 31, 0, 0, 0)
#     @test epc == Epoch(2018, 12, 31, 0, 0, 0)
#     for i in 1:86400
#         epc += 1.0
#     end
#     @test epc == Epoch(2019, 1, 1, 0, 0, 0)

#     epc = Epoch(2019, 1, 1, 0, 0, 0)
#     epc += 365.0*86400.0 + 0.001
#     @test epc == Epoch(2020, 1, 1, 0, 0, 0)
# end