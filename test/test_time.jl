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

let
    @test SatelliteDynamics.Time.valid_time_system("GPS") == true
    @test SatelliteDynamics.Time.valid_time_system("TAI") == true
    @test SatelliteDynamics.Time.valid_time_system("TT") == true
    @test SatelliteDynamics.Time.valid_time_system("UTC") == true
    @test SatelliteDynamics.Time.valid_time_system("UT1") == true

    @test SatelliteDynamics.Time.valid_time_system("BLAH") == false
end

let
    jd = caldate_to_jd(2018, 6, 1)

    dutc = -37.0

    # GPS
    @test time_system_offset(jd, 0, "GPS", "GPS") == 0
    @test time_system_offset(jd, 0, "GPS", "TT")  == TT_GPS
    @test time_system_offset(jd, 0, "GPS", "UTC") == -18
    @test isapprox(time_system_offset(jd, 0, "GPS", "UT1"), -17.92267, atol=1e-4)
    @test time_system_offset(jd, 0, "GPS", "TAI") == TAI_GPS

    # TT
    @test time_system_offset(jd, 0, "TT", "GPS") == GPS_TT
    @test time_system_offset(jd, 0, "TT", "TT")  == 0
    @test time_system_offset(jd, 0, "TT", "UTC") == dutc + TAI_TT
    @test isapprox(time_system_offset(jd, 0, "TT", "UT1"), -69.10667, atol=1e-4)
    @test time_system_offset(jd, 0, "TT", "TAI") == TAI_TT

    # UTC
    @test time_system_offset(jd, 0, "UTC", "GPS") == 18
    @test time_system_offset(jd, 0, "UTC", "TT")  == -dutc + TT_TAI
    @test time_system_offset(jd, 0, "UTC", "UTC") == 0.0
    @test isapprox(time_system_offset(jd, 0, "UTC", "UT1"), 0.0769968, atol=1e-4)
    @test time_system_offset(jd, 0, "UTC", "TAI") == -dutc

    # UT1
    @test isapprox(time_system_offset(jd, 0, "UT1", "GPS"), 17.9230032, atol=1e-4)
    @test isapprox(time_system_offset(jd, 0, "UT1", "TT"), 69.1070032, atol=1e-4)
    @test isapprox(time_system_offset(jd, 0, "UT1", "UTC"), -0.0769968, atol=1e-4)
    @test time_system_offset(jd, 0, "UT1", "UT1") == 0
    @test isapprox(time_system_offset(jd, 0, "UT1", "TAI"), 36.9230032, atol=1e-4)

    # TAI
    @test time_system_offset(jd, 0, "TAI", "GPS") == GPS_TAI
    @test time_system_offset(jd, 0, "TAI", "TT")  == TT_TAI
    @test time_system_offset(jd, 0, "TAI", "UTC") == dutc
    @test isapprox(time_system_offset(jd, 0, "TAI", "UT1"), -36.92267, atol=1e-4)
    @test time_system_offset(jd, 0, "TAI", "TAI") == 0
end

let
    epc = Epoch(2018, 12, 20, 0, 0, 0, 0.0, tsys="TAI")
    @test epc.days        == 2458472
    @test epc.seconds     == 43200
    @test epc.nanoseconds == 0.0
    @test epc.tsys        == "TAI"

    epc = Epoch(2018, 12, 20, 0, 0, .5, 1.0001, tsys="TAI")
    @test epc.days        == 2458472
    @test epc.seconds     == 43200
    @test epc.nanoseconds == 500000001.0001
    @test epc.tsys        == "TAI"
end

# String constructor
let
    epc = Epoch("2018-12-20")
    year, month, day, hour, minute, seconds, nanoseconds = caldate(epc)
    @test year        == 2018
    @test month       == 12
    @test day         == 20
    @test hour        == 0
    @test minute      == 0
    @test seconds     == 0
    @test nanoseconds == 0.0
    @test epc.tsys    == "UTC"

    epc = Epoch("2018-12-20T16:22:19.0Z")
    year, month, day, hour, minute, seconds, nanoseconds = caldate(epc)
    @test year        == 2018
    @test month       == 12
    @test day         == 20
    @test hour        == 16
    @test minute      == 22
    @test seconds     == 19
    @test nanoseconds == 0.0
    @test epc.tsys    == "UTC"

    epc = Epoch("2018-12-20T16:22:19.123Z")
    year, month, day, hour, minute, seconds, nanoseconds = caldate(epc)
    @test year        == 2018
    @test month       == 12
    @test day         == 20
    @test hour        == 16
    @test minute      == 22
    @test seconds     == 19
    @test nanoseconds == 123000000
    @test epc.tsys    == "UTC"

    epc = Epoch("2018-12-20T16:22:19.123456789Z")
    year, month, day, hour, minute, seconds, nanoseconds = caldate(epc)
    @test year        == 2018
    @test month       == 12
    @test day         == 20
    @test hour        == 16
    @test minute      == 22
    @test seconds     == 19
    @test nanoseconds == 123456789
    @test epc.tsys    == "UTC"

    epc = Epoch("2018-12-20T16:22:19Z")
    year, month, day, hour, minute, seconds, nanoseconds = caldate(epc)
    @test year        == 2018
    @test month       == 12
    @test day         == 20
    @test hour        == 16
    @test minute      == 22
    @test seconds     == 19
    @test nanoseconds == 0.0
    @test epc.tsys    == "UTC"

    epc = Epoch("20181220T162219Z")
    year, month, day, hour, minute, seconds, nanoseconds = caldate(epc)
    @test year        == 2018
    @test month       == 12
    @test day         == 20
    @test hour        == 16
    @test minute      == 22
    @test seconds     == 19
    @test nanoseconds == 0.0
    @test epc.tsys    == "UTC"

    epc = Epoch("2018-12-01 16:22:19 GPS")
    year, month, day, hour, minute, seconds, nanoseconds = caldate(epc)
    @test year        == 2018
    @test month       == 12
    @test day         == 1
    @test hour        == 16
    @test minute      == 22
    @test seconds     == 19
    @test nanoseconds == 0.0
    @test epc.tsys    == "GPS"

    epc = Epoch("2018-12-01 16:22:19.0 GPS")
    year, month, day, hour, minute, seconds, nanoseconds = caldate(epc)
    @test year        == 2018
    @test month       == 12
    @test day         == 1
    @test hour        == 16
    @test minute      == 22
    @test seconds     == 19
    @test nanoseconds == 0.0
    @test epc.tsys    == "GPS"

    epc = Epoch("2018-12-01 16:22:19.123 GPS")
    year, month, day, hour, minute, seconds, nanoseconds = caldate(epc)
    @test year        == 2018
    @test month       == 12
    @test day         == 1
    @test hour        == 16
    @test minute      == 22
    @test seconds     == 19
    @test nanoseconds == 123000000
    @test epc.tsys    == "GPS"

    epc = Epoch("2018-12-01 16:22:19.123456789 GPS")
    year, month, day, hour, minute, seconds, nanoseconds = caldate(epc)
    @test year        == 2018
    @test month       == 12
    @test day         == 1
    @test hour        == 16
    @test minute      == 22
    @test seconds     == 19
    @test nanoseconds == 123456789
    @test epc.tsys    == "GPS"

end

# Addition
let
    epc = Epoch("2019-01-01 12:00:00 TAI")
    @test epc.days        == 2458485
    @test epc.seconds     == 0
    @test epc.nanoseconds == 0

    epc += 1.0e-9
    @test epc.days        == 2458485
    @test epc.seconds     == 0
    @test epc.nanoseconds == 1

    epc -= 2.0e-9
    @test epc.days        == 2458484
    @test epc.seconds     == 86399
    @test epc.nanoseconds == 999999999

    for i in 1:365*86400
        epc += 1
    end
    @test epc.days        == 2458849
    @test epc.seconds     == 86399
    @test epc.nanoseconds == 999999999

    year, month, day, hour, minute, second, nanoseconds = caldate(epc)
    @test year        == 2020
    @test month       == 1
    @test day         == 1
    @test hour        == 11
    @test minute      == 59
    @test second      == 59
    @test nanoseconds == 999999999
end

# Test long term addition
let
    epc = Epoch(1980, 1, 1, 0, 0, 0, 1)
    year, month, day, hour, minute, second, nanoseconds = caldate(epc)
    @test year        == 1980
    @test month       == 1
    @test day         == 1
    @test hour        == 0
    @test minute      == 0
    @test second      == 0
    @test nanoseconds == 1
    @test epc.tsys    == "TAI"

    for i in 1:40*365.25*86400
        epc += 1
    end

    year, month, day, hour, minute, second, nanoseconds = caldate(epc)
    @test year        == 2020
    @test month       == 1
    @test day         == 1
    @test hour        == 0
    @test minute      == 0
    @test second      == 0
    @test nanoseconds == 1
    @test epc.tsys    == "TAI"
end

# Direct conversion
let
    epc = Epoch("2018-12-20T16:22:19.123456789Z")
    year, month, day, hour, minute, second, nanoseconds = caldate(epc)
    @test year        == 2018
    @test month       == 12
    @test day         == 20
    @test hour        == 16
    @test minute      == 22
    @test second      == 19
    @test nanoseconds == 123456789.0
end

# Conversion into different time system
let
    epc = Epoch("2018-12-20T16:22:19.123456789Z")
    year, month, day, hour, minute, second, nanoseconds = caldate(epc, tsys="TAI")
    @test year        == 2018
    @test month       == 12
    @test day         == 20
    @test hour        == 16
    @test minute      == 22
    @test second      == 19 + 37
    @test nanoseconds == 123456789.0

    epc = Epoch("2018-12-20T16:22:19.123456789Z")
    year, month, day, hour, minute, second, nanoseconds = caldate(epc, tsys="GPS")
    @test year        == 2018
    @test month       == 12
    @test day         == 20
    @test hour        == 16
    @test minute      == 22
    @test second      == 19 + 37 - 19
    @test nanoseconds == 123456789.0
end

let
    epc = Epoch(2000, 1, 1)
    @test jd(epc) == 2451544.5
    @test jd(epc, tsys="UTC") < 2451544.5
end

let
    epc = Epoch(2000, 1, 1, 0)
    @test mjd(epc) == 51544
    @test mjd(epc, tsys="UTC") < 51544
end

let
    epc = Epoch(2000, 1, 1)
    @test day_of_year(epc) == 1

    epc = Epoch(2000, 1, 1, 12)
    @test day_of_year(epc) == 1.5

    epc = Epoch(2000, 12, 31)
    @test day_of_year(epc) == 366

    epc = Epoch(2001, 1, 1)
    @test day_of_year(epc) == 1

    epc = Epoch(2001, 12, 31)
    @test day_of_year(epc) == 365
end

let
    epc = Epoch(2000, 1, 1)
    @test isapprox(gmst(epc, use_degrees=true), 99.835, atol=1e-3)

    epc = Epoch(2000, 1, 1)
    @test isapprox(gmst(epc, use_degrees=false), 1.742, atol=1e-3)
end

let
    epc = Epoch(2000, 1, 1)
    @test isapprox(gast(epc, use_degrees=true), 99.832, atol=1e-3)

    epc = Epoch(2000, 1, 1)
    @test isapprox(gast(epc, use_degrees=false), 1.742, atol=1e-3)
end

# Arithmetic Comparisons
let
    @test Epoch(2000, 1, 1, 12, 23, 59, 123456789) == Epoch(2000, 1, 1, 12, 23, 59, 123456789)
    @test Epoch(2000, 1, 1, 12, 23, 59, 123456789) != Epoch(2000, 1, 1, 12, 23, 59, 123456788)
end

let
    @test !(Epoch(2000, 1, 1, 12, 23, 59, 123456789) < Epoch(2000, 1, 1, 12, 23, 59, 123456789))
    @test !(Epoch(2000, 1, 1, 12, 23, 59, 123456789) <= Epoch(2000, 1, 1, 12, 23, 59, 123456788))
    @test Epoch(2000, 1, 1, 12, 23, 59, 123456788) < Epoch(2000, 1, 1, 12, 23, 59, 123456789)
    @test Epoch(2000, 1, 1, 12, 23, 59, 123456788) <= Epoch(2000, 1, 1, 12, 23, 59, 123456789)
end

let
    @test !(Epoch(2000, 1, 1, 12, 23, 59, 123456789) > Epoch(2000, 1, 1, 12, 23, 59, 123456789))
    @test !(Epoch(2000, 1, 1, 12, 23, 59, 123456788) >= Epoch(2000, 1, 1, 12, 23, 59, 123456789))
    @test Epoch(2000, 1, 1, 12, 23, 59, 123456789) > Epoch(2000, 1, 1, 12, 23, 59, 123456788)
    @test Epoch(2000, 1, 1, 12, 23, 59, 123456789) >= Epoch(2000, 1, 1, 12, 23, 59, 123456788)
end