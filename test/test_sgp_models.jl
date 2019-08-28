# Test TLE checksum
let
    line1 = "1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753"
    line2 = "2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667"

    csum1 = tle_checksum(line1[1:end-1])
    csum2 = tle_checksum(line2[1:end-1])

    @test csum1 == parse(Int, line1[end])
    @test csum2 == parse(Int, line2[end])
end

# Test TLE initialization
let
    # Bad line 1 checksum
    line1 = "1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4752"
    line2 = "2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667"

    @test_throws ArgumentError TLE(line1, line2)

    # Bad line 2 checksum
    line1 = "1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753"
    line2 = "2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413668"

    @test_throws ArgumentError TLE(line1, line2)

    # Bad line 1 length
    line1 = "1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  759"
    line2 = "2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667"

    @test_throws ArgumentError TLE(line1, line2)

    # Bad line 1 length
    line1 = "1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753"
    line2 = "2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.8241915741361"

    @test_throws ArgumentError TLE(line1, line2)
end

let
    # Test Spacecraft Report 3 TLE 1
    line1 = "1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753"
    line2 = "2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667"

    tle = TLE(line1, line2)
    rv  = sgp4(tle, tle.epoch)
    array_isapprox(rv[1:3], [7022.46647249, -1400.06656182, 0.05106558], atol=1e-8)
    # array_isapprox(rv[4:6], [1.893831081, 6.405894873, 4.534806701], atol=1e-8)

    # Test Spacecraft Report 3 TLE 2
    line1 = "1 04632U 70093B   04031.91070959 -.00000084  00000-0  10000-3 0  9955"
    line2 = "2 04632  11.4628 273.1101 1450506 207.6000 143.9350  1.20231981 44145"

    tle = TLE(line1, line2)
    rv  = sgp4(tle, tle.epoch)
    # array_isapprox(rv[1:3], [2334.12240168, -41920.42940755, -0.03695311], atol=1e-8)
    # array_isapprox(rv[4:6], [2.826320327, -0.065091320, 0.570935915], atol=1e-8)

    # Test Spacecraft Report 3 TLE 3
    line1 = "1 06251U 62025E   06176.82412014  .00008885  00000-0  12808-3 0  3985"
    line2 = "2 06251  58.0579  54.0425 0030035 139.1568 221.1854 15.56387291  6774"

    tle = TLE(line1, line2)
    rv  = sgp4(tle, tle.epoch)
    array_isapprox(rv[1:3], [3988.29488121, 5498.97562026, 0.92898335], atol=1e-8)
    # array_isapprox(rv[4:6], [-3.290042995, 2.357636871, 6.496621836], atol=1e-8)
end

let
    # Test Spacecraft Report 3 TLE 1
    line1 = "1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753"
    line2 = "2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667"

    tle = TLE(line1, line2)
    rv  = ecef(tle, tle.epoch)
end

let
    # Test Spacecraft Report 3 TLE 1
    line1 = "1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753"
    line2 = "2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667"

    tle = TLE(line1, line2)
    rv  = eci(tle, tle.epoch)
end