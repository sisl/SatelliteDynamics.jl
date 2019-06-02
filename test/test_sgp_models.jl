# Test TLE checksum
let
    line1 = "1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753"
    line2 = "2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667"

    csum1 = SGPModels.tle_checksum(line1[1:end-1])
    csum2 = SGPModels.tle_checksum(line2[1:end-1])

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
    line1 = "1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753"
    line2 = "2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667"

    tle = TLE(line1, line2)

    println(tle.line1)
    println(tle.line2)

    println(tle.epoch)
    # println(sgp4(tle, tle.epoch))
    # println(state(tle, tle.epoch))
    for i in 0:12
        println(360.0*i)
        println(state(tle, tle.epoch+360.0*i))
        println()
    end
end