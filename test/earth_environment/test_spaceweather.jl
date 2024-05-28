

let
    load_space_weather_data(TEST_SW_DATA)
end

##########################
# Geomagnetic Index Data #
##########################

let
    @test length(KpIndices(53941)) == 8
    @test isapprox(KpDailyIndex(53941), 12.3333, atol=1e-4)
    @test length(ApIndices(53941)) == 8
    @test isapprox(ApDailyIndex(53941), 6.0, atol=1e-4)
end

###################
# Solar Flux Data #
###################

let
    @test f107Observed(53941) == 75.5
    @test f107Adjusted(53941) == 77.9
    @test isapprox(f107ObservedAvg(53941), 77.3, atol=1e-4)
    @test isapprox(f107AdjustedAvg(53941), 79.6, atol=1e-4)
    
    data = f107Data(53941)
    @test data[1] == 75.5
    @test data[2] == 77.9
    @test isapprox(data[3], 77.3, atol=1e-4)
    @test isapprox(data[4], 79.6, atol=1e-4)
end