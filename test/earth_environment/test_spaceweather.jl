
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
    @test f107Observed(53941) == 75.9
    @test f107Adjusted(53941) == 78.3
    @test isapprox(f107ObservedAvg(53941), 72.2951, atol=1e-4)
    @test isapprox(f107AdjustedAvg(53941), 74.5593, atol=1e-4)
    
    data = f107Data(53941)
    @test data[1] == 75.9
    @test data[2] == 78.3
    @test isapprox(data[3], 72.2951, atol=1e-4)
    @test isapprox(data[4], 74.5593, atol=1e-4)
end