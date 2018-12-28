using Test
using Random
using LinearAlgebra
using Logging

using SatelliteDynamics

# Set logging level
global_logger(SimpleLogger(stderr, Logging.Debug))

# Fix randomness during tests
Random.seed!(0)

@time @testset "SatelliteDynamics Package Tests" begin
    testdir = joinpath(dirname(@__DIR__), "test")
    @time @testset "SatelliteDynamics.Contants" begin
        include(joinpath(testdir, "test_constants.jl"))
    end
    @time @testset "SatelliteDynamics.Universe" begin
        include(joinpath(testdir, "test_univ.jl"))
    end
    @time @testset "SatelliteDynamics.Time" begin
        include(joinpath(testdir, "test_time.jl"))
    end
    @time @testset "SatelliteDynamics.Refsys" begin
        include(joinpath(testdir, "test_reference_systems.jl"))
    end
    @time @testset "SatelliteDynamics.Attitude" begin
        include(joinpath(testdir, "test_attitude.jl"))
    end
    @time @testset "SatelliteDynamics.Coordinates" begin
        include(joinpath(testdir, "test_coordinates.jl"))
    end
    @time @testset "SatelliteDynamics.Astrodynamics" begin
        include(joinpath(testdir, "test_astrodynamics.jl"))
    end
    @time @testset "SatelliteDynamics.OrbitDynamics" begin
        include(joinpath(testdir, "test_orbitdynamics.jl"))
    end
    @time @testset "SatelliteDynamics.Simulation" begin
        include(joinpath(testdir, "test_simulation.jl"))
    end
end