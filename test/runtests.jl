using Test
using Random
using LinearAlgebra

using SatelliteDynamics

Random.seed!(0)

@testset "SatSchedule" begin
    testdir = joinpath(dirname(@__DIR__), "test")
    @testset "constants" begin
        include(joinpath(testdir, "test_constants.jl"))
    end
    @testset "universe" begin
        include(joinpath(testdir, "test_univ.jl"))
    end
    # @testset "time" begin
    #     include(joinpath(testdir, "test_time.jl"))
    # end
    # @testset "refsys" begin
    #     include(joinpath(testdir, "test_refsys.jl"))
    # end
    # @testset "coordinates" begin
    #     include(joinpath(testdir, "test_coordinates.jl"))
    # end
    # @testset "astrodynamics" begin
    #     include(joinpath(testdir, "test_astrodynamics.jl"))
    # end
    # @testset "orbit_dynamics" begin
    #     include(joinpath(testdir, "test_orbit_dynamics.jl"))
    # end
    # @testset "simulation" begin
    #     include(joinpath(testdir, "test_simulation.jl"))
    # end
    # @testset "data_structures" begin
    #     include(joinpath(testdir, "test_data_structures.jl"))
    # end
end