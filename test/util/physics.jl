using RvSpectMLBase
using Test

#=
 This only tests functions from src/utils that don't depend on data structures.
 test/data_structure.jl should test functions that depend on the data structure.
 Functions using multiple data structures should be tested in the more complex data structure's file.
=#

@testset "Simple Physics Functions" begin

    @testset "Util" begin
        @test RvSpectMLBase.calc_doppler_factor(10) ≈ 1.0000000333564094
        @test RvSpectMLBase.predict_intrinsic_stellar_line_width(5780) ≈ 9883.42046054907
    end

    @testset "calc_snr" begin
        @test RvSpectMLBase.calc_snr(8.0,4.0) ≈ 4
        @test RvSpectMLBase.calc_snr(8*ones(100),4*ones(100)) ≈ 40
    end

end
