using RvSpectMLBase
using Test

#=
 This only tests functions from src/utils that don't depend on data structures.
 test/data_structure.jl should test functions that depend on the data structure.
 Functions using multiple data structures should be tested in the more complex data structure's file.
=#

@testset "Computer Science Algorithms" begin

    @testset "allequal" begin
    @test RvSpectMLBase.allequal([1,1,1,1])
    @test !RvSpectMLBase.allequal([1,2,3,4])
    end

    @testset "findargminmax" begin
        res = RvSpectMLBase.findargminmax([10.0, 30.0, -20.0, 40.0])
        @test res.min == -20.0
        @test res.max ==  40.0
        @test res.argmin ==  3
        @test res.argmax == 4
    end
    @testset "searchsortednearest" begin
        @test RvSpectMLBase.searchsortednearest(-cos.(π*(1:16)./16), 0.0) == 8
    end
end
