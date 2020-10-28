using RvSpectMLBase
using Test

#=
 This only tests functions from src/utils that don't depend on data structures.
 test/data_structure.jl should test functions that depend on the data structure.
 Functions using multiple data structures should be tested in the more complex data structure's file.
=#

@testset "File Utilities" begin
    @test_nowarn read_data_paths(paths_to_search=[pwd(),joinpath(pkgdir(RvSpectMLBase),"test")])
    @test_nowarn code_to_include_param_jl(paths_to_search=[pwd(),joinpath(pkgdir(RvSpectMLBase),"test")])
end
