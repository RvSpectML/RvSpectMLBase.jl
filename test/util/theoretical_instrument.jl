import RvSpectMLBase
using Test

@testset "Theoretical Instrument" begin
    using RvSpectMLBase.TheoreticalInstrument
    num_pixels = 100
    num_orders = 10
    resolution = 1e5

    n1 = TheoreticalInstrument1D()
    lambda2d = collect(reshape(range(5000.0,stop=6000.0,length=num_pixels*num_orders),num_pixels,num_orders))
    @test_nowarn TheoreticalInstrument.TheoreticalInstrument2D(lambda2d)
    n2 = TheoreticalInstrument.TheoreticalInstrument2D(lambda2d)
    @testset "Calc λs" begin
            @test_nowarn TheoreticalInstrument.calc_λs(n1)
            @test_nowarn TheoreticalInstrument.calc_λs(n2)
    end
end
