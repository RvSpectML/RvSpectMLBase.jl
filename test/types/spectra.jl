import RvSpectMLBase
using Test

@testset "Spectra1DBasic" begin
    using RvSpectMLBase.TheoreticalInstrument
    n1 = TheoreticalInstrument1D()
    lambdas = TheoreticalInstrument.calc_λs(n1)
    num_pixels = length(lambdas)

    @testset "Create Spectra1DBasic" begin
        @test_nowarn Spectra1DBasic(lambdas, ones(num_pixels), ones(num_pixels), n1)
    end
    spec = Spectra1DBasic(lambdas, ones(num_pixels), ones(num_pixels), n1)
end

@testset "Spectra2DBasic" begin
    num_pixels = 100
    num_orders = 10
    lambda2d = collect(reshape(range(5000.0,stop=6000.0,length=num_pixels*num_orders),num_pixels,num_orders))
    n2 = TheoreticalInstrument2D(lambda2d)
    @test_nowarn TheoreticalInstrument.demo_generate_spectrum_line(n2)
    spec = TheoreticalInstrument.demo_generate_spectrum_line(n2)
end
