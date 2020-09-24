import RvSpectMLBase
using Test

@testset "Theoretical Instrument" begin
    using RvSpectMLBase.TheoreticalInstrument
    num_pixels = 100
    num_orders = 10
    #resolution = 1e5

    @test_nowarn TheoreticalInstrument1D()
    n1 = TheoreticalInstrument1D()

    @testset "1D Extracted Spectra Traits" begin
            @test typeof(n1) <: AbstractInstrument
            @test typeof(n1) <: AnyTheoreticalInstrument
            @test min_pixel(n1) == 1
            @test max_pixel(n1) == 128*8192
            @test first(pixels_all(n1)) == min_pixel(n1)
            @test last(pixels_all(n1)) == max_pixel(n1)
            @test default_ccf_mask_v_width(n1) > 0
            @test get_inst_module(n1) == TheoreticalInstrument
    end

    lambda2d = collect(reshape(range(5000.0,stop=6000.0,length=num_pixels*num_orders),num_pixels,num_orders))
    @test_nowarn TheoreticalInstrument2D(lambda2d)
    n2 = TheoreticalInstrument2D(lambda2d)

    @testset "2D Extracted Spectra Traits" begin
            @test typeof(n2) <: AbstractInstrument
            @test typeof(n2) <: AnyTheoreticalInstrument
            @test min_order(n2) == 1
            @test max_order(n2) == num_orders
            @test first(orders_all(n2)) == min_order(n2)
            @test last(orders_all(n2)) == max_order(n2)
            @test orders_to_use_default(n2) == min_order(n2):max_order(n2)
            @test min_pixel_in_order(n2) == 1
            @test max_pixel_in_order(n2) == num_pixels
            @test first(pixels_all(n2)) == min_pixel_in_order(n2)
            @test last(pixels_all(n2)) == max_pixel_in_order(n2)
            @test default_ccf_mask_v_width(n2) > 0
            @test get_inst_module(n2) == TheoreticalInstrument
    end
end
