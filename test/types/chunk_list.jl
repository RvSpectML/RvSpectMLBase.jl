using RvSpectMLBase
using Test


@testset "ChunkList" begin
    num_pixels = 100
    num_orders = 10
    resolution = 1e5
    lambda2d = collect(reshape(range(5000.0,stop=6000.0,length=num_pixels*num_orders),num_pixels,num_orders))
    inst = TheoreticalInstrument.TheoreticalInstrument2D(lambda2d)
    spectrum = TheoreticalInstrument.demo_generate_spectrum_line(inst)
    @test_nowarn make_orders_into_chunks(spectrum, inst)
    @test_nowarn make_orders_into_chunks(spectrum, orders_to_use=1:5, pixels_to_use=fill(20:80,5) )
end
