using RvSpectMLBase
using Test

@testset "ChunkOfSpectrum" begin
        num_pixels = 100
        num_orders = 10
        lambda = reshape(range(5000.0, stop=6000, length=num_pixels*num_orders), num_pixels,num_orders)
        data = ones(num_pixels,num_orders)
        order = 5
        pixels = 20:30
        @test_nowarn ChunkOfSpectrum(lambda, data, data, (pixels=pixels, order=order) )
        ch = ChunkOfSpectrum(lambda, data, data, (pixels=pixels, order=order) )
        @test RvSpectMLBase.calc_normalization(ch) ≈ 1
        @test RvSpectMLBase.calc_normalization_var_weighted(ch) ≈ 1
        @test RvSpectMLBase.find_orders_with_line(5421.0,lambda)[1] == order
        @test RvSpectMLBase.find_pixels_for_line_in_chunk(ch, 5420, 5422) == 2:3
end
