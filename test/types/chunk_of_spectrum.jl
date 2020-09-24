using RvSpectMLBase
using Test

@testset "ChunkOfSpectrum" begin
        num_pixels = 100
        num_orders = 10
        @test_nowarn ChunkOfSpectrum(ones(num_pixels),ones(num_pixels),ones(num_pixels))

        lambda = reshape(range(5000.0, stop=6000, length=num_pixels*num_orders), num_pixels,num_orders)
        data = ones(100,10)
        order = 5
        pixels = 20:30
        @test_nowarn ChunkOfSpectrum(view(lambda,pixels,order), view(data,pixels,order), view(data,pixels,order))
        @test_nowarn ChunkOfSpectrum(lambda, data, data, order,pixels)
        @test_nowarn ChunkOfSpectrum(lambda, data, data, (pixels=pixels, order=order) )
end
