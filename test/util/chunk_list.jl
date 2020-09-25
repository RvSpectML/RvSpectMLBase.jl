using RvSpectMLBase
using Test

@testset "ChunkList" begin
    using DataFrames
    num_pixels = 1000
    num_orders = 10
    lambda2d = collect(reshape(range(5000.0,stop=6000.0,length=num_pixels*num_orders),num_pixels,num_orders))
    inst = TheoreticalInstrument2D(lambda2d)
    spectrum = TheoreticalInstrument.demo_generate_spectrum_line(inst)
    @test_nowarn make_orders_into_chunks(spectrum, inst)
    @test_nowarn make_orders_into_chunks(spectrum, orders_to_use=1:5, pixels_to_use=fill(20:80,5) )

    df = DataFrame(:lambda_lo=>[5501], :lambda_hi=>[5505])
    @test_nowarn make_chunk_list(spectrum, df)
    cl = make_chunk_list(spectrum, df)

    num_obs = 3
    times = sort(rand(num_obs))
    @test_skip ChunkListTimeseries(times, [cl, cl, cl])
end
