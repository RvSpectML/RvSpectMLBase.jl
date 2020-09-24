import RvSpectMLBase
using Test

@testset "Spectra Timeseries" begin
    using RvSpectMLBase.TheoreticalInstrument
    using DataFrames
    num_pixels = 1000
    num_orders = 10
    num_obs = 3
    lambda2d = collect(reshape(range(5000.0,stop=6000.0,length=num_pixels*num_orders),num_pixels,num_orders))
    inst = TheoreticalInstrument2D(lambda2d)
    times = sort(rand(num_obs))
    rvs = range(0.0, stop=10.0, length=length(times))
    line_list = DataFrame(:lambda=>[5501], :weight=>[0.5])
    @test_nowarn TheoreticalInstrument.generate_spectra_timeseries(times,line_list,inst, rvs)
    spectra = TheoreticalInstrument.generate_spectra_timeseries(times,line_list,inst, rvs)
    chunk_list = DataFrame(:lambda_lo=>[5501], :lambda_hi=>[5505])
    @test_nowarn make_chunk_list_timeseries(spectra,chunk_list)
    clt = make_chunk_list_timeseries(spectra,chunk_list)
    @test num_times(clt) == num_obs
    @test num_chunks(clt) == 1
    @test_nowarn normalize_spectra!(clt,spectra )
    @test all(map(cl->calc_normalization(cl),clt.chunk_list) .≈ 1 )
    @test_nowarn make_order_list_timeseries(spectra)
    @test_nowarn filter_bad_chunks(clt)
end
