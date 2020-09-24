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
    df = DataFrame(:lambda=>[5501], :weight=>[0.5])
    @test_nowarn TheoreticalInstrument.generate_spectra_timeseries(times,df,inst, rvs)
    spectra = TheoreticalInstrument.generate_spectra_timeseries(times,df,inst, rvs)
end
