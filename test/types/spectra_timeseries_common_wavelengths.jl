import RvSpectMLBase
using Test

#=
@testset "Spectra Timeseries Common Wavelengths" begin
    # TODO Write test for Spectra Timeseries Common Wavelengths type
    # Won't actually test it here, since code to make this from a timeseries is not in Base package
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
    clt = spectra  # Is this right type?
    chunk_map = map(i->make_grid_for_chunk(clt,i,oversample_factor=oversample_factor), 1:num_chunks(clt))
    # need to figure out dimensions for arguments below
    @test_nowarn stscw = SpectralTimeSeriesCommonWavelengths(λ, flux, var, chunk_map, TheoreticalInstrument1D() )
end
=#
