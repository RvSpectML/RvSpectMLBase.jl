import RvSpectMLBase
using Test

@testset "Spectra1DBasic" begin
    using RvSpectMLBase.TheoreticalInstrument
    n1 = TheoreticalInstrument1D()
    lambdas = TheoreticalInstrument.calc_λs(n1)
    num_pixels = length(lambdas)
    spec = Spectra1DBasic(lambdas, ones(num_pixels), ones(num_pixels), n1)
end

@testset "Spectra2DBasic" begin
    num_pixels = 1000
    num_orders = 10
    lambda2d = collect(reshape(range(5000.0,stop=6000.0,length=num_pixels*num_orders),num_pixels,num_orders))
    n2 = TheoreticalInstrument2D(lambda2d)
    spec = TheoreticalInstrument.demo_generate_spectrum_line(n2)
    @test get_λ_range(spec).min ≈ 5000
    @test get_λ_range(spec).max ≈ 6000

    df = DataFrame(:lambda_lo=>[5501], :lambda_hi=>[5505])
    cl = make_chunk_list_around_lines(spec, df)
    @test_nowarn calc_normalization(cl)
    norm = calc_normalization(cl)
    @test_nowarn normalize_spectrum!(spec,norm )
    @test calc_normalization(cl) ≈ 1
    @test_nowarn calc_normalization_var_weighted(cl)
    norm = calc_normalization(cl)
    @test_nowarn normalize_spectrum!(spec,norm )
    @test calc_normalization_var_weighted(cl) ≈ 1
    @test_nowarn apply_doppler_boost!(spec,1.001)
    @test_nowarn RvSpectMLBase.discard_large_metadata(spec)
end
