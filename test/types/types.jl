using RvSpectMLBase
using Test

@testset "Constructors and bare-bones functions" begin
    include("chunk_of_spectrum.jl")
    include("theoretical_instrument.jl")
    include("spectra.jl")
    include("chunk_list.jl")
    include("spectra_timeseries.jl")    
    #include("spectra_timeseries_common_wavelengths.jl")  # TODO: Add tests
end
