using RvSpectMLBase
using Test

@testset "Utility functions that don't use package's types." begin
    include("physics.jl")
    include("cs_algs.jl")
    include("files.jl")
end

@testset "Utility functions for RvSpecMLBase's types" begin
    include("chunk_of_spectrum.jl")
    include("theoretical_instrument.jl")
    #include("find_line.jl")                                # TODO: Add tests
    include("chunk_list.jl")
    include("spectra.jl")
    include("spectra_timeseries.jl")
    #include("spectra_timeseries_common_wavelengths.jl")    # TODO: Add tests
end
