""" Delegates loading of code specifying types essential to the package and other packages that use it.  """

include("instruments.jl")
export AbstractInstrument, AbstractInstrument1D, AbstractInstrument2D
#export Generic1D, Generic2D   # WARNING: Might remove these in future

include("spectra.jl")
export AbstractSpectra, AbstractSpectra1D, AbstractSpectra2D
export Spectra1DBasic, Spectra2DBasic

include("chunks.jl")
export AbstractChunkOfSpectrum, ChunkOfSpectrum

include("chunk_list.jl")
export AbstractChunkList, ChunkList
export AbstractChunkListTimeseries, ChunkListTimeseries
export length, num_chunks, num_times
export extract_chunklist_timeseries_with_subset_obs

include("spectra_timeseries_common_wavelengths.jl")
export AbstractSpectralTimeSeriesCommonWavelengths
export SpectralTimeSeriesCommonWavelengths
#=
# WARNING: Still experimental
#export make_vec_metadata_from_spectral_timeseries
export make_spectral_time_series_common_wavelengths_with_selected_times
=#
