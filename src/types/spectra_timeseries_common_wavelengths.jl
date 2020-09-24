"""
Code declaring SpectralTimeSeriesCommonWavelengths and their abstract versions.

Author: Eric Ford
Created: August 2020
"""

""" Abstract type for a time series of spectra that share a common wavelength grid. """
abstract type AbstractSpectralTimeSeriesCommonWavelengths <: AbstractSpectra1D   end

""" Time series of spectra that share a common wavelength grid. """
struct SpectralTimeSeriesCommonWavelengths{T1<:Real,T2<:Real,T3<:Real,AA1<:AbstractArray{T1,1},AA2<:AbstractArray{T2,2},AA3<:AbstractArray{T3,2},
            AA4<:AbstractArray{UnitRange{Int64},1}, InstT<:AbstractInstrument } <: AbstractSpectralTimeSeriesCommonWavelengths
    λ::AA1
    flux::AA2
    var::AA3
    chunk_map::AA4
    inst::InstT
    metadata::MetadataT
end

function SpectralTimeSeriesCommonWavelengths(λ::A1, flux::A2, var::A3, chunk_map::A4, inst::InstT;
        metadata::MetadataT = MetadataT() ) where {
          T1<:Real, T2<:Real, T3<:Real, A1<:AbstractArray{T1,1}, A2<:AbstractArray{T2,2}, A3<:AbstractArray{T3,2},
          A4<:AbstractArray{UnitRange{Int64},1}, InstT<:AbstractInstrument }
    #println("len(λ) = ", length(λ))
    #println("size(flux) = ", size(flux))
    @assert length(λ) == size(flux,1)
    @assert length(λ) == size(var,1)
    @assert 1 <= length(λ)
    SpectralTimeSeriesCommonWavelengths{eltype(λ),eltype(flux),eltype(var),typeof(λ),typeof(flux),typeof(var),typeof(chunk_map),typeof(inst)}(λ,flux,var,chunk_map,inst,metadata)
end

"""   make_spectral_time_series_common_wavelengths_with_selected_times( input, time_idx )
Returns a SpectralTimeSeriesCommonWavelengths, retaining only those times and spectra specified by time_idx.
"""
function make_spectral_time_series_common_wavelengths_with_selected_times(input::STSCWT, time_idx::AbstractVector{T1} ) where { STSCWT<:AbstractSpectralTimeSeriesCommonWavelengths, T1<:Integer }
    @assert minimum(time_idx) >= 1
    @assert maximum(time_idx) <= size(input.flux,2)
    metadata = length(input.metadata) == size(input.flux,2) ? input.metadata[time_idx] : MetadataT()
    SpectralTimeSeriesCommonWavelengths(input.λ, input.flux[:,time_idx], input.var[:,time_idx], input.chunk_map, input.inst, metadata=metadata )
end
