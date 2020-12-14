"""
Code declaring ChunkOfSpectrum and abstract version.

Author: Eric Ford
Created: August 2020
"""

""" Abstract type for any ChunkOfSpectrum """
abstract type AbstractChunkOfSpectrum end

""" ChunkOfSpectrum for views into Spectra1DBasic or Spectra2DBasic """
struct ChunkOfSpectrum{T1<:Real,T2<:Real,T3<:Real,AA1<:AbstractArray{T1,1},AA2<:AbstractArray{T2,1},AA3<:AbstractArray{T3,1}} <: AbstractChunkOfSpectrum
    λ::AA1
    flux::AA2
    var::AA3
end

function ChunkOfSpectrum{T1,T2,T3}(λ::A1, flux::A2, var::A3) where {  T1<:Real, T2<:Real, T3<:Real, A1<:AbstractArray{T1,1}, A2<:AbstractArray{T2,1}, A3<:AbstractArray{T3,1} }
    @assert size(λ) == size(flux)
    @assert size(λ) == size(var)
    min_pixels_in_chunk = 3
    max_pixels_in_chunk = 10000
    #@assert min_pixels_in_chunk <= length(λ) <= max_pixels_in_chunk
    ChunkOfSpectrum{eltype(λ),eltype(flux),eltype(var),typeof(λ),typeof(flux),typeof(var)}(λ,flux,var)
end

function ChunkOfSpectrum(λ::A1, flux::A2, var::A3, order::Integer, pixels::AUR) where {  T1<:Real, T2<:Real, T3<:Real, A1<:AbstractArray{T1,2}, A2<:AbstractArray{T2,2}, A3<:AbstractArray{T3,2}, AUR<:AbstractUnitRange }
    @assert size(λ) == size(flux)
    @assert size(λ) == size(var)
    @assert 1 <= order <= size(λ,2)
    #println("order = ", order, "   pixels = ",pixels, "  size(λ,1) = ", size(λ,1))
    #flush(stdout)
    @assert 1 <= first(pixels) <= last(pixels) <= size(λ,1)
    ChunkOfSpectrum{T1,T2,T3}(view(λ,pixels,order),view(flux,pixels,order),view(var,pixels,order))
end

function ChunkOfSpectrum(λ::A1, flux::A2, var::A3, pixels::AUR) where {  T1<:Real, T2<:Real, T3<:Real, A1<:AbstractArray{T1,1}, A2<:AbstractArray{T2,1}, A3<:AbstractArray{T3,1}, AUR<:AbstractUnitRange }
    @assert size(λ) == size(flux)
    @assert size(λ) == size(var)
    @assert 1 <= first(pixels) <= last(pixels) <= size(λ,1)
    ChunkOfSpectrum{T1,T2,T3}(view(λ,pixels),view(flux,pixels),view(var,pixels))
end

function ChunkOfSpectrum(λ::A1, flux::A2, var::A3, loc::NamedTuple{(:pixels, :order),Tuple{AUR,I1}}) where {  T1<:Real, T2<:Real, T3<:Real, A1<:AbstractArray{T1,2}, A2<:AbstractArray{T2,2}, A3<:AbstractArray{T3,2}, AUR<:AbstractUnitRange, I1<:Integer }
    ChunkOfSpectrum(λ,flux,var,loc.order,loc.pixels)
end

function ChunkOfSpectrum(spectra::AS, order::Integer, pixels::AUR) where { AS<:AbstractSpectra2D, AUR<:AbstractUnitRange }
    ChunkOfSpectrum(spectra.λ,spectra.flux,spectra.var,order,pixels)
end

function ChunkOfSpectrum(spectra::AS, loc::NamedTuple{(:pixels, :order),Tuple{AUR,I1}}) where {  AS<:AbstractSpectra2D, AUR<:AbstractUnitRange, I1<:Integer }
    ChunkOfSpectrum(spectra.λ,spectra.flux,spectra.var,loc.order,loc.pixels)
end

function ChunkOfSpectrum(spectra::AS, pixels::AUR) where { AS<:AbstractSpectra1D, AUR<:AbstractUnitRange }
    ChunkOfSpectrum(spectra.λ,spectra.flux,spectra.var,pixels)
end

function ChunkOfSpectrum(spectra::AS, loc::NamedTuple{(:pixels, :order),Tuple{AUR,I1}}) where {  AS<:AbstractSpectra1D, AUR<:AbstractUnitRange, I1<:Integer }
    @assert loc.order == 1
    ChunkOfSpectrum(spectra.λ,spectra.flux,spectra.var,loc.pixels)
end

function empty_chunk_of_spectrum()
    ChunkOfSpectrum{Float64,Float64,Float64}(Float64[],Float64[],Float64[])
end

# TODO: Generalize in case not standard view's
get_order_index(chunk::AC) where { AC<:AbstractChunkOfSpectrum } = chunk.flux.indices[2]
get_pixels_range(chunk::AC) where { AC<:AbstractChunkOfSpectrum } = chunk.flux.indices[1]
