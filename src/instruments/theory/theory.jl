"""
   Delegates loading functions & traits for a generic theoretical spectrograph

Author: Eric Ford and collaborators
Created: August 2020
"""

"""
Module providing types and traits and customized functions for a generic theoretical spectrograph.
"""
module TheoreticalInstrument
using ..RvSpectMLBase
#import ..RvSpectML: AbstractInstrument, AbstractInstrument1D, AbstractInstrument2D
using DataFrames # , FITSIO


default_theoretical_instrument_resolution = 137500
default_theoretical_instrument_λ_min = 4700.0
default_theoretical_instrument_λ_max = 6000.0

""" Trait for 1D spectra from a theoretical instrument """
struct TheoreticalInstrument1D <: AbstractInstrument1D
    λ_min::Float64
    λ_max::Float64
    resolution::Float64

    function TheoreticalInstrument1D(λl::Real, λh::Real, res::Real)
        @assert 3000 <= λl < λh
        @assert λl < λh <= 10000
        @assert 10000 <= res <= 1000000
        new(λl,λh,res)
    end
end

function TheoreticalInstrument1D(;λ_min::Real = default_theoretical_instrument_λ_min, λ_max::Real = default_theoretical_instrument_λ_max, resolution::Real = default_theoretical_instrument_resolution )
    TheoreticalInstrument1D(λ_min,λ_max,resolution)
end


λ_min(x::TheoreticalInstrument1D) = x.λ_min
λ_max(x::TheoreticalInstrument1D) = x.λ_max

#=  TODO: Write once figure out 1d case first
struct TheoreticalInstrument2D <: AbstractInstrument2D
end
=#

""" Trait for a 2D spectra from a theoretical instrument """
struct TheoreticalInstrument2D <: AbstractInstrument2D
    λ_min::AbstractArray{Float64,1}
    λ_max::AbstractArray{Float64,1}
    pixels_per_order::Int64

    function TheoreticalInstrument2D(λ_in::AbstractArray{Float64,2}; pixels_per_order::Integer = size(λ_in,1))
        @assert size(λ_in,1) >= 6
        @assert size(λ_in,2) >= 1
        @assert 1 <= pixels_per_order <= 16386

        lambda_min = vec(minimum(λ_in,dims=1))
        lambda_max = vec(maximum(λ_in,dims=1))
        new(lambda_min, lambda_max, pixels_per_order)
    end
end

# Trait for any spectra from TheoreticalInstrument (could improve by using SimpleTraits)
const AnyTheoreticalInstrument = Union{TheoreticalInstrument1D,TheoreticalInstrument2D}
export TheoreticalInstrument, TheoreticalInstrument1D, TheoreticalInstrument2D, AnyTheoreticalInstrument

λ_min(x::TheoreticalInstrument2D) = minimum(x.λ_min)
λ_max(x::TheoreticalInstrument2D) = maximum(x.λ_max)

include("traits.jl")
export min_pixel, max_pixel, min_order, max_order, min_pixel_in_order, max_pixel_in_order
export orders_to_use_default, min_col_default, max_col_default
export metadata_symbols_default, metadata_strings_default
export default_ccf_mask_v_width

# Need to think about what to do for IO, if anything.
#include("io.jl")
#export make_manifest, read_metadata, read_data, read_solar_data

include("util.jl")

include("gen_data.jl")
export generate_spectrum, generate_spectra_timeseries

end
