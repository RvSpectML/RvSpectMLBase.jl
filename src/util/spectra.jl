"""
Various utilities for manipulating spectra

Author: Eric Ford
Created: August 2020
Contact: https://github.com/eford/
"""

import NaNMath

"""
    `apply_doppler_boost!(spectrum, doppler_factor)`
    `apply_doppler_boost!(spectra, df)`

Apply Doppler boost to spectra's λ's and update its metadata[:doppler_factor], so it will know how to undo the transform.
# Arguments:
- `spectrum::AbstractSpectra`: spectrum to be boosted
- `doppler_factor`: boost factor (1 = noop)

or:
- `spectra`: spectra to be boosted
- `df`: DataFrame provides `:drift` and `:ssb_rv` (in m/s) for calculating the Doppler boost for each spectrum

Returns spectrum/spectra with λ boosted
TODO: Improve documentation formatting.  This can serve as a template.
"""
function apply_doppler_boost! end

function apply_doppler_boost!(spectra::AS,doppler_factor::Real) where {AS<:AbstractSpectra}
    if doppler_factor == one(doppler_factor) return spectra end
    #println("# t= ",time, " doppler_factor= ",doppler_factor)
    spectra.λ .*= doppler_factor
    if hasproperty(spectra.metadata,:doppler_factor)
        spectra.metadata[:doppler_factor] /= doppler_factor
    else
        spectra.metadata[:doppler_factor] = 1/doppler_factor
    end
    return spectra
end

function apply_doppler_boost!(spectra::AbstractArray{AS}, df::DataFrame ) where { AS<:AbstractSpectra }
    @assert size(spectra,1) == size(df,1)
    local doppler_factor = ones(size(spectra))
    if !hasproperty(df,:drift) @info "apply_doppler_boost! didn't find :drift to apply."   end
    if  hasproperty(df,:drift)        doppler_factor .*= calc_doppler_factor.(df[!,:drift])          end
    if !hasproperty(df,:drift) @info "apply_doppler_boost! didn't find :ssb_rv to apply."  end
    if  hasproperty(df,:ssb_rv)       doppler_factor   .*= calc_doppler_factor.(df[!,:ssb_rv])       end
    if !hasproperty(df,:drift) @info "apply_doppler_boost! didn't find :diff_ext_rv to apply."  end
    if  hasproperty(df,:diff_ext_rv)  doppler_factor   .*= calc_doppler_factor.(df[!,:diff_ext_rv])  end
    map(x->apply_doppler_boost!(x[1],x[2]), zip(spectra,doppler_factor) );
end

"""
`calc_snr(flux, var)`
`calc_snr(spectrum, pixels, order_idx)`

Calculate total SNR in (region of) spectra.
"""
function calc_snr end

function calc_snr(flux::AbstractArray{T1},var::AbstractArray{T2}) where {T1<:Real, T2<:Real}
    @assert size(flux) == size(var)
    sqrt(NaNMath.sum( flux.^2 ./ var))
end

function calc_snr(flux::Real,var::Real)
    flux / sqrt(var)
end

function calc_snr(spectrum::ST, pixels::AR, order::Integer) where { ST<:AbstractSpectra2D, AR<:AbstractRange{Int64}, AA1<:AbstractArray{Int64,1} } #, AAR<:AbstractArray{AR,1} }
    flux = view(spectrum.flux,pixels,order)
    var = view(spectrum.var,pixels,order)
    #=
    #if all(isnan.(flux) .|| isnan.(var) )
    if all( isnan.( flux.^2 ./ var ) )
        println("# all NaNs in order = ", order, " pixels = ", pixels)
        return 0
    end
    =#
    calc_snr( view(spectrum.flux,pixels,order), view(spectrum.var,pixels,order) )
end

""" Normalize spectrum, multiplying fluxes by scale_fac. """
function normalize_spectrum!(spectrum::ST, scale_fac::Real) where { ST<:AbstractSpectra }
    @assert 0 < scale_fac < Inf
    @assert !isnan(scale_fac^2)
    spectrum.flux .*= scale_fac
    spectrum.var .*= scale_fac^2
    return spectrum
end


""" Normalize each spectrum based on sum of fluxes in chunk_timeseries region of each spectrum. """
function normalize_spectra!(chunk_timeseries::ACLT, spectra::AS) where { ACLT<:AbstractChunkListTimeseries, ST<:AbstractSpectra, AS<:AbstractArray{ST} }
    @assert length(chunk_timeseries) == length(spectra)
    for t in 1:length(chunk_timeseries)
        scale_fac = calc_normalization(chunk_timeseries.chunk_list[t])
        # println("# t= ",t, " scale_fac= ", scale_fac)
        normalize_spectrum!(spectra[t], scale_fac)
    end
    return chunk_timeseries
end

""" Return instrument associated with spectrum """
function get_inst(spectrum::AS) where { AS<:AbstractSpectra2D }
    return spectrum.inst
end

""" Return instrument associated with first spectrum in array """
function get_inst(spectra::AAS) where { AS<:AbstractSpectra2D, AAS<:AbstractVector{AS} }
    @assert length(spectra)>=1
    get_inst(first(spectra))
end

""" Return the largest minimum wavelength and smallest maximum wavelength of a spectrum.
"""
function get_λ_range(data::CLT) where { CLT<:AbstractSpectra }
   λmin = minimum(data.λ)
   λmax = maximum(data.λ)
   return (min = λmin, max = λmax)
end

""" Return the largest minimum wavelength and smallest maximum wavelength across an array of spectra.
Calls get_λ_range(AbstractSpectra2D) that should be specialized for each instrument. """
function get_λ_range(data::ACLT) where { CLT<:AbstractSpectra, ACLT<:AbstractArray{CLT} }
   λminmax = get_λ_range.(data)
   λmin = maximum(map(p->p[1],λminmax))
   λmax = minimum(map(p->p[2],λminmax))
   return (min = λmin, max = λmax)
end

""" Extract the metadata from a time series of spectra and return it as an array. """
function make_vec_metadata_from_spectral_timeseries(spec_arr::AA) where { AS<:AbstractSpectra, AA<:AbstractArray{AS,1} }
    map(s->s.metadata,spec_arr)
end

""" Return DataFrame with information about which pixels and wavelengths to use from each order"""
function get_order_info(all_spectra::AAS; orders_to_use::Union{UnitRange,AbstractVector{T1} } = min_order(get_inst(all_spectra)):max_order(get_inst(all_spectra))) where { AS<:AbstractSpectra2D, AAS<:AbstractVector{AS}, T1<:Integer }
  @assert length(orders_to_use) >= 1

  inst = get_inst(all_spectra)
  num_orders = length(orders_to_use)
  #pixel_ranges = Array{UnitRange,1}(undef,num_orders)
  order_info = Any[]
  for order in orders_to_use #min_order(inst):max_order(inst)
    this_order_info = Dict{Symbol,Any}(:order=>order)
    pixel_range = get_pixel_range(inst,order)
    min_pixels_in_order = 60
    #=
    pixel_range = min_col_default(inst,order):max_col_default(inst,order)
    if typeof(inst) <: AnyEXPRES &&
        ( EXPRES.max_col_excalibur(inst, order) - EXPRES.min_col_excalibur(inst, order) >= min_pixels_in_order )

        minc = max(first(pixel_range), EXPRES.min_col_excalibur(inst,order), EXPRES.min_col_nonnan(inst,order))
        maxc = min(last(pixel_range), EXPRES.max_col_excalibur(inst,order), EXPRES.max_col_nonnan(inst,order))
        pixel_range = minc:maxc
      end
    =#
    if length(pixel_range) > min_pixels_in_order
        this_order_info[:pixel_range] = pixel_range
        (this_order_info[:λ_min], this_order_info[:λ_max], this_order_info[:mean_Δv]) = get_shared_wavelength_range_for_order(all_spectra,order, pixels_to_use=pixel_range )
        push!(order_info,this_order_info)
    end
  end
  return DataFrame(order_info)
end

function calc_mean_Δv(λ::AbstractVector{T1}, λ_min::Real, λ_max::Real ) where { T1<:Real }
    idx_min = searchsortedfirst(λ,λ_min)
    idx_max = searchsortedlast(λ,λ_max)
    @assert 1 <= idx_min <= length(λ)
    @assert 1 <= idx_max <= length(λ)
    Δv = log(λ[idx_max]/λ[idx_min])/(idx_max-idx_min) * RvSpectMLBase.speed_of_light_mps
end

function get_shared_wavelength_range_for_order(spectra::AbstractVector{AS}, order::Integer; pixels_to_use::AR = min_col_default(get_inst(spectra),order):max_col_default(get_inst(spectra),order),
             boost_factor::AA1 = ones(length(spectra)), verbose::Bool=false ) where {
             AS<:AbstractSpectra2D, AR<:AbstractRange{Int64}, T1<:Real, AA1<:AbstractVector{T1} }
        num_obs = length(spectra)
        @assert 1 <= num_obs # <= max_num_spectra
        if verbose   println("# order = ", order, " pixels_to_use = ", pixels_to_use)   end
        λ_min = maximum(minimum(spectra[t].λ[pixels_to_use,order])/boost_factor[t] for t in 1:num_obs)
        λ_max = minimum(maximum(spectra[t].λ[pixels_to_use,order])/boost_factor[t] for t in 1:num_obs)
        mean_Δv = sum(map(t->calc_mean_Δv(view(spectra[t].λ, pixels_to_use, order), λ_min, λ_max),1:num_obs)) / num_obs
        return (λ_min=λ_min, λ_max=λ_max, mean_Δv=mean_Δv)
        #=
        @assert spacing == :Log || spacing == :Linear
        if spacing == :Log
            @warn "There's some issues with end points exceeding the bounds.  Round off error?  May cause bounds errors."
        end
        =#
end

function discard_large_metadata(data::Union{ACLT,AS, AAS}) where { ACLT<:AbstractChunkListTimeseries, AS<:AbstractSpectra, AAS<:AbstractArray{AS} }
    if typeof(data.inst) <: AnyTheoreticalInstrument
        # Nothing to do
    else
        inst_module = get_inst_module(data.inst)
        inst_module.discard_large_metadata(data)
    end

    #=
    if typeof(data.inst) <: AnyEXPRES
        inst_module.discard_blaze(data)
        inst_module.discard_continuum(data)
        inst_module.discard_tellurics(data)
        inst_module.discard_pixel_mask(data)
        inst_module.discard_excalibur_mask(data)
    end
    if typeof(data.inst) <: AnyNEID
        #discard_blaze(data)
        #discard_continuum(data)
        #discard_tellurics(data)
        #discard_pixel_mask(data)
        #discard_excalibur_mask(data)
    end
    =#
end

function discard_blaze(metadata::Dict{Symbol,Any} )
    delete!(metadata,:blaze)
end

function discard_blaze(data::AST) where { AST<:AbstractSpectra }
   discard_blaze(data.metadata)
end

function discard_blaze(data::ACLT) where { CLT<:AbstractSpectra, ACLT<:AbstractArray{CLT} }
   map(spectra->discard_blaze(spectra.metadata),data)
end

function discard_blaze(data::CLT) where { CLT<:AbstractChunkListTimeseries }
   map(discard_blaze,data.metadata)
end

function discard_continuum(metadata::Dict{Symbol,Any} )
    delete!(metadata,:continuum)
end

function discard_continuum(data::AST) where { AST<:AbstractSpectra }
   discard_continuum(data.metadata)
end

function discard_continuum(data::ACLT) where { CLT<:AbstractSpectra, ACLT<:AbstractArray{CLT} }
   map(spectra->discard_continuum(spectra.metadata),data)
end

function discard_continuum(data::CLT) where { CLT<:AbstractChunkListTimeseries }
   map(discard_continuum,data.metadata)
end

function discard_tellurics(metadata::Dict{Symbol,Any} )
    delete!(metadata,:tellurics)
end

function discard_tellurics(data::AST) where { AST<:AbstractSpectra }
   discard_tellurics(data.metadata)
end

function discard_tellurics(data::ACLT) where { CLT<:AbstractSpectra, ACLT<:AbstractArray{CLT} }
   map(spectra->discard_tellurics(spectra.metadata),data)
end

function discard_tellurics(data::CLT) where { CLT<:AbstractChunkListTimeseries }
   map(discard_tellurics,data.metadata)
end

function discard_pixel_mask(metadata::Dict{Symbol,Any} )
    delete!(metadata,:pixel_mask)
end

function discard_pixel_mask(data::AST) where { AST<:AbstractSpectra }
   discard_pixel_mask(data.metadata)
end

function discard_pixel_mask(data::ACLT) where { CLT<:AbstractSpectra, ACLT<:AbstractArray{CLT} }
   map(spectra->discard_pixel_mask(spectra.metadata),data)
end

function discard_pixel_mask(data::CLT) where { CLT<:AbstractChunkListTimeseries }
   map(discard_pixel_mask,data.metadata)
end

function discard_excalibur_mask(metadata::Dict{Symbol,Any} )
    delete!(metadata,:excalibur_mask)
end

function discard_excalibur_mask(data::AST) where { AST<:AbstractSpectra }
   discard_excalibur_mask(data.metadata)
end

function discard_excalibur_mask(data::ACLT) where { CLT<:AbstractSpectra, ACLT<:AbstractArray{CLT} }
   map(spectra->discard_excalibur_mask(spectra.metadata),data)
end

function discard_excalibur_mask(data::CLT) where { CLT<:AbstractChunkListTimeseries }
   map(discard_excalibur_mask,data.metadata)
end
