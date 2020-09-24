"""
Code for creating and manipulating chunklists (i.e., list of views into a spectrum).
For example, creating a list of views into the orders of a spectrum to be analyzed or
creating a list of views into chunks around entries in a line list.

Author: Eric Ford
Created: August 2020
Contact: https://github.com/eford/
"""

""" make_orders_into_chunks

Return a ChunkList with a region of spectrum from each order in orders_to_use.
# Arguments
- spectra<:AbstractSpectra
- inst:  Instrument trait that provides default values
# Optional arguments
- orders_to_use: Range or Array (orders_to_use(inst))
- pixels_to_use: Array of Ranges (each from min_col to max_col)
or
- min_col: (min_col_default(inst,order)) and
- max_col: (max_col_default(inst,order))
"""
function make_orders_into_chunks
end

function make_orders_into_chunks(spectra::AS, inst::AbstractInstrument2D;
        orders_to_use = orders_to_use_default(spectra.inst) # 1:size(spectra.flux,2),
        ) where {AS<:AbstractSpectra }
    @assert eltype(orders_to_use) <: Integer
    @assert all( min_order(spectra.inst) .<= orders_to_use .<= max_order(spectra.inst) )
    pixels_to_use = map(ord->min_col_default(spectra.inst,ord):max_col_default(spectra.inst,ord),orders_to_use)
    make_orders_into_chunks(spectra,orders_to_use=orders_to_use, pixels_to_use=pixels_to_use)
end

function make_orders_into_chunks(spectra::AS;
        orders_to_use::Union{AR,AA1}, pixels_to_use::AAR ) where {
         AS<:AbstractSpectra2D, AR<:AbstractRange{Int64}, AA1<:AbstractArray{Int64,1}, AAR<:AbstractArray{AR,1} }
    @assert eltype(orders_to_use) <: Integer
    #@assert all( min_order(inst) .<= orders_to_use .<= max_order(inst) )
    #@assert minimum(pixels_to_use) >= min_pixel_in_order(inst)
    #@assert maximum(pixels_to_use) <= max_pixel_in_order(inst)
    ChunkList( map(order_idx->
                    ChunkOfSpectrum(spectra,(pixels=pixels_to_use[order_idx],order=orders_to_use[order_idx])),
                    1:length(orders_to_use) ), 1:length(orders_to_use) )
end

#=
# Doesn't really make sense for 1D spectra
function make_orders_into_chunks(spectra::AS;
        orders_to_use::Integer, pixels_to_use::AR ) where {
         AS<:AbstractSpectra, AR<:AbstractRange{Int64} }
    @assert eltype(orders_to_use) <: Integer
    #@assert all( min_order(inst) .<= orders_to_use .<= max_order(inst) )
    #@assert minimum(pixels_to_use) >= min_pixel_in_order(inst)
    #@assert maximum(pixels_to_use) <= max_pixel_in_order(inst)
    ChunkList( [ChunkOfSpectrum(spectra,(pixels=pixels_to_use,order=orders_to_use))], [orders_to_use] )
end
=#

""" make_grid_for_chunk
Create a range with equal spacing between points with end points set based on union of all chunks in timeseries.
# Arguments:
- timeseries: ChunkListTimeseries
- chunk index:
- oversample_factor: (1)
"""
function make_grid_for_chunk(timeseries::ACLT, c::Integer; oversample_factor::Real = 1.0, spacing::Symbol = :Linear, remove_rv_est::Bool = false ) where { ACLT<:AbstractChunkListTimeseries }
    num_obs = length(timeseries.chunk_list)
    @assert num_obs >= 1
    @assert 1<= c <= length(first(timeseries.chunk_list).data)
    @assert allequal(map(chunk->length(chunk.data),timeseries.chunk_list))
    @assert spacing == :Log || spacing == :Linear
    if spacing == :Log
        @warn "There's some issues with end points exceeding the bounds.  Round off error?  May cause bounds errors."
    end
    # Create grid, so that chunks at all times include the grid's minimum and maximum wavelength.
    if remove_rv_est   @assert haskey(first(timeseries.metadata),:rv_est)   end
    boost_factor = [ remove_rv_est ? calc_doppler_factor(timeseries.metadata[t][:rv_est]) : 1 for t in 1:num_obs ]
    λ_min = maximum(minimum(timeseries.chunk_list[t].data[c].λ)/boost_factor[t] for t in 1:num_obs)
    λ_max = minimum(maximum(timeseries.chunk_list[t].data[c].λ)/boost_factor[t] for t in 1:num_obs)
    Δλ_grid_obs = mean(log(timeseries.chunk_list[t].data[c].λ[end]/
                           timeseries.chunk_list[t].data[c].λ[1]   )/
                         (length(timeseries.chunk_list[t].data[c].λ)-1) for t in 1:num_obs)
    num_pixels_obs = log(λ_max/λ_min)/Δλ_grid_obs
    num_pixels_gen = (num_pixels_obs-1) * oversample_factor + 1
    if spacing == :Log
        Δlnλ_grid_obs = mean(log(timeseries.chunk_list[t].data[c].λ[end]/
                                 timeseries.chunk_list[t].data[c].λ[1]   )/
                                 (length(timeseries.chunk_list[t].data[c].λ)-1) for t in 1:num_obs)
        num_pixels_obs = log(λ_max/λ_min)/Δlnλ_grid_obs
        num_pixels_gen = (num_pixels_obs-1) * oversample_factor + 1
            Δlnλ_grid_gen = log(λ_max/λ_min)/ (num_pixels_gen-1)
        return exp.(range(log(λ_min),stop=log(λ_max),step=Δlnλ_grid_gen))
    elseif spacing == :Linear
        Δλ_grid_obs = mean((timeseries.chunk_list[t].data[c].λ[end]-timeseries.chunk_list[t].data[c].λ[1] )/
                            (length(timeseries.chunk_list[t].data[c].λ)-1) for t in 1:num_obs)
        num_pixels_obs = (λ_max-λ_min)/Δλ_grid_obs
        num_pixels_gen = (num_pixels_obs-1) * oversample_factor + 1
        Δλ_grid_gen = (λ_max-λ_min)/ (num_pixels_gen-1)
        return range(λ_min,stop=λ_max,step=Δλ_grid_gen)
    end
end

""" Return a ChunkList of best regions of spectrum with lines in line_line.
    line_list is a DataFrame containing :lambda_lo and :lambda_hi.
    Pads edges by Δ.
"""
function make_chunk_list(spectra::AS, line_list::DataFrame; rv_shift::Real = 0, Δ::Real=Δλoλ_edge_pad_default) where { AS<:AbstractSpectra }
    @assert hasproperty(line_list,:lambda_lo)
    @assert hasproperty(line_list,:lambda_hi)
    boost_factor = calc_doppler_factor(rv_shift)
    if rv_shift != 0
        @warn("I haven't tested this yet, especially the sign.")  # TODO
    end
    line_locs = map(row->find_line_best(row.lambda_lo*boost_factor,row.lambda_hi*boost_factor,spectra,Δ=Δ), eachrow(line_list) )
    cl = ChunkList(map(loc->ChunkOfSpectrum(spectra,loc),line_locs), map(loc->loc[2], line_locs) )
end

function make_chunk_list_timeseries(spectra::AS,chunk_list_df::DataFrame; rv_shift::Real = 0) where {ST<:AbstractSpectra, AS<:AbstractArray{ST,1} }
    times = map(s->s.metadata[:bjd],spectra)
    metadata = make_vec_metadata_from_spectral_timeseries(spectra)
    time_series_of_chunk_lists = map(spec->RvSpectMLBase.make_chunk_list(spec,chunk_list_df, rv_shift=rv_shift),spectra)
    chunk_list_timeseries = ChunkListTimeseries(times, time_series_of_chunk_lists, inst=first(spectra).inst, metadata=metadata )
end

function extract_chunk_list_timeseries_for_order(clt::AbstractChunkListTimeseries, order::Integer) # where {ST<:AbstractSpectra, AS<:AbstractArray{ST,1} }
    @assert min_order(clt.inst) <= order <= max_order(clt.inst)
    num_chunks_to_search = num_chunks(clt)
    chunk_order = map(ch->first(clt.chunk_list).data[ch].λ.indices[2], 1:num_chunks_to_search )
    chunks_in_order = findall(chunk_order .== order)
    num_chunks_in_order = sum(chunks_in_order)
    if !(num_chunks_in_order >= 1)
        return nothing
    end
    new_chunk_list_timeseries = Vector{ChunkList}(undef, length(clt.times))
    for i in 1:length(clt.times)
        chunk_list = ChunkList(clt.chunk_list[i].data[chunks_in_order], fill(order,num_chunks_in_order) )
        new_chunk_list_timeseries[i] = chunk_list
    end
    output = ChunkListTimeseries(clt.times, new_chunk_list_timeseries, inst=clt.inst, metadata=clt.metadata )
    return output
end

#=
function make_order_list_timeseries(spectra::AS) #= , order_list::AOL ) =# where {ST<:AbstractSpectra, AS<:AbstractArray{ST,1} #=, CLT<:AbstractChunkList, AOL::AbstractArray{CLT,1} =# }
    times = map(s->s.metadata[:bjd],spectra)
    inst = first(spectra).inst
    metadata = make_vec_metadata_from_spectral_timeseries(spectra)
    order_list = map( spec->RvSpectMLBase.make_orders_into_chunks(spec,inst), spectra)
    chunk_list_timeseries = ChunkListTimeseries(times, order_list, inst=first(spectra).inst, metadata=metadata )
end
=#

function make_order_list_timeseries(spectra::AS; orders_to_use = orders_to_use_default(first(spectra).inst) ) #= , order_list::AOL ) =# where {ST<:AbstractSpectra, AS<:AbstractArray{ST,1} #=, CLT<:AbstractChunkList, AOL::AbstractArray{CLT,1} =# }
    times = map(s->s.metadata[:bjd],spectra)
    inst = first(spectra).inst
    metadata = make_vec_metadata_from_spectral_timeseries(spectra)
    order_list = map( spec->RvSpectMLBase.make_orders_into_chunks(spec,inst, orders_to_use=orders_to_use), spectra)
    chunk_list_timeseries = ChunkListTimeseries(times, order_list, inst=first(spectra).inst, metadata=metadata )
end


function make_chunk_list_expr(spectra::AS, range_list::DataFrame ) where { AS<:AbstractSpectra }
    @assert hasproperty(range_list,:lambda_lo)
    @assert hasproperty(range_list,:lambda_hi)
    cl = typeof(ChunkOfSpectrum(spectra,1,1:size(spectra.flux,1)))[]
    ord_list = Int64[]
    for row in eachrow(range_list)
        orders = find_orders_in_range(row.lambda_lo,row.lambda_hi, spectra.λ)
        if length(orders) == 1
            order = orders[1]
            pixels = find_cols_to_fit(spectra.λ[:,order],row.lambda_lo,row.lambda_hi)
            if 1 <= first(pixels) <= size(spectra.λ,1) && 1 <= last(pixels) <= size(spectra.λ,1) &&  # valid location
                calc_snr(spectra.flux[pixels,order],spectra.var[pixels,order]) > 0     # not just NaN's
                println("# Found one order (", order, ") for range ", row.lambda_lo, " - ", row.lambda_hi )
                push!(cl, ChunkOfSpectrum(spectra,order,pixels) )
                push!(ord_list, order)
            else
                println("# Found one order (", order, ") for range ", row.lambda_lo, " - ", row.lambda_hi, " but not valid or all NaNs.")
                continue
            end
        elseif length(orders) > 1
            pixels = map(ord->find_cols_to_fit(spectra.λ[:,ord],row.lambda_lo,row.lambda_hi), orders)
            scores = map( i->calc_snr(spectra.flux[pixels[i],orders[i]],spectra.var[pixels[i],orders[i]]), 1:length(orders) )
            idx_best = findmax(scores)[2]
            println("# Found ", length(orders), " orders for range ", row.lambda_lo, " - ", row.lambda_hi, " picking ", orders[idx_best])
            push!(cl, ChunkOfSpectrum(spectra,orders[idx_best],pixels[idx_best]) )
            push!(ord_list, orders[idx_best])
            idx_all = findall(scores)

        else # length(orders) == 0
            #println("# Found no orders for range ", row.lambda_lo, " - ", row.lambda_hi )

        end
    end
    return ChunkList(cl,ord_list)
end


""" Return (chunk_timeseries, line_list) that have been trimmed of any chunks that are bad based on any spectra in the chunk_timeseries.
    chunk_timeseries: ChunkListTimeseries
    line_linst:  DataFrame w/ lambda_lo, lambda_hi
    verbose: print debugging info (false)
"""
#= Is there any reason to keep this version?
function filter_bad_chunks(chunk_list_timeseries::ACLT, line_list::DataFrame; verbose::Union{Int,Bool} = false) where { ACLT<:AbstractChunkListTimeseries }
    @assert(length(chunk_list_timeseries)>=1)
    @assert(hasproperty(line_list,:lambda_lo))
    @assert(hasproperty(line_list,:lambda_hi))
    inst = chunk_list_timeseries.inst
    idx_keep = trues(num_chunks(chunk_list_timeseries))
    for t in 1:length(chunk_list_timeseries)
        idx_bad_λ = findall(c->any(isnan.(chunk_list_timeseries.chunk_list[t].data[c].λ)),1:num_chunks(chunk_list_timeseries))
        idx_bad_flux = findall(c->any(isnan.(chunk_list_timeseries.chunk_list[t].data[c].flux)),1:num_chunks(chunk_list_timeseries))
        idx_bad_var = findall(c->any(isnan.(chunk_list_timeseries.chunk_list[t].data[c].var)),1:num_chunks(chunk_list_timeseries))
        idx_not_sorted = findall(c->!issorted(chunk_list_timeseries.chunk_list[t].data[c].λ),1:num_chunks(chunk_list_timeseries))
        idx_keep[idx_bad_λ] .= false
        idx_keep[idx_bad_flux] .= false
        idx_keep[idx_bad_var] .= false
        idx_keep[idx_not_sorted] .= false
        if verbose && (length(idx_bad_λ)+length(idx_bad_flux)+length(idx_bad_var)+length(idx_not_sorted))
               flush(stdout)
               println("# Removing chunks", vcat(idx_bad_λ,idx_bad_flux,idx_bad_var), " at time ", t, " due to NaNs (",
                            length(idx_bad_λ),",",length(idx_bad_flux),",",length(idx_bad_var),") and ",
                            length(idx_not_sorted), " for not being sorted.")
        end
        #=
        # TODO: Move these kinds of checks to a traits/plan-based system
        if hasproperty(chunk_list_timeseries.metadata[t],:pixel_mask)
            idx_keep[.!chunk_list_timeseries.metadata[t][:pixel_mask]] .= false
        end
        if hasproperty(chunk_list_timeseries.metadata[t],:excalibur_mask)
            idx_keep[.!chunk_list_timeseries.metadata[t][:excalibur_mask]] .= false
        end
        =#
    end
    chunks_to_remove = findall(.!idx_keep)
    if length(chunks_to_remove) == 0
        println("# No lines to remove.")
        return (chunk_timeseries=chunk_list_timeseries, line_list=line_list)
    else
        println("# Removing ", length(chunks_to_remove), " chunks.")
        map(c->println("# ",c,": ",line_list.lambda_lo[c]," - ",line_list.lambda_hi[c]),chunks_to_remove)
        new_line_list = line_list[findall(idx_keep),:]
        new_chunk_list_timeseries = [ChunkList(chunk_list_timeseries.chunk_list[t].data[idx_keep]) for t in 1:length(chunk_list_timeseries) ]
        return (chunk_timeseries=ChunkListTimeseries(chunk_list_timeseries.times,new_chunk_list_timeseries, inst=inst, metadata=chunk_list_timeseries.metadata), line_list=new_line_list)
    end
end
=#
function filter_bad_chunks(chunk_list_timeseries::ACLT; verbose::Bool = false) where { ACLT<:AbstractChunkListTimeseries }
    @assert(length(chunk_list_timeseries)>=1)
    idx_keep = trues(num_chunks(chunk_list_timeseries))
    for t in 1:length(chunk_list_timeseries)
        idx_bad_λ = findall(c->any(isnan.(chunk_list_timeseries.chunk_list[t].data[c].λ)),1:num_chunks(chunk_list_timeseries))
        idx_bad_flux = findall(c->any(isnan.(chunk_list_timeseries.chunk_list[t].data[c].flux)),1:num_chunks(chunk_list_timeseries))
        idx_bad_var = findall(c->any(isnan.(chunk_list_timeseries.chunk_list[t].data[c].var)),1:num_chunks(chunk_list_timeseries))
        idx_keep[idx_bad_λ] .= false
        idx_keep[idx_bad_flux] .= false
        idx_keep[idx_bad_var] .= false
        if verbose && (length(idx_bad_λ)+length(idx_bad_flux)+length(idx_bad_var) > 0)
               println("# Removing chunks", vcat(idx_bad_λ,idx_bad_flux,idx_bad_var), " at time ", t, " due to NaNs (",length(idx_bad_λ),",",length(idx_bad_flux),",",length(idx_bad_var),").")
        end
    end
    chunks_to_remove = findall(.!idx_keep)
    if length(chunks_to_remove) == 0
        if verbose   println("# Nothing to remove.")   end
        return chunk_list_timeseries
    else
        if verbose
            println("# Removing ", length(chunks_to_remove), " chunks due to NaNs.")
            map(c->println("# ",c,": ",findall(.!idx_keep)))
        end
        new_chunk_list_timeseries = [ChunkList(chunk_list_timeseries.chunk_list[t].data[idx_keep], chunk_list_timeseries.chunk_list[t].order[idx_keep]) for t in 1:length(chunk_list_timeseries) ]
        #return ChunkListTimeseries(chunk_list_timeseries.times[idx_keep],new_chunk_list_timeseries, inst=chunk_list_timeseries.inst, metadata=chunk_list_timeseries.metadata[idx_keep])
        return ChunkListTimeseries(chunk_list_timeseries.times,new_chunk_list_timeseries, inst=chunk_list_timeseries.inst, metadata=chunk_list_timeseries.metadata)
    end
end


""" Calc normalization of chunk based on average flux in a ChunkOfSpectrum. """
function calc_normalization(chunk::AC) where { AC<:AbstractChunkOfSpectrum}
    total_flux = NaNMath.sum(Float64.(chunk.flux))
    num_pixels = length(chunk.flux)
    scale_fac = num_pixels / total_flux
end

""" Calc normalization of chunk based on average flux in a ChunkOfSpectrum using inverse variance weighting. """
function calc_normalization_var_weighted(chunk::AC) where { AC<:AbstractChunkOfSpectrum}
    sum_weighted_flux = NaNMath.sum(Float64.(chunk.flux) ./ chunk.var )
    sum_weights = NaNMath.sum(1.0 ./ chunk.var)
    scale_fac = sum_weights / sum_weighted_flux
end

""" Calc normalization of spectra based on average flux in a ChunkList. """
function calc_normalization(chunk_list::ACL) where { ACL<:AbstractChunkList}
    #total_flux = NaNMath.sum(NaNMath.sum(Float64.(chunk_list.data[c].flux))
    #                       for c in 1:length(chunk_list) )
    total_flux = sum(sum(Float64.(chunk_list.data[c].flux))
                        for c in 1:length(chunk_list) )
    num_pixels = sum( length(chunk_list.data[c].flux) for c in 1:length(chunk_list) )
    scale_fac = num_pixels / total_flux
end

""" Calc normalization of spectra based on average flux in a ChunkList using inverse variance weighting. """
function calc_normalization_var_weighted(chunk_list::ACL) where { ACL<:AbstractChunkList}
    #total_flux = NaNMath.sum(NaNMath.sum(Float64.(chunk_list.data[c].flux))
    #                    for c in 1:length(chunk_list) )
    total_flux = sum(sum(Float64.(chunk_list.data[c].flux))
                    for c in 1:length(chunk_list) )
    num_pixels = sum( length(chunk_list.data[c].flux) for c in 1:length(chunk_list) )
    scale_fac = num_pixels / total_flux
end
