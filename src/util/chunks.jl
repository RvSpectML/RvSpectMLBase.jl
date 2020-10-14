"""
Code for creating and manipulating chunklists (i.e., list of views into a spectrum).
For example, creating a list of views into the orders of a spectrum to be analyzed or
creating a list of views into chunks around entries in a line list.

Author: Eric Ford
Created: August 2020
Contact: https://github.com/eford/
"""

function make_chunk_list_from_loc_df(spectra::AS, inst::AbstractInstrument2D, df_λ_good::DataFrame ) where {AS<:AbstractSpectra }
    chunk_locs_df = make_order_pixel_list_for_chunks_in_order(spectra, inst, df_λ_good)
    ChunkList( map(r-> ChunkOfSpectrum(spectra,(pixels=r.pixels,order=r.order)),  eachrow(chunk_locs_df) ),
            map(r->r.order,eachrow(chunk_locs_df) ) )

end

function make_chunk_list_timeseries_telluric_free(spectra::AS, df_λ_good::DataFrame; min_pixels_in_chunk::Integer = 250, verbose::Bool = false ) where {ST<:AbstractSpectra, AS<:AbstractArray{ST,1} }
    max_Δv_match_chunk = 2*max_bc
    inst = first(spectra).inst
    times = map(s->s.metadata[:bjd],spectra)
    metadata = make_vec_metadata_from_spectral_timeseries(spectra)
    time_series_of_chunk_lists = map(spec->RvSpectMLBase.make_chunk_list_from_loc_df(spec, inst, df_λ_good), spectra)
    @assert all(map(cl->length(cl.data),time_series_of_chunk_lists) .>= 1)
    (min_num_chunks, obs_idx_min_num_chunks) = findmin(map(obs_idx->length(time_series_of_chunk_lists[obs_idx].data),1:length(time_series_of_chunk_lists) ))
    obs_idx_min_num_chunks = findall(map(obs_idx->length(time_series_of_chunk_lists[obs_idx].data)==min_num_chunks,1:length(time_series_of_chunk_lists) ))

    λ_min_targets = map(ch_idx->median(map(obs_idx->minimum(time_series_of_chunk_lists[obs_idx][ch_idx].λ),obs_idx_min_num_chunks)),1:min_num_chunks)
    λ_max_targets = map(ch_idx->median(map(obs_idx->maximum(time_series_of_chunk_lists[obs_idx][ch_idx].λ),obs_idx_min_num_chunks)),1:min_num_chunks)
    order_idx_targets = map(ch_idx->round(Int,median(map(obs_idx->time_series_of_chunk_lists[obs_idx][ch_idx].λ.indices[2],obs_idx_min_num_chunks))),
                                1:min_num_chunks)
    similarity_of_chunks = zeros(min_num_chunks)
    similarity_chunk_idx = zeros(Int64,min_num_chunks)
    for obs_idx in 1:length(time_series_of_chunk_lists)
        #=
        if length(time_series_of_chunk_lists[obs_idx].data) == min_num_chunks
            continue  # TODO: Technically should still check that they all match
        end
        =#
        ch_perm_idx = zeros(Int64,min_num_chunks)
        for ch_idx in 1:min_num_chunks
            if verbose    println("# Trying to match ch_idx= ",ch_idx, " lmin= ", λ_min_targets[ch_idx], " lmax= ", λ_max_targets[ch_idx])   end
            idx_order_match = findall(map(ch->time_series_of_chunk_lists[obs_idx][ch].λ.indices[2], 1:length(time_series_of_chunk_lists[obs_idx].data)).==order_idx_targets[ch_idx])
            if verbose  &&  ch_idx >= 110
                println("# idx_order_match= ", idx_order_match)
            end

            if !(length(idx_order_match) >= 1)   continue    end
            # Find if any ch_idx2 is good match to ch_idx's wavelength range
            similarity_of_chunks .= Inf
            similarity_chunk_idx .= 0
            for (i,ch_idx2) in enumerate(idx_order_match)
                delta_λmax = maximum(time_series_of_chunk_lists[obs_idx][ch_idx2].λ)-λ_max_targets[ch_idx]
                delta_λmin = minimum(time_series_of_chunk_lists[obs_idx][ch_idx2].λ)-λ_min_targets[ch_idx]
                if verbose    && ch_idx >= 110
                    println("# ch_idx2 = ", ch_idx2, " dlmin= ", delta_λmin/λ_max_targets[ch_idx]*speed_of_light_mps, " dlmax= ", delta_λmax/λ_max_targets[ch_idx]*speed_of_light_mps)
                end

                similarity_of_chunks[i] = ((abs(delta_λmax)/λ_max_targets[ch_idx]*speed_of_light_mps) / max_Δv_match_chunk)^2 +
                                          ((abs(delta_λmin)/λ_max_targets[ch_idx]*speed_of_light_mps) / max_Δv_match_chunk)^2
                similarity_chunk_idx[i] = ch_idx2
                #= if (abs(delta_λmax)/λ_max_targets[ch_idx]*speed_of_light_mps < max_Δv_match_chunk) &&
                    (abs(delta_λmax)/λ_max_targets[ch_idx]*speed_of_light_mps < max_Δv_match_chunk)
                    ch_perm_idx[ch_idx] = ch_idx2
                end =#
            end # for over ch_idx2 to match target

            ch_idx2 = similarity_chunk_idx[findmin(similarity_of_chunks)[2]]
            delta_λmax = maximum(time_series_of_chunk_lists[obs_idx][ch_idx2].λ)-λ_max_targets[ch_idx]
            delta_λmin = minimum(time_series_of_chunk_lists[obs_idx][ch_idx2].λ)-λ_min_targets[ch_idx]
            if (abs(delta_λmax)/λ_max_targets[ch_idx]*speed_of_light_mps < max_Δv_match_chunk) &&
                (abs(delta_λmin)/λ_max_targets[ch_idx]*speed_of_light_mps < max_Δv_match_chunk)
                ch_perm_idx[ch_idx] = ch_idx2
            #= else
                flush(stdout)
                println("# obs_idx= ", obs_idx, " ch_idx= ",ch_idx, " ch_idx2= ",ch_idx2, " delta_λmax=",delta_λmax, " delta_λmin=",delta_λmin, " similarity_of_chunks=",similarity_of_chunks[ch_idx2])
                flush(stdout)
                =#
            end
        end # for over ch_idx that is target
        if verbose    println("# obs_idx = ", obs_idx)   end
        if verbose    println("# ch_perm_idx = ",ch_perm_idx)   end
        @assert !any(ch_perm_idx.==0)
        if !all(ch_perm_idx.==1:length(ch_perm_idx))
            time_series_of_chunk_lists[obs_idx] = ChunkList( map(ch_idx-> time_series_of_chunk_lists[obs_idx].data[ch_idx],  ch_perm_idx),
                                                             map(ch_idx-> time_series_of_chunk_lists[obs_idx].order[ch_idx], ch_perm_idx) )
        end
    end # for over obs_idx
    min_pixels_in_chunk_in_any_obs = map(ch_idx->minimum(map(obs_idx->length(time_series_of_chunk_lists[obs_idx][ch_idx].λ),1:length(time_series_of_chunk_lists) )),1:length(first(time_series_of_chunk_lists)) )
    idx_keep = findall(x->x>min_pixels_in_chunk,min_pixels_in_chunk_in_any_obs)

    time_series_of_nonsmall_chunk_lists = map(obs_idx -> ChunkList( map(ch_idx-> time_series_of_chunk_lists[obs_idx].data[ch_idx],  idx_keep),
                                                                    map(ch_idx-> time_series_of_chunk_lists[obs_idx].order[ch_idx], idx_keep) ), 1:length(time_series_of_chunk_lists) )

    chunk_list_timeseries = ChunkListTimeseries(times, time_series_of_nonsmall_chunk_lists, inst=inst, metadata=metadata )
    #chunk_list_timeseries = ChunkListTimeseries(times, time_series_of_chunk_lists, inst=inst, metadata=metadata )
end


function make_order_pixel_list_for_chunks_in_order(spectra::AS, inst::AbstractInstrument2D, df_λ_good::DataFrame;
            orders_to_use = orders_to_use_default(spectra.inst), verbose::Bool = false
            ) where {AS<:AbstractSpectra }
    @assert eltype(orders_to_use) <: Integer
    @assert all( min_order(spectra.inst) .<= orders_to_use .<= max_order(spectra.inst) )
    @assert hasproperty(df_λ_good, :lambda_lo)
    @assert hasproperty(df_λ_good, :lambda_hi)
    min_pixels_per_chunk = 5  # arbitrary for now

    df_out = DataFrame(:pixels=>UnitRange[], :order=>Int[])
    for order_idx in orders_to_use
        pixels_to_use = min_col_default(spectra.inst,order_idx):max_col_default(spectra.inst,order_idx)
        if length(pixels_to_use) < min_pixels_per_chunk
            if verbose   println("# Skipping order_idx= ", order_idx, " due to pixel length")   end
            continue
        end
        if verbose    println("# Order idx= ",order_idx)   end
        # TODO OPT: could speed up by scanning through wavelengths and list together) or avoiding this allocation all together
        idx_good = is_in_wavelength_range_list.(spectra.λ[pixels_to_use,order_idx], list=df_λ_good)

        idx_start = findfirst(idx_good)
        idx_stop = findlast(idx_good)
        if isnothing(idx_start) || isnothing(idx_stop)
            if verbose   println("# Skipping order_idx= ", order_idx, " due to finding no good pixels")   end
            continue
        end
        idx_start = idx_start
        idx_stop = idx_stop
        if verbose   println("# Searching pixels:", pixels_to_use[idx_start]:pixels_to_use[idx_stop], " or λ= ",spectra.λ[pixels_to_use[idx_start],order_idx], " - ", spectra.λ[pixels_to_use[idx_stop],order_idx])   end
        idx_lo = idx_start
        idx_hi = idx_start
        while idx_hi < idx_stop
            idx_lo = findfirst( view(idx_good,idx_hi:idx_stop) )
            if isnothing(idx_lo)
                if verbose   println("# Searching pixels:", pixels_to_use[idx_hi]:pixels_to_use[idx_stop], " or λ= ",spectra.λ[pixels_to_use[idx_start],order_idx], " - ", spectra.λ[pixels_to_use[idx_stop],order_idx])   end
                break
            end
            idx_lo += idx_hi-1
            idx_hi = findfirst( view(.!idx_good,idx_lo:idx_stop) )
            if isnothing(idx_hi)
                idx_hi = idx_stop
            else
                idx_hi += idx_lo-2
            end
            if idx_hi-idx_lo > min_pixels_per_chunk # TODO Make min_pixels_per_chunk
                if verbose   println("# Found pixels:", pixels_to_use[idx_lo]:pixels_to_use[idx_hi], " or λ= ",spectra.λ[pixels_to_use[idx_lo],order_idx], " - ", spectra.λ[pixels_to_use[idx_hi],order_idx])   end
                push!(df_out,Dict(:pixels=>pixels_to_use[idx_lo]:pixels_to_use[idx_hi], :order=>order_idx))
                idx_hi += 1
            end
        end # while still in order
    end  # for over order_to_use
    return df_out
end

""" `make_orders_into_chunks`
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
                    1:length(orders_to_use) ), orders_to_use )
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

""" `make_grid_for_chunk`
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

""" `make_chunk_list_around_lines( spectra, line_list)`
Return a ChunkList of best regions of spectrum with lines in line_line.
    line_list is a DataFrame containing :lambda_lo and :lambda_hi.
    Pads edges by Δ.
"""
function make_chunk_list_around_lines(spectra::AS, line_list::DataFrame; rv_shift::Real = 0, Δ::Real=Δλoλ_edge_pad_default) where { AS<:AbstractSpectra }
    @assert hasproperty(line_list,:lambda_lo)
    @assert hasproperty(line_list,:lambda_hi)
    boost_factor = calc_doppler_factor(rv_shift)
    if rv_shift != 0
        @warn("I haven't tested this yet, especially the sign.")  # TODO
    end
    line_locs = map(row->find_line_best(row.lambda_lo*boost_factor,row.lambda_hi*boost_factor,spectra,Δ=Δ), eachrow(line_list) )
    # TODO: Update to add order to ChunkList
    cl = ChunkList(map(loc->ChunkOfSpectrum(spectra,loc),line_locs), map(loc->loc[2], line_locs) )
end

function make_chunk_list_timeseries_around_lines(spectra::AS,line_list_df::DataFrame; rv_shift::Real = 0) where {ST<:AbstractSpectra, AS<:AbstractArray{ST,1} }
    times = map(s->s.metadata[:bjd],spectra)
    metadata = make_vec_metadata_from_spectral_timeseries(spectra)
    time_series_of_chunk_lists = map(spec->RvSpectMLBase.make_chunk_list_around_lines(spec,line_list_df, rv_shift=rv_shift), spectra)
    chunk_list_timeseries = ChunkListTimeseries(times, time_series_of_chunk_lists, inst=first(spectra).inst, metadata=metadata )
end

function extract_chunk_list_timeseries_for_order(clt::AbstractChunkListTimeseries, order::Integer) # where {ST<:AbstractSpectra, AS<:AbstractArray{ST,1} }
    # TODO: Decide/document if using order or order_idx for selecting chunks
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
    # TODO: Decide/document if using order or order_idx for selecting chunks
    times = map(s->s.metadata[:bjd],spectra)
    inst = first(spectra).inst
    metadata = make_vec_metadata_from_spectral_timeseries(spectra)
    order_list = map( spec->RvSpectMLBase.make_orders_into_chunks(spec,inst, orders_to_use=orders_to_use), spectra)
    chunk_list_timeseries = ChunkListTimeseries(times, order_list, inst=first(spectra).inst, metadata=metadata )
end


function make_chunk_list_expr(spectra::AS, range_list::DataFrame ) where { AS<:AbstractSpectra }
    # TODO: Decide/document if using order or order_idx for selecting chunks
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


""" Return (chunk_timeseries) that have been trimmed of any chunks that are bad based on any spectra in the chunk_timeseries.
For now just checks for NaNs.  Instruments can provide their own checks.
Inputs:
- `chunk_timeseries`: ChunkListTimeseries
Optional arguemnts:
- `verbose`: print debugging info (false)
Returns:
- ChunkListTimeseries that has removed problematic chunks.
"""
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



#= This version has some EXPRES-specific code.  Could incorporate into to EchelleInstruments.EXPRES module.
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
