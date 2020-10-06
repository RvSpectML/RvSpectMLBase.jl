"""
Code for finding a line or wavelength range in spectra.

Author: Eric Ford
Created: August 2020
Contact: https://github.com/eford/
"""

""" Return list of all order indices that contain a pixel with wavelength lambda """
function find_orders_with_line(goal::Real,lambda::AbstractArray{T,2}) where T<:Real
   order_min(i) = lambda[1,i]
   order_max(i) = lambda[end,i]
#=
    for i in 1:10
       println("# i= ",i," order_min= ",order_min(i)," order_max= ",order_max(i), "   goal= ",goal)
   end
   flush(stdout)
=#
   findall(i->order_min(i)<=goal<=order_max(i), 1:size(lambda,2) )
end


""" Return list of all order indices that contain all pixels with wavelengths between goal_lo and goal_hi """
function find_orders_with_line(goal_lo::Real,goal_hi::Real,lambda::AbstractArray{T,2}) where T<:Real
   order_min(i) = lambda[1,i]
   order_max(i) = lambda[end,i]
   findall(i->order_min(i)<=goal_lo && goal_hi<=order_max(i), 1:size(lambda,2) )
end

""" Return list of all order indices that include any wavelengths between goal_lo and goal_hi """
function find_orders_in_range(goal_lo::Real,goal_hi::Real,lambda::AbstractArray{T,2}) where T<:Real
   order_min(i) = lambda[1,i]
   order_max(i) = lambda[end,i]
   findall(i-> (goal_lo<=order_min(i)<=goal_hi) || (goal_lo<=order_max(i)<=goal_hi), 1:size(lambda,2) )
end


# Find indicies for pixels around lines
const Δλoλ_fit_line_default = 5*(1.8*1000/speed_of_light_mps)
const Δλoλ_edge_pad_default = 0*(1.8*1000/speed_of_light_mps)

""" Return a range of columns indices with wavelengths within Δ of line_center """
function find_cols_to_fit(wavelengths::AbstractArray{T,1}, line_center::Real; Δ::Real = Δλoλ_fit_line_default) where T<:Real
    @assert Δ >= zero(Δ)
    first = findfirst(x->x>=line_center*(1-Δ),wavelengths)
    last = findlast(x->x<=line_center*(1+Δ),wavelengths)
    if isnothing(first) || isnothing(last)   return 0:0   end
    if last<first    return last:first    end
    return first:last
end

""" Return a range of columns indices with wavelengths between line_lo and line_hi """
function find_cols_to_fit(wavelengths::AbstractArray{T,1}, line_lo::Real, line_hi::Real; Δ::Real = Δλoλ_edge_pad_default) where T<:Real
    @assert line_lo < line_hi
    first = findfirst(x->x>=line_lo*(1-Δ),wavelengths)
    last = findlast(x->x<=line_hi*(1+Δ),wavelengths)
    if isnothing(first) || isnothing(last)   return 0:0   end
    if last<first    return last:first    end
    return first:last
end

""" Return list of (pixels, order_idx) pairs that contain pixels with desireed wavelengths.
    Excludes locations that contain any pixels with var == NaN.
"""
function findall_line end

function findall_line(goal::Real,lambda::AbstractArray{T1,2},var::AbstractArray{T2,2}; Δ::Real = Δλoλ_fit_line_default) where {T1<:Real, T2<:Real}
    @assert lambda[1,1] <= goal <= lambda[end,end]
    @assert size(lambda) == size(var)
    @assert Δ >= zero(Δ)
    orders = find_orders_with_line(goal,lambda)
    @assert length(orders) >= 1
    locs = map(o->(pixels=find_cols_to_fit(lambda[:,o],goal,Δ=Δ),order=o), orders)
    locs_good_idx = findall(t->!any(isnan.(var[t[1],t[2]])),locs)
    #locs_good_idx = findall(t-> !(first(t.pixels)==0 || last(t.pixels)==0 || t.order==0) && (!any(isnan.(var[t.pixels,t.order]))) ,locs)
    if length(locs) != length(locs_good_idx)
        locs = locs[locs_good_idx]
    end
    return locs
end

function findall_line(goal_lo::Real,goal_hi::Real, lambda::AbstractArray{T1,2},var::AbstractArray{T2,2}; Δ::Real = Δλoλ_edge_pad_default, verbose::Bool = false) where {T1<:Real, T2<:Real}
    @assert lambda[1,1] <= goal_lo < goal_hi <= lambda[end,end]
    orders = find_orders_with_line(goal_lo,goal_hi,lambda)
    #if ! (length(orders) >= 1) return end
    #=
    if verbose
        for i in 1:5
        println("# i= ",i," min(order)= ",minimum(lambda[:,i])," max(order)= ",maximum(lambda[:,i]), "   goal_lo= ",goal_lo, " goal_hi = ",goal_hi)
        end
    end
    flush(stdout)
    =#
    @assert length(orders) >= 1
    locs = map(o->(pixels=find_cols_to_fit(lambda[:,o],goal_lo, goal_hi,Δ=Δ),order=o), orders)
    #locs_good_idx = findall(t->!any(isnan.(var[t[1],t[2]])),locs)
    locs_good_idx = findall(t-> !(first(t.pixels)==0 || last(t.pixels)==0 || t.order==0) && (!any(isnan.(var[t.pixels,t.order]))) ,locs)
    if length(locs) != length(locs_good_idx)
        locs = locs[locs_good_idx]
    end
    return locs
end

function findall_line(goal::Real,lambda::AbstractArray{T1,1},var::AbstractArray{T2,1}; Δ::Real = Δλoλ_fit_line_default) where {T1<:Real, T2<:Real}
    @assert lambda[1] <= goal <= lambda[end]
    @assert size(lambda) == size(var)
    @assert Δ >= zero(Δ)
    locs = find_cols_to_fit(lambda,goal,Δ=Δ)
    locs_good_idx = findall(t->!any(isnan.(var[t])),locs)
    #locs_good_idx = findall(t-> !(first(t.pixels)==0 || last(t.pixels)==0 || t.order==0) && (!any(isnan.(var[t.pixels,t.order]))) ,locs)
    if length(locs) != length(locs_good_idx)
        locs = locs[locs_good_idx]
    end
    return locs
end

function findall_line(goal_lo::Real,goal_hi::Real, lambda::AbstractArray{T1,1},var::AbstractArray{T2,1}; Δ::Real = Δλoλ_edge_pad_default, verbose::Bool = false) where {T1<:Real, T2<:Real}
    @assert lambda[1] <= goal_lo < goal_hi <= lambda[end]
#=    if verbose
        for i in 1:5
        println("# i= ",i," min(order)= ",minimum(lambda[:,i])," max(order)= ",maximum(lambda[:,i]), "   goal_lo= ",goal_lo, " goal_hi = ",goal_hi)
        end
    end
    flush(stdout)
    =#
    locs = find_cols_to_fit(lambda,goal_lo, goal_hi,Δ=Δ)
    locs_good_idx = findall(t-> !any(isnan.(var[t])) ,locs)
    #locs_good_idx = findall(t-> !(first(t.pixels)==0 || last(t.pixels)==0 || t.order==0) && (!any(isnan.(var[t.pixels,t.order]))) ,locs)
    if length(locs) != length(locs_good_idx)
        locs = locs[locs_good_idx]
    end
    return locs
end

function findall_line(goal::Real,spectra::AS; Δ::Real = Δλoλ_fit_line_default) where {AS<:AbstractSpectra}
    findall_line(goal,spectra.λ,spectra.var, Δ=Δ)
end

function findall_line(goal_lo::Real,goal_hi::Real,spectra::AS; Δ::Real = Δλoλ_edge_pad_default) where {AS<:AbstractSpectra}
    findall_line(goal_lo,goal_hi,spectra.λ,spectra.var, Δ=Δ)
end

""" Return (pixels, order_idx) pair that contain "best" region of spectra, based on highest SNR. """
function find_line_best end

function find_line_best(goal::Real,lambda::AbstractArray{T1,2},flux::AbstractArray{T2,2},var::AbstractArray{T3,2}; Δ::Real = Δλoλ_fit_line_default) where {T1<:Real, T2<:Real, T3<:Real}
    locs = findall_line(goal,lambda,var,Δ=Δ)
    if length(locs) == 0   return  missing end
    #scores = map( t->sum( flux[t[1],t[2]] ./ var[t[1],t[2]])/sum( 1.0 ./ var[t[1],t[2]]), locs)
    scores = map( t->calc_snr(flux[t[1],t[2]],var[t[1],t[2]]), locs)
    idx_best = findmax(scores)
    locs[idx_best[2]]
end

function find_line_best(goal_lo::Real,goal_hi::Real, lambda::AbstractArray{T1,2},flux::AbstractArray{T2,2},var::AbstractArray{T3,2}; Δ::Real = Δλoλ_edge_pad_default) where {T1<:Real, T2<:Real, T3<:Real}
    locs = findall_line(goal_lo,goal_hi,lambda,var,Δ=Δ)
    if length(locs) == 0
        println("=>(",goal_lo, ", ",goal_hi, ")  Δ=",Δ)
        return  missing
    end
    #scores = map( t->sum( flux[t[1],t[2]] ./ var[t[1],t[2]])/sum( 1.0 ./ var[t[1],t[2]]), locs)
    scores = map( t->calc_snr(flux[t[1],t[2]],var[t[1],t[2]]), locs)
    idx_best = findmax(scores)
    locs[idx_best[2]]
end

function find_line_best(goal::Real,lambda::AbstractArray{T1,1},flux::AbstractArray{T2,1},var::AbstractArray{T3,1}; Δ::Real = Δλoλ_fit_line_default) where {T1<:Real, T2<:Real, T3<:Real}
    cols = find_cols_to_fit(lambda,goal, Δ=Δ)
    @assert( ( first(cols)==0 && last(cols)==0)  || !any(isnan.(var[cols])) )
    return cols
    #=
    locs = findall_line(goal,lambda,var,Δ=Δ)
    if length(locs) == 0   return  missing end
    #scores = map( t->sum( flux[t[1],t[2]] ./ var[t[1],t[2]])/sum( 1.0 ./ var[t[1],t[2]]), locs)
    return locs
    scores = map( t->calc_snr(flux[t[1],t[2]],var[t[1],t[2]]), locs)
    idx_best = findmax(scores)
    locs[idx_best[2]]
    =#
end

function find_line_best(goal_lo::Real,goal_hi::Real, lambda::AbstractArray{T1,1},flux::AbstractArray{T2,1},var::AbstractArray{T3,1}; Δ::Real = Δλoλ_edge_pad_default) where {T1<:Real, T2<:Real, T3<:Real}
    cols = find_cols_to_fit(lambda,goal_lo, goal_hi, Δ=Δ)
    @assert( ( first(cols)==0 && last(cols)==0)  || !any(isnan.(var[cols])) )
    return cols
    #=
    locs = findall_line(goal_lo,goal_hi,lambda,var,Δ=Δ)
    if length(locs) == 0
        println("=>(",goal_lo, ", ",goal_hi, ")  Δ=",Δ)
        return  missing
    end
    return locs
    #scores = map( t->sum( flux[t[1],t[2]] ./ var[t[1],t[2]])/sum( 1.0 ./ var[t[1],t[2]]), locs)
    scores = map( t->calc_snr(flux[t],var[t]), locs)
    idx_best = findmax(scores)
    locs[idx_best[2]]
    =#
end

function find_line_best(goal::Real,spectra::AS; Δ::Real = Δλoλ_fit_line_default) where {AS<:AbstractSpectra}
    find_line_best(goal,spectra.λ,spectra.flux,spectra.var, Δ=Δ)
end

function find_line_best(goal_lo::Real,goal_hi::Real,spectra::AS; Δ::Real = Δλoλ_edge_pad_default) where {AS<:AbstractSpectra}
    find_line_best(goal_lo,goal_hi,spectra.λ,spectra.flux,spectra.var, Δ=Δ)
end


""" Find pixels included in a range of wavelengths """
function find_pixels_for_line_in_chunk( chunk::AbstractChunkOfSpectrum, λ_min::Real, λ_max::Real )# ; plan::LineFinderPlan = LineFinderPlan() )
  idx_lo = searchsortedfirst(chunk.λ, λ_min, by=x->x>=λ_min)
  idx_tmp = searchsortedlast(chunk.λ[idx_lo:end], λ_max, by=x->x<=λ_max, rev=true)
  idx_hi = idx_lo + idx_tmp - 1
  return idx_lo:idx_hi
end

function find_pixels_for_line_in_chunklist( chunk_list::AbstractChunkList, λ_min::Real, λ_max::Real; verbose::Bool = true )
  ch_idx_all = findall(c-> (λ_min <= minimum(chunk_list.data[c].λ)) && (maximum(chunk_list.data[c].λ) <= λ_max) ,1:length(chunk_list))
  println("Hello")
  #map(c->(chunk_idx=c, pixels=find_pixels_for_line_in_chunk(chunk_list.data[c], λ_min, λ_max) ), ch_idx)
  ch_idx = 0
  if length(ch_idx_all) > 1
    snr_of_chunks_with_line = map(c->RvSpectMLBase.calc_snr(chunk_list.data[c].flux, chunk_list.data[c].var), ch_idx_all)
    ch_idx_to_keep = argmax(snr_of_chunks_with_line)
    ch_idx = ch_idx_all[ch_idx_to_keep]
    if verbose
      println(" Found λ=",λ_min,"-",λ_max," in chunks: ", ch_idx_all, " containing ", length.(ch_idx_all), " pixels.")
      println(" SNRs = ", snr_of_chunks_with_line)
      println(" Keeping chunk #",ch_idx)
    end
  elseif length(ch_idx_all) == 1
    ch_idx = first(ch_idx_all)
    if verbose
      println(" Found λ=",λ_min,"-",λ_max," in chunk: ", ch_idx, " containing ", length(ch_idx), " pixels.")
      snr_of_chunk_with_line = RvSpectMLBase.calc_snr(chunk_list.data[ch_idx].flux, chunk_list.data[ch_idx].var)
      println(" SNRs = ", snr_of_chunk_with_line)
    end

  end
  if ch_idx == 0
    error("Didn't find λ = " *string(λ_min)*" - " *string(λ_max)* " in chunklist.")

  end
  return (chunk_idx=ch_idx, pixels=find_pixels_for_line_in_chunk(chunk_list.data[ch_idx], λ_min, λ_max) )
end

function find_pixels_for_line_in_chunklist( chunk_list::AbstractChunkList, λ_min::Real, λ_max::Real, chunk_id::Integer)
  return (chunk_idx=chunk_id, pixels=find_pixels_for_line_in_chunk(chunk_list.data[chunk_id], λ_min, λ_max) )
end





""" `is_in_wavelength_range_list(λ; list )`
Return true if λ is between lambda_lo and lambda_hi for any row in list
"""
function is_in_wavelength_range_list(λ::Real; list::DataFrame  )
    @assert hasproperty(list, :lambda_lo)
    @assert hasproperty(list, :lambda_hi)
    idx =  searchsortedfirst(list[:,:lambda_hi], λ)
    return idx>size(list,1) || !(list[idx,:lambda_lo]<=λ<=list[idx,:lambda_hi]) ?  false : true
end


""" `is_in_wavelength_range_list(λ_lo, λ_hi; list )`
Return true if there is overlap between (λ_lo, λ_hi) and lambda_lo and lambda_hi for any row in list
# TODO: test
"""
function is_in_wavelength_range_list(λ_lo::Real, λ_hi::Real; list::DataFrame  )
    @assert λ_lo < λ_hi
    @assert hasproperty(list, :lambda_lo)
    @assert hasproperty(list, :lambda_hi)
    idx =  searchsortedfirst(list[:,:lambda_hi], λ_lo)
    if idx>size(list,1)    return false  end
    if λ_lo<=list[idx,:lambda_hi] &&  λ_hi>=list[idx,:lambda_lo]
        return true
    else
        return false
    end
end
