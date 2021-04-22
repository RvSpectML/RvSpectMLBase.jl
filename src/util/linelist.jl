function assign_lines_to_orders(line_list::DataFrame, order_info::DataFrame; v_center::Real=0, Δv_to_avoid_tellurics::Real = RvSpectMLBase.max_bc )
  @assert hasproperty(line_list, :lambda)
  @assert hasproperty(line_list, :lambda_lo)
  @assert hasproperty(line_list, :lambda_hi)
  ll1 = copy(line_list)
  ll2 = copy(line_list)
  ll3 = copy(line_list)
  ll1[!,:order] .= 0
  ll2[!,:order] .= 0
  ll3[!,:order] .= 0
  if !(hasproperty(line_list,:lambda_lo) && hasproperty(line_list,:lambda_hi))
    add_line_boundaries_to_line_list(line_list; Δv_to_avoid_tellurics=Δv_to_avoid_tellurics, v_center_to_avoid_tellurics=v_center )
  end
  for (i,line) in enumerate(eachrow(line_list))
    order_indices = findall( map(order-> order.λ_min <= line.lambda_lo <= order.λ_max && order.λ_min <= line.lambda_hi <= order.λ_max, eachrow(order_info) ) )
    if length(order_indices) == 1
      ll1[i,:order] = order_info.order[order_indices[1]]
    elseif length(order_indices) == 2
      ll1[i,:order] = order_info.order[order_indices[1]]
      ll2[i,:order] = order_info.order[order_indices[2]]
    elseif length(order_indices) == 3
      ll1[i,:order] = order_info.order[order_indices[1]]
      ll2[i,:order] = order_info.order[order_indices[2]]
      ll3[i,:order] = order_info.order[order_indices[3]]
    elseif length(order_indices) > 3
      @warn "One line appears in >3 orders.  Not implemented yet."
      println("line = ", line, "  order_indices = ", order_indices, "  line.order = ",order_info.order[order_indices])
      #=  df_tmp = DataFrame()
      for key in keys(eachcol(line_list))
        df_tmp[!,key] = fill(line[key], length(order_indices)-1 )
      end
      df_tmp[:order] .= order_indices[2:end]
      append!(line_list,df_tmp)
      =#
      #break
    end
  end
  append!(ll1,ll2)
  append!(ll1,ll3)
  ll1 |> @filter( _.order >= 1 ) |> @orderby(_.order) |> @thenby(_.lambda_lo) |> DataFrame
  #ll1 |> @filter( _.order >= 1 ) |> @orderby(_.lambda_lo) |> @thenby(_.order) |> DataFrame
end

function add_line_boundaries_to_line_list(line_list::DataFrame; Δv_to_avoid_tellurics::Real = 0.0, v_center_to_avoid_tellurics::Real = 0.0 )
    @assert 0.5*RvSpectMLBase.max_bc_earth_rotation <= Δv_to_avoid_tellurics <= 4*RvSpectMLBase.max_bc
    @assert hasproperty(line_list,:lambda)
    @assert !hasproperty(line_list,:lambda_lo)
    @assert !hasproperty(line_list,:lambda_hi)
    line_list_to_search_for_tellurics = copy(line_list)
    line_list_to_search_for_tellurics.lambda    = line_list_to_search_for_tellurics.lambda.*(calc_doppler_factor(v_center_to_avoid_tellurics))
    line_list_to_search_for_tellurics.lambda_lo = line_list_to_search_for_tellurics.lambda.*(calc_doppler_factor(v_center_to_avoid_tellurics)/calc_doppler_factor(Δv_to_avoid_tellurics))
    line_list_to_search_for_tellurics.lambda_hi = line_list_to_search_for_tellurics.lambda.*(calc_doppler_factor(v_center_to_avoid_tellurics)*calc_doppler_factor(Δv_to_avoid_tellurics))
    return line_list_to_search_for_tellurics
end

function expand_line_boundaries_in_line_list(line_list::DataFrame; Δv_to_avoid_tellurics::Real = 0.0, v_center_to_avoid_tellurics::Real = 0.0 )
    @assert 0.5*RvSpectMLBase.max_bc_earth_rotation <= Δv_to_avoid_tellurics <= 4*RvSpectMLBase.max_bc
    @assert hasproperty(line_list,:lambda)
    @assert hasproperty(line_list,:lambda_lo)
    @assert hasproperty(line_list,:lambda_hi)
    line_list_to_search_for_tellurics = copy(line_list)
    line_list_to_search_for_tellurics.lambda    = line_list_to_search_for_tellurics.lambda.*(calc_doppler_factor(v_center_to_avoid_tellurics))
    line_list_to_search_for_tellurics.lambda_lo = line_list_to_search_for_tellurics.lambda_lo.*(calc_doppler_factor(v_center_to_avoid_tellurics)/calc_doppler_factor(Δv_to_avoid_tellurics))
    line_list_to_search_for_tellurics.lambda_hi = line_list_to_search_for_tellurics.lambda_hi.*(calc_doppler_factor(v_center_to_avoid_tellurics)*calc_doppler_factor(Δv_to_avoid_tellurics))
    return line_list_to_search_for_tellurics
end
