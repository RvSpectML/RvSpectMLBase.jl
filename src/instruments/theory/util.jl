import ..RvSpectMLBase.InstrumentsCommon: get_λ_range

function get_λ_range(data::ST) where { T1<:Real, T2<:Real, T3<:Real, A1<:AbstractArray{T1,2}, A2<:AbstractArray{T2,2}, A3<:AbstractArray{T3,2},
                                       IT<:AnyTheoreticalInstrument, ST<:Spectra2DBasic{T1,T2,T3,A1,A2,A3,IT} }
   (λmin, λmax) = extrema(data.λ)
   return (min=λmin, max=λmax)
end

""" `calc_λs`
Generate an array of n wavelengths uniformly spaced in log(lambda)
"""
function calc_λs end

""" `calc_λs(λ_min, λ_max, n)`
Generate an array of n wavelengths uniformly spaced in log(lambda) from λ_min to λ_max.
"""
function calc_λs(λ_min::Real, λ_max::Real, n::Integer)
   r = range(log(λ_min), stop=log(λ_max), length=n)
   return exp.(r)
end

""" `calc_λs(inst)`
Generate an array of wavelengths based on a theoretical instrument's properties
"""
function calc_λs(inst::TheoreticalInstrument1D)
   Δlnλ = 1.0/inst.resolution
   logλ_min = log(λ_min(inst))
   logλ_max = log(λ_max(inst))
   n = ceil(Int,(logλ_max-logλ_min) / Δlnλ )
   @assert 1 <= n <= max_pixel(inst)
   r = range(logλ_min, stop=logλ_max, length=n)
   return exp.(r)
end

""" `calc_λs(inst)`
Compute the wavelength for one pixel of a theoretical 1d instrument
"""
function calc_λ(pixel, inst::TheoreticalInstrument1D)
   @assert min_pixel(inst) 1<= pixel <= max_pixel(inst)
   Δlnλ = 1.0/inst.resolution
   logλ_min = log(λ_min(inst))
   logλ_max = log(λ_max(inst))
   n = ceil(Int,(logλ_max-logλ_min) / Δlnλ )
   @assert 1 <= n <= max_pixel(inst)
   r = range(logλ_min, stop=logλ_max, length=n)
   @assert 1<= pixel <= length(r)
   exp(r[pixel])
end


""" `calc_λs(inst)`
Generate an array of wavelengths based on a theoretical instrument's properties
"""
function calc_λs(inst::TheoreticalInstrument2D)
   num_orders = length(inst.λ_min)
   num_pixels = inst.pixels_per_order
   λ = zeros(num_pixels,num_orders)
   for i in 1:num_orders
      #λ_max_order = exp(logλ_min + (logλ_max-logλ_min)*(i  )/num_orders)
      #λ_min_order = exp(logλ_min + (logλ_max-logλ_min)*(i-1)/num_orders)
      λ[:,i] .= calc_λs(inst.λ_min[i], inst.λ_max[i], num_pixels)
   end
   return λ
end
