
#import ..RvSpectML:   # Already in module
import ..RvSpectMLBase: predict_intrinsic_stellar_line_width, speed_of_light_mps

""" calc_λs(λ_min, λ_max, n)
Generate an array of n wavelengths uniformly spaced in log(lambda) from λ_min to λ_max.
"""
function calc_λs(λ_min::Real, λ_max::Real, n::Integer)
   r = range(log(λ_min), stop=log(λ_max), length=n)
   return exp.(r)
end

""" calc_λs(inst)
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

""" calc_λs(inst)
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


""" calc_λs(inst)
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

"""  generate_spectrum(line_list, inst; time, rv, ssbz, snr_per_pixel, line_width, add_noise )
Generate a spectrum using a line_list (as DataFrame with columns lambda and weight) and theoretical instrument.
Optionally, specify times, rv's, barycentric corrections, etc.
"""
function generate_spectrum(line_list::DataFrame, inst::AnyTheoreticalInstrument; time::Real = 0.0, rv::Real = 0.0, ssbz::Real = 0.0,
               snr_per_pixel::Real = 1000.0, line_width::Real = predict_intrinsic_stellar_line_width(5780), add_noise::Bool = true ) # where { LLT<:AbstractLineList }
   @assert hasproperty(line_list,:lambda)
   @assert hasproperty(line_list,:weight)
   λ = calc_λs(inst)
   flux = ones(size(λ))
   #var = ones(size(λ))
   doppler_factor = calc_doppler_factor(rv) * (1+ssbz)  # TODO: Fix problem... Is ssbz in m/s or v/c?
   for i in 1:size(line_list,1)
      #width = doppler_factor* line_list.lambda[i]*line_width*1000/speed_of_light_mps   #TODO PUT PACK just for testing
      width = doppler_factor*line_list.lambda[i]*line_width*1000/speed_of_light_mps
      flux .*= RvSpectMLBase.absorption_line.(λ, mid=line_list.lambda[i]*doppler_factor, depth=line_list.weight[i], width=width)
   end
   if add_noise
      flux .+= randn(size(flux)).*sqrt.(flux)./snr_per_pixel
      flux[flux .< 0.0] .= 0.0
   end
   var = flux ./ snr_per_pixel^2
   metadata = Dict{Symbol,Any}( :rv_true=>rv, :snr_per_pixel=>snr_per_pixel, :line_width=>line_width, :bjd=>time, :target=>"Simulation", :ssbz=>ssbz)
   if typeof(inst) <: AbstractInstrument1D

      return Spectra1DBasic(λ,flux,var,inst, metadata)
   elseif typeof(inst) <: AbstractInstrument2D
      return Spectra2DBasic(λ,flux,var,inst, metadata)
   else
      @error("Invalid instrument type: ", typeof(inst))
   end
end

"""  generate_spectra_timeseries(line_list, inst; time, rv, ssbz, snr_per_pixel, line_width, add_noise )
Generate a time series of spectra using times, a line_list (as DataFrame with columns lambda and weight) and a theoretical instrument.
Optionally, specify times, rv's, barycentric corrections, etc.
"""
function generate_spectra_timeseries(times::AbstractArray, line_list::DataFrame, inst::AnyTheoreticalInstrument, rvs::AbstractArray; ssbzs::AbstractArray = zeros(length(rvs)),
               snr_per_pixel::Real = 1000.0, line_width::Real = predict_intrinsic_stellar_line_width(5780) ) # where { LLT<:AbstractLineList }
  @assert length(times) == length(rvs) == length(ssbzs)
  @assert hasproperty(line_list,:lambda)
  @assert hasproperty(line_list,:weight)
  spectra = map(i->generate_spectrum(line_list, inst, snr_per_pixel=snr_per_pixel, line_width=line_width, time=times[i], rv=rvs[i], ssbz=ssbzs[i] ), 1:length(times) )
end


""" demo_generate_spectrum_line
Generate spectrum with one line for testing purposes.
"""
function demo_generate_spectrum_line(inst::AnyTheoreticalInstrument)
   ll = DataFrame(:lambda=>[5500.0], :weight=>[0.6] )
   generate_spectrum(ll, inst)
end
