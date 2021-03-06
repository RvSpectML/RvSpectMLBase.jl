
#import ..RvSpectML:   # Already in module
import ..RvSpectMLBase: predict_intrinsic_stellar_line_width, speed_of_light_mps

"""  `generate_spectrum(line_list, inst; time, rv, ssbz, snr_per_pixel, line_width, add_noise )`
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
   metadata[:normalization] = :continuum
   if typeof(inst) <: AbstractInstrument1D

      return Spectra1DBasic(λ,flux,var,inst, metadata)
   elseif typeof(inst) <: AbstractInstrument2D
      return Spectra2DBasic(λ,flux,var,inst, metadata)
   else
      @error("Invalid instrument type: ", typeof(inst))
   end
end

"""  `generate_spectra_timeseries(line_list, inst; time, rv, ssbz, snr_per_pixel, line_width, add_noise )`
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


""" `demo_generate_spectrum_line( inst )`
Generate spectrum with one line for testing purposes.
"""
function demo_generate_spectrum_line(inst::AnyTheoreticalInstrument)
   ll = DataFrame(:lambda=>[5500.0], :weight=>[0.6] )
   generate_spectrum(ll, inst)
end
