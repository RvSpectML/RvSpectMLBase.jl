# Default values shared across instruments
#const teff_solar = 5780.0
#const vrot_solar = 1800.0
#const default_line_width_mps = predict_intrinsic_stellar_line_width(teff_solar,v_rot=vrot_solar)  # m/s

default_chunk_size_factor = 3       # For default_calc_chunk_width TODO: Figure out what value to use.  Ask Alex
# Is line below still used?
#default_calc_chunk_width = ChunkWidthFixedΔlnλ(default_chunk_size_factor*default_line_width_mps/speed_of_light_mps)  # Used by read_mask_espresso and read_vald for assigning lambda_lo and lambda_hi
default_min_chunk_Δv = 20000           # m/s  for ChunkWidthFixedΔlnλ

default_Δv_to_avoid_tellurics = 15000.0  # m/s

 # TODO: OPT: Make above const once settle on good values
