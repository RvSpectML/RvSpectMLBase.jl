"""
Authors: Various
Compiled by: Eric Ford (see each function for credits)
Created: August 2020
Contact: https://github.com/eford/
"""

# util/physics.jl included separatesly from RvSpectMLBase since may be used early
#include("util/physics.jl")

include("find_line.jl")
# intended for internal use, so functions not exported

include("chunks.jl")
export make_orders_into_chunks, make_grid_for_chunk,  extract_chunk_list_timeseries_for_order
export make_chunk_list_around_lines, make_chunk_list_timeseries_around_lines, make_chunk_list_tellric_free, make_chunk_list_timeseries_tellric_free
export make_order_list_timeseries
export find_orders_with_line, find_orders_in_range
export findall_line, find_line_best
export filter_bad_chunks, find_pixels_for_line_in_chunk, find_pixels_for_line_in_chunklist
export calc_normalization, calc_normalization_var_weighted

include("spectra.jl")
export apply_doppler_boost!, normalize_spectrum!, normalize_spectra!
export get_λ_range, get_inst, make_vec_metadata_from_spectral_timeseries, discard_large_metadata
#export calc_snr

include("cs_algs.jl")
# Don't export to avoid risk of conflicts
# export searchsortednearest

include("files.jl")
export read_data_paths, code_to_include_param_jl
#export make_manifest

include("pipeline.jl")
using .Pipeline
export PipelinePlan
export make_plot,  save_plot, save_data, need_to, has_cache, read_cache, set_cache!  # Query pipeline
export need_to!, dont_need_to!, reset_all_needs!             # Write to pipeline
export make_plot!, dont_make_plot!, make_all_plots!,  make_no_plots!
