"""
Delegates loading of code with functions and parameters specific to different instruments.
Subdirectories of src/instruments include provide functions specialized for each instrument,
typically the file I/O and pre-processing, so data ends up in a common format.
src/instruments/common.jl provides routines that can be shared by instruments.
"""

module InstrumentsCommon

using ..RvSpectMLBase
using DataFrames

include("common.jl")   # Mostly trait functions to be specialized by instruments

#include("io.jl")
#export read_manifest, read_header, read_metadata_from_fits

#include("linelists.jl")
#export read_linelist_espresso, read_linelist_vald

#include("masks.jl")
#export ChunkWidthFixedΔlnλ
#export read_mask_espresso, read_mask_vald

#export predict_intrinsic_stellar_line_width
#export λ_vac_to_air, λ_air_to_vac
#export find_overlapping_chunks
#export merge_chunks

include("param.jl")

export min_order, max_order, min_pixel_in_order, max_pixel_in_order
export orders_to_use_default, min_col_default, max_col_default, get_pixel_range
export orders_all, pixels_all, max_pixels_in_spectra       # generic implementations avaliable
export metadata_symbols_default, metadata_strings_default  # need to specialize
export default_ccf_mask_v_width
export get_inst_module

end # module InstrumentsCommon

# Now returning to module EchelleInstrments

#=
For each instrument/data file type, users need to create a sub-type of either AbstractInstrument2D or AbstractInstrument1D and
to implement the functions in neid/traits.jl for each instrument trait.
=#

include("theory/theory.jl")
import .TheoreticalInstrument: TheoreticalInstrument1D, TheoreticalInstrument2D, AnyTheoreticalInstrument
import .TheoreticalInstrument: get_inst_module
#import .TheoreticalInstrument: filter_line_list, find_worst_telluric_in_each_chunk
export TheoreticalInstrument, TheoreticalInstrument1D, TheoreticalInstrument2D, AnyTheoreticalInstrument
#=
import .TheoreticalInstrument: min_order, max_order, min_pixel_in_order, max_pixel_in_order
import .TheoreticalInstrument: orders_to_use_default, min_col_default, max_col_default
import .TheoreticalInstrument: orders_all, pixels_all, max_pixels_in_spectra       # generic implementations avaliable
import .TheoreticalInstrument: metadata_symbols_default, metadata_strings_default  # need to specialize
import .TheoreticalInstrument: default_ccf_mask_v_width
import .TheoreticalInstrument: get_inst_module
=#

#=
include("tellurics.jl")
export find_worst_telluric_in_each_chunk, make_clean_line_list_from_tellurics_expres
export filter_line_list
=#

import .InstrumentsCommon: pixels_all, min_pixel, max_pixel, max_pixels_in_spectra, min_col_default, max_col_default, get_pixel_range
import .InstrumentsCommon: orders_all, min_order, max_order, orders_to_use_default
import .InstrumentsCommon: min_pixel_in_order, max_pixel_in_order, min_pixels_in_chunk
import .InstrumentsCommon: metadata_symbols_default, metadata_strings_default
import .InstrumentsCommon: default_ccf_mask_v_width
import .InstrumentsCommon: make_clean_line_list_from_tellurics, choose_obs_idx_for_init_guess

export min_order, max_order, min_pixel_in_order, max_pixel_in_order, get_pixel_range
export orders_to_use_default, min_col_default, max_col_default
export orders_all, pixels_all, max_pixels_in_spectra       # generic implementations avaliable
export metadata_symbols_default, metadata_strings_default  # need to specialize
export default_ccf_mask_v_width
export make_clean_line_list_from_tellurics
export choose_obs_idx_for_init_guess
export get_inst_module
