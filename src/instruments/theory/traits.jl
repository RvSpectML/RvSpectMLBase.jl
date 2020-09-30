"""
   Traits for a theoretical spectrograph
Author: Eric Ford
Created: August 2020
"""

""" Delegates loading of code specifying types essential to the package.  """

import ..RvSpectMLBase.InstrumentsCommon: min_pixel, max_pixel, min_order, max_order, min_pixel_in_order, max_pixel_in_order
import ..RvSpectMLBase.InstrumentsCommon: orders_to_use_default, min_col_default, max_col_default

min_pixel(::TheoreticalInstrument1D) = 1
max_pixel(::TheoreticalInstrument1D) = 128*8192
#min_order(::TheoreticalInstrument1D) = 1
#max_order(::TheoreticalInstrument1D) = 1
min_pixel_in_order(::TheoreticalInstrument1D) = 1
max_pixel_in_order(::TheoreticalInstrument1D) = 8192
#orders_to_use_default(inst::TheoreticalInstrument1D) = 1
#min_col_default(::TheoreticalInstrument1D ) = 1
#max_col_default(inst::TheoreticalInstrument1D ) = max_pixel_in_order(inst)
#min_col_default(::TheoreticalInstrument1D, ord::Integer) = min_col_default(inst)
#max_col_default(inst::TheoreticalInstrument1D, ord::Integer) = max_col_default(inst)

min_order(::TheoreticalInstrument2D) = 1
max_order(inst::TheoreticalInstrument2D) = length(inst.λ_min)
min_pixel(::TheoreticalInstrument2D) = 1
max_pixel(inst::TheoreticalInstrument2D) = inst.pixels_per_order*length(inst.λ_min)
min_pixel_in_order(::TheoreticalInstrument2D) = 1
max_pixel_in_order(inst::TheoreticalInstrument2D) = inst.pixels_per_order

orders_to_use_default(inst::TheoreticalInstrument2D) = min_order(inst):max_order(inst)
min_col_default(::TheoreticalInstrument2D, ord::Integer) = 1
max_col_default(inst::TheoreticalInstrument2D, ord::Integer) = max_pixel_in_order(inst)

#import ..RvSpectMLBase: metadata_symbols_default, metadata_strings_default
#metadata_symbols_default(::AnyD) = Symbol[:bjd, :target, :ssbz]
#metadata_strings_default(::AnyD) = String["OBSJD", "SKY-OBJ", "SSBZ000"]

import ..RvSpectMLBase.InstrumentsCommon: default_ccf_mask_v_width
default_ccf_mask_v_width(::AnyTheoreticalInstrument) = 500.0  #
export default_ccf_mask_v_width

import ..RvSpectMLBase.InstrumentsCommon: get_inst_module
get_inst_module(::AnyTheoreticalInstrument) = TheoreticalInstrument
export get_inst_module

#import ..RvSpectMLBase: filter_line_list, find_worst_telluric_in_each_chunk
# No tellurics in theoretical spectra
#export filter_line_list, find_worst_telluric_in_each_chunk
