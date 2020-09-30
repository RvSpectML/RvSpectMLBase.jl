```@meta
CurrentModule = RvSpectMLBase.InstrumentsCommon
```
# Functions declared in RvSpectMLBase and specialized for each instrument

```@contents
Pages = ["abstract_functions.md"]
Depth = 2
```

!!! todo
    Figure out what Documenter can't find the docstrings for these.
    They're in src/instruments/common.jl.

## Instrument-specific traits
```@docs
min_order
max_order
orders_all
orders_to_use_default
min_pixel_in_order
max_pixel_in_order
min_pixel
max_pixel
pixels_all
min_pixels_in_chunk
max_pixels_in_spectra
min_col_default
max_col_default
metadata_symbols_default
metadata_strings_default
default_ccf_mask_v_width
get_inst_module
```

## Instrument-specific computation
```@docs
get_λ_range
filter_line_list
```
