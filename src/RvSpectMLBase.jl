"""
Author: Eric Ford and collaborators (see header of each file)
Created: August 2020
Contact: https://github.com/RvSpectML/RvSpectMLBase.jl
"""

"""
[**RvSpectMLBase.jl**](https://github.com/RvSpectML/RvSpectMLBase.jl) is the base package for a set of complementary Julia package's for measuring radial velocities (RVs) from stellar spectroscopic timeseries via machine learning (ML).

"""

module RvSpectMLBase


# Packages we'll use in many places, but try to keep these minimal.
using LinearAlgebra, Statistics
using DataFrames, Query
using NaNMath

# types.jl is responsible for exporting its own types and functions
include("types/types.jl")

# Load basic physics early, since may be used by instruments
include("util/physics.jl")
export calc_doppler_factor, apply_doppler_boost!
export absorption_line
export speed_of_light_mps

# Load basic CS algorithms early?
#include("util/cs_algs.jl")
#export searchsortednearest


# instruments.jl & the instruments it contains are responsible for exporting their own functions & modules
include("instruments/instruments.jl")

# Load the rest of the utility functions
include("util/util.jl")
#export searchsortednearest

end
