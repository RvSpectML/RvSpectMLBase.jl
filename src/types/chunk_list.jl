"""
Code declaring ChunkList, ChunkListTimeseries, and their abstract versions.

Author: Eric Ford
Created: August 2020
"""

abstract type AbstractChunkList end

""" Struct containing an array of ChunkOfSpectrum """
struct ChunkList{CT<:AbstractChunkOfSpectrum, ACT<:AbstractArray{CT,1}, OT<:Integer, AOT<:AbstractArray{OT,1} } <: AbstractChunkList
      data::ACT
      order::AOT

      function ChunkList{CT,ACT,OT,AOT}(chunks::ACT, orders::AOT) where {CT<:AbstractChunkOfSpectrum, ACT<:AbstractArray{CT,1}, OT<:Integer, AOT<:AbstractArray{OT,1} }
          #@assert length(chunks)>=1
          @assert length(chunks) == length(orders)
          new(chunks,orders)
      end
end

function ChunkList(chunks::ACT, orders::AOT) where {CT<:AbstractChunkOfSpectrum, ACT<:AbstractArray{CT,1}, OT<:Integer, AOT<:AbstractArray{OT,1} }
        ChunkList{CT,ACT,OT,AOT}(chunks,orders)
end
import Base.getindex, Base.setindex!, Base.append!

""" Provide direct access to data, an AbstractArray of ChunkOfSpectrum's via [] operator """
function getindex(x::CLT,idx) where {CLT<:AbstractChunkList}
    x.data[idx]
end

""" Provide direct access to data, an AbstractArray of ChunkOfSpectrum's via [] operator """
function setindex!(x::CLT, idx, y::CLT) where {CLT<:AbstractChunkList}
    x.data[idx] .= y.data[idx]
end

""" Append second chunk list to the first.  """
function append!(x::CLT, y::CLT) where {CLT<:AbstractChunkList}
    append!(x.data,y.data)
end

import Base.length
""" Return number of chunks in ChunkList """
length(cl::CLT) where {CLT<:AbstractChunkList} = length(cl.data)



abstract type AbstractChunkListTimeseries end

""" Mtching lists of times and array of ChunkLists """
struct ChunkListTimeseries{ TT<:Real, AT<:AbstractArray{TT,1}, CLT<:AbstractChunkList, ACLT<:AbstractArray{CLT,1}, InstT<:AbstractInstrument, AVMT<:AbstractVector{MetadataT} } <: AbstractChunkListTimeseries
    times::AT
    chunk_list::ACLT
    inst::InstT
    metadata::AVMT
end

function ChunkListTimeseries(times::AT, clt::ACLT; inst::InstT = TheoreticalInstrument1D(),
            metadata::AbstractVector{MetadataT} = fill(MetadataT,length(times)) ) where {
                TT<:Real, AT<:AbstractArray{TT,1}, CLT<:AbstractChunkList, ACLT<:AbstractArray{CLT,1},  InstT<:AbstractInstrument, AVMT<:AbstractVector{MetadataT} }
    #ChunkListTimeseries{eltype(times),typeof(times),eltype(clt),typeof(clt),typeof(inst),typeof(metadata)}(times,clt,inst,metadata)
    ChunkListTimeseries(times,clt,inst,metadata)
end

#=
function ChunkListTimeseries(t::AT, cl::ACLT) where {CLT<:AbstractChunkList, ACLT<:AbstractArray{CLT,1}, TT<:Real, AT<:AbstractArray{TT,1} }
    ChunkListTimeseries{CLT,ACLT,TT,AT}(t,cl)
end
=#
import Base.getindex, Base.setindex!

""" Allow direct access to data, an AbstractArray of ChunkOfSpectrum's via [] operator """
function getindex(x::CLT,idx) where {CLT<:AbstractChunkListTimeseries}
    x.chunk_list[idx]
end

""" Allow direct access to data, an AbstractArray of ChunkOfSpectrum's via [] operator """
function setindex!(x::CLT, idx, y::CLT) where {CLT<:AbstractChunkListTimeseries}
    x.times[idx] .= y.times[idx]
    x.chunk_list[idx] .= y.chunk_list[idx]
    x.metadata[idx] .= y.metadata[idx]
end

""" Return number of times/ChunkLists in a ChunkListTimeseries """
length(clts::ACLT) where {ACLT<:AbstractChunkListTimeseries} = length(clts.chunk_list)
""" Return number of times/ChunkLists in a ChunkListTimeseries """
num_times(clts::ACLT) where {ACLT<:AbstractChunkListTimeseries} = length(clts.chunk_list)

""" Number of chunks in first chunk_list in chunk_list_timeseries. """
num_chunks(chunk_list_timeseries::ACLT) where {ACLT<:AbstractChunkListTimeseries} = length(first(chunk_list_timeseries.chunk_list))

""" Get time of ith observation in chunk_list_timeseries.  """
time(chunk_list_timeseries::ACLT, i::Integer) where {ACLT<:AbstractChunkListTimeseries} = chunk_list_timeseries.times[i]

""" Get metadata[:key] of ith observation in chunk_list_timeseries.  """
metadata(chunk_list_timeseries::ACLT, i::Integer, key::Symbol) where {ACLT<:AbstractChunkListTimeseries} = chunk_list_timeseries.metadata[i][key]

""" Get instrument used for chunk_list_timeseries.  """
instrument(chunk_list_timeseries::ACLT) where {ACLT<:AbstractChunkListTimeseries} = chunk_list_timeseries.inst

""" Make a ChunkListTimeseries containing a subset of observations from an existing  ChunkListTimeseries """
function extract_chunklist_timeseries_with_subset_obs( clt::AbstractChunkListTimeseries, ref_obs_idx::Union{Integer,AbstractVector{Integer}} ) where { CLT<:AbstractChunkListTimeseries }
  if typeof(ref_obs_idx) <: Integer
      @assert 1 <= ref_obs_idx <= num_times(clt)
  else
      @assert all(1 .<= ref_obs_idx .<= num_times(clt))
  end
  ChunkListTimeseries([clt.times[ref_obs_idx]], [clt.chunk_list[ref_obs_idx]], clt.inst, [clt.metadata[ref_obs_idx]]   )
end
