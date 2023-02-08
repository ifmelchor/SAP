#!/usr/local/bin julia
# coding=utf-8

# GNU GPL v2 licenced to I. Melchor and J. Almendros 08/2022

# AbstractArrayResponse and ArrayResponse copied from github.com/anowacki/Beamforming.jl

mutable struct Bounds{T<:Real}
  azimin :: T
  azimax :: T
  slomin :: T
  slomax :: T
end


struct xySta{T<:Real}
  x :: T
  y :: T
end


struct BaseParams{T<:Real,J<:Int}
  data    :: Array{T}
  stalist :: Vector{xySta}
  xref    :: T
  yref    :: T
  ccerr   :: T                      #  --> correlation level
  fsem    :: J                      #  --> sampling rate
  lwin    :: J                      #  --> time window length
  nwin    :: J                      #  --> number of time windows
  nadv    :: T                      #  --> percentage of time advance (0--1)
  pmax    :: Vector{T}              #  --> maximum slownes for the grid
  pinc    :: Vector{T}              #  --> slownes interval for the grid
  citer   :: Vector{Tuple{J,J}}
  fqbands :: Vector{Tuple{T,T}}
end


"""
    AbstractArrayResponse
Supertype of all array response types.  All subtypes should contain the following
fields:
- `freqs`: Frequencies at which the array response was evaluated (Hz)
- `sx`: Horizontal slowness in s/° in x-direction of grid (first dimension)
- `sy`: Horizontal slowness in s/° in y-direction of grid (second dimension)
- `power`: Power of array response at each `(sx, sy)` point.
"""
abstract type AbstractArrayResponse end

"""
    ArrayResponse{T}
Struct containing the array response function for a set of stations.
Slownesses in an `ArrayResponse` are always in s/°, and coordinates are
in km.
This can be plotted using [`Plots.plot`](@ref).
"""
struct ArrayResponse{T} <: AbstractArrayResponse
    "Frequencies at which the array response was evaluated (Hz)"
    freqs::Vector{T}
    "Cartesian x-coordinates in each station (km)"
    x::Vector{T}
    "Cartesian y-coordinates in each station (km)"
    y::Vector{T}
    "Horizontal slownesses in s/km in x-direction of beamforming grid (first dimension)"
    sx::Vector{T}
    "Horizontal slownesses in s/km in y-direction of beamforming grid (second dimension)"
    sy::Vector{T}
    "Power of array response at each (sx, sy) point"
    power::Array{T,2}
end

