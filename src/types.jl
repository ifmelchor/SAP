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
  stalist :: Vector{xySta{T}}
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
  toff    :: J
end


struct ArrayResponse{T<:Real}
    freqs :: AbstractArray{T}                #  Frequencies at which the array response was evaluated (Hz)
    x     :: AbstractArray{T}                #  Cartesian x-coordinates in each station (km)
    y     :: AbstractArray{T}                #  Cartesian y-coordinates in each station (km)
    sx    :: AbstractArray{T}                #  Horizontal slownesses in s/km in x-direction of beamforming grid (first dimension)
    sy    :: AbstractArray{T}                #  Horizontal slownesses in s/km in y-direction of beamforming grid (second dimension)
    power :: AbstractArray{T,2}               #  Power of array response at each (sx, sy) point
    
end

