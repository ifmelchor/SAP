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


struct Base{J<:Integer}
  nite    :: J              #  --> slowness grid nite x nite
  nwin    :: J              #  --> number of time windows
  nsta    :: J              #  --> number of stations
  lwin    :: J              #  --> time window length
  citer   :: Vector{Tuple{J,J}}
  slow2   :: Bool
end


mutable struct ATF{T<:Real}
    freqs :: AbstractArray{T}
    x     :: AbstractArray{T}
    y     :: AbstractArray{T}
    sx    :: AbstractArray{T} # Horizontal slownesses in s/km in x-direction of beamforming grid (first dimension)
    sy    :: AbstractArray{T} # Horizontal slownesses in s/km in y-direction of beamforming grid (second dimension)
    power :: AbstractArray{T,2} # Power of array response at each (sx, sy) point
end

