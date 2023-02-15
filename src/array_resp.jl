#!/usr/local/bin julia
# coding=utf-8

# Utility functions for .jl

# These functions were taken from Beamforming julia package
# see more in github.com/anowacki/Beamforming.jl

# Ivan Melchor 2023


function array_response(x::AbstractArray{T}, y::AbstractArray{T}, smax::T, ds::T, f1::T, f2::T, df::T) where T<:Real

   freqs = f1:df:f2
   s     = -smax:ds:smax
   power = zeros(T, length(s), length(s))
   arf = ArrayResponse(freqs, x, y, s, s, power)
    
   return _compute_arf!(arf)
end


"""
    _compute_arf!(arf)
Do the actual ARF calculation for an `ArrayResponse`.
"""
function _compute_arf!(arf::ArrayResponse{T}) where T
    @inbounds for (j, slowy) in enumerate(arf.sy), (i, slowx) in enumerate(arf.sx)
        for f in arf.freqs
            ω = T(2)*π*f
            pow = zero(T)*im
            @simd for istat in eachindex(arf.x, arf.y)
                pow += cis(ω*(slowx*arf.x[istat] + slowy*arf.y[istat]))
            end
            arf.power[i,j] += abs2(pow)
        end
    end
    arf
end


