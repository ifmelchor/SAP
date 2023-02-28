#!/usr/local/bin julia
# coding=utf-8

# Utility functions for .jl

# These functions were modified from Beamforming julia package
# see more in github.com/anowacki/Beamforming.jl

# Ivan Melchor 2023


function array_response(x::Array{T}, y::Array{T}, smax::T, ds::T, f1::T, f2::T, df::T) where T<:Real

   freqs = f1:df:f2
   s     = -smax:ds:smax
   power = zeros(T, length(s), length(s))
   
   @inbounds for (j, slowy) in enumerate(s), (i, slowx) in enumerate(s)
        for f in freqs
            ω = T(2)*π*f
            pow = zero(T)*im
            for istat in eachindex(x, y)
                pow += cis(ω*(slowx*x[istat] + slowy*y[istat]))
            end
            power[i,j] += abs2(pow)
        end

    end
   return ArrayResponse(freqs, x, y, s, s, power)
end


"""
    _compute_arf!(arf)
Do the actual ARF calculation for an `ArrayResponse`.
"""
# function _compute_arf!(arf::ArrayResponse)
#     @inbounds for (j, slowy) in enumerate(arf.sy), (i, slowx) in enumerate(arf.sx)
#         for f in arf.freqs
#             ω = Float64(2)*π*f
#             pow = zero(Float64)*im
#             @simd for istat in eachindex(arf.x, arf.y)
#                 pow += cis(ω*(slowx*arf.x[istat] + slowy*arf.y[istat]))
#             end
#             arf.power[i,j] += abs2(pow)
#         end
#     end
#     arf
# end


