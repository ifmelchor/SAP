#!/usr/local/bin julia
# coding=utf-8

# Utility functions for .jl

# Ivan Melchor 2023

"""
    array_transfunc(arf)
    Compute the array transfer function.
    slomax/slinc in s/km
    x/y positions of the sensors in km
"""
function array_transfunc(x::Array{T}, y::Array{T}, slomax::T, sloinc::T, fmin::T, fmax::T, finc::T) where T<:Real

   freqs = fmin:finc:fmax
   s     = -slomax:sloinc:slomax
   power = zeros(T, length(s), length(s))
   nsta  = length(x) # == length(y)
   
   @inbounds for (j, sy) in enumerate(s), (i, sx) in enumerate(s)
        for f in freqs
            ω = T(2)*π*f
            pow = zero(T)*im
            for j in 1:nsta
                pow += cis(ω*(sx*x[j] + sy*y[j]))
            end
            power[i,j] += abs2(pow/nsta)
        end
        power[i,j] /= length(freqs)
    end

   return ATF(freqs, x, y, s, s, power)
end


