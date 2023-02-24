#!/usr/local/bin julia
# coding=utf-8

# GNU GPL v2 licenced to I. Melchor and J. Almendros 08/2022


function RelSE(data::Array{T}, xsta::Vector{T}, ysta::Vector{T}, sta_ref::J, data_master::Array{T}, s_master::Vector{T}, pmax::Vector{T}, pinc::Vector{T}, fsem::J, lwin::J, nwin::J, nadv::T, toff::J) where {T<:Real, J<:Integer}

    nsta   = length(xStaUTM)
    slow_m, baz_m = s_master
    pxy_m = [slow_m*sin(baz_m/pi), slow_m*cos(baz_m/pi)]

    # delays of master event arrivals to array stations
    dtime_func = _dtimefunc(xsta, ysta, fsem)
    dt_m = dtime_func(pxy_m)

    # define a window ?? bow, hanning, etc... using DPS.jl

    # iterate over time
    for nk in 1:nwin
        n0  = 1 + toff + floor(Int64, lwin*nadv*(nk-1))

        # nini1 and nini2 are init time for two files
        ni1 = n0 + nini1
        ni2 = n0 + nini2

        # cross-correlation

        for ii in 1:nsta
            mi1 = ni1 + dt_m[ii]
            mi2 = ni2 + dt_m[ii]
            
            d11 = data1[ii,mi1:mi1+lwin]
            d22 = data2[ii,mi2:mi2+lwin]
            
            c11 = dot(d11,d11)
            c12 = dot(d11,d22)
            c22 = dot(d22,d22)

            cc = c12/sqrt(c11*c22)












end