#!/usr/local/bin julia
# coding=utf-8

# GNU GPL v2 licenced to I. Melchor and J. Almendros 08/2022

__precompile__()

module SAP

    using LinearAlgebra
    using Statistics
    
    export zlcc, get_dtimes, array_transfunc, p2r

    include("types.jl")

    include("utils.jl")

    include("array_resp.jl")

    include("cc8mre.jl")

end
