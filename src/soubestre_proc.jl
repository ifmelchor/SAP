

"""
   net_lte_run(*args)

Funcion LTE-network para calcular espectrograma cruzado segun Soubestre et al.
"""
# function net_lte_run(s_data::Array{Float64,2}, d_data::Union{Array{Float64,1}, Nothing}, fs::Integer, nwin::Integer, lwin::Integer, nsubwin::Union{Integer,Nothing}, lswin::Union{Integer,Nothing}, nadv::Float64, fq_band::Tuple{Float64,Float64}, NW::Float64, pad::Float64, add_param::Bool, polar::Bool, pe_order::Integer, pe_delta::Integer, ap_twin::Float64, ap_th::Float64)

#     # comp info:
#     #  1 --> Z
#     #  2 --> 
#     #  3 --> 

#     # compute the frequency domain
#     if lswin
#         freq, fqr = _fqbds(lswin, fs, fq_band, pad=pad)
#     else
#         freq, fqr = _fqbds(lwin, fs, fq_band, pad=pad)
#     end

#     # define base
#     nfs = size(freq, 1)
#     base = LTEBase(s_data, d_data, fs, freq, fq_band, fqr, nfs, nwin, lwin, nadv, nini, NW, pad, add_param, polar, pe_order, pe_delta, ap_twin, ap_th)
    
#     # run lte
#     lte = _run(base)

#     return lte
# end
