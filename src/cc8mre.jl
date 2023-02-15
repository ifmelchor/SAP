#!/usr/local/bin julia
# coding=utf-8

# GNU GPL v2 licenced to I. Melchor and J. Almendros 08/2022


# me gustaría añadir varias cosas
#  1: loop sobre rangos de frequencia. Esto hace que debamos filtrar directamente en julia y no en seisvo, como esta hecho. 
#  2: guardar info para poder calcular la funcion de transferencia por ventana de calculo


function CC8(data::Union{Array{T}, Nothing}, xStaUTM::Vector{T}, yStaUTM::Vector{T}, pmax::Vector{T}, pinc::Vector{T}, fsem::J, lwin::J, nwin::J, nadv::T, ccerr::T, toff::J) where T<:Real where J<:Integer
    
    length(xStaUTM) == length(yStaUTM) ||
        throw(DimensionMismatch("xStaUTM and yStaUTM coordinate vectors not the same length"))

    # define baseparams
    nsta = length(xStaUTM)
    xysta = [xySta(xStaUTM[s], yStaUTM[s]) for s in 1:nsta]
    refxy = refsta(xysta)

    # array resp
    # if aresp
    #     pmax_f = findmax(pmax)
    #     smax   = pmax_f[1]
    #     ds  = pinc[pmax_f[2]]
    #     f1  = findmin([x[1] for x in fqband])[1]
    #     f2  = findmax([x[2] for x in fqband])[1]
    #     arf = array_response(xStaUTM, yStaUTM, smax, ds, f1, f2, 0.1)
    # else
    #     arf = nothing
    # end

    cciter = Vector{Tuple{Int,Int}}()
    for ii in 1:nsta-1
        for jj in ii+1:nsta
            push!(cciter, (ii, jj))
        end
    end

    base = BaseParams(data, xysta, refxy[1], refxy[2], ccerr, fsem, lwin, nwin, nadv, pmax, pinc, cciter, toff*fsem)
    res = _core(base)

    return res
end


function _idx(ip::J, nite::J, pxy0::Vector{T}, base::BaseParams) where J<:Integer where T<:Real 
    pxx = Array{Float64}(undef, nite)
    pyy = Array{Float64}(undef, nite)
    for i in 1:nite
        # for x
        px = pxy0[1] - base.pmax[ip] + base.pinc[ip]*(i-1)
        pxx[i] = base.pinc[ip] * floor(Int64, px/base.pinc[ip])
        # for y
        py = pxy0[2] - base.pmax[ip] + base.pinc[ip]*(i-1)
        pyy[i] = base.pinc[ip] * floor(Int64, py/base.pinc[ip])
    end

    return [[px, py] for px in pxx for py in pyy]
end


"""
   _pccorr(*args)

Calcula la correlacion cruzada para un slowness y banda de frequencia
"""
function _pccorr(nkk::J, pxy::Vector{T}, base::BaseParams) where J<:Integer where T<:Real
    # compute delat times for each station
    dtimes = [pxy[1]*(sta.x-base.xref) + pxy[2]*(sta.y-base.yref) for sta in base.stalist]
    
    # filter data
    # fq_band = base.fqbands[nfq]

    # build cc matrix
    nsta = length(base.stalist)
    cc = zeros(Real, nsta, nsta)
    for ii in 1:nsta
        mii = nkk + floor(Int64, base.fsem*dtimes[ii])
        dii = @view base.data[ii, 1+mii:base.lwin+mii]
        
        for jj in ii:nsta
            mjj = nkk + floor(Int64, base.fsem*dtimes[jj])
            djj = @view base.data[jj, 1+mjj:base.lwin+mjj]
            cc[ii,jj] += dot(dii,djj)
        end
    end

    # computes crosscorr coefficient
    suma = sum([cc[ii,jj]/sqrt(cc[ii,ii]*cc[jj,jj]) for (ii, jj) in base.citer])
    return (2*suma + nsta) / nsta^2
end


function _rms(nkk::J, pxy::Vector{T}, base::BaseParams) where J<:Integer where T<:Real

    dtimes = [pxy[1]*(sta.x-base.xref) + pxy[2]*(sta.y-base.yref) for sta in base.stalist]

    erg = 0.
    nsta = length(base.stalist)
    for ii in 1:nsta
        mii = nkk + floor(Int64, base.fsem*dtimes[ii])
        dii = @view base.data[ii, 1+mii:base.lwin+mii]
        erg += sqrt(mean(dii.*dii))
    end
    
    return erg /= nsta
end


"""
   _empty_dict(*args)

Genera un dict vacio para llenar durante el procesado.
"""
function _empty_dict(base::BaseParams)
    dict = Dict()
    nip = length(base.pmax)

    for ip in 1:nip
        dict[ip] = Dict()
        nite = 1 + 2*floor(Int64, base.pmax[ip]/base.pinc[ip])

        for attr in ("slow", "maac", "bazm", "rms")
            dict[ip][attr] = Array{Float64}(undef, base.nwin)
        end

        dict[ip]["slomap"] = Array{Float64}(undef, base.nwin, nite, nite)
    end


    return dict
end


function _core(base::BaseParams)
    # initialize variables
    dict = _empty_dict(base)

    pxy0 = [0.,0.]
    for ip in 1:length(base.pmax)  
        nite = 1 + 2*floor(Int64, base.pmax[ip]/base.pinc[ip])
        pxylist = _idx(ip, nite, pxy0, base)
        
        best_maac = -1.
        best_pxy  = [0.,0.]
        for nk in 1:base.nwin
            n0  = 1 + base.toff + floor(Int64, base.lwin*base.nadv*(nk-1))
            ccmap = map(pxyl->_pccorr(n0, pxyl, base), pxylist)

            # find max
            ccmax = findmax(ccmap)
            maac  = ccmax[1]
            px = pxylist[ccmax[2]][1]
            py = pxylist[ccmax[2]][2]
            
            slow, bazm = r2p(-px, -py)
            dict[ip]["rms"][nk] = _rms(n0, [px,py], base)
            dict[ip]["slow"][nk] = slow
            dict[ip]["bazm"][nk] = bazm
            dict[ip]["maac"][nk] = maac
            dict[ip]["slomap"][nk,:,:] = reshape(ccmap, nite, nite)

            # get the best pxy
            if maac > best_maac
                best_maac = maac
                best_pxy = [px, py]
            end
        end

        # change px0 py0 and keep computing with next slowness domain only if
        # best macc > ccerr
        if best_maac > base.ccerr
            pxy0 = best_pxy
        else
            break
        end
        
    end

    return dict
end