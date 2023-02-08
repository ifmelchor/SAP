#!/usr/local/bin julia
# coding=utf-8

# GNU GPL v2 licenced to I. Melchor and J. Almendros 08/2022


# me gustaría añadir varias cosas
#  1: loop sobre rangos de frequencia. Esto hace que debamos filtrar directamente en julia y no en seisvo, como esta hecho. 
#  2: guardar info para poder calcular la funcion de transferencia por ventana de calculo


function CC8(data::Union{Array{T}, Nothing}, xStaUTM::Vector{T}, yStaUTM::Vector{T}, pmax::Vector{T}, pinc::Vector{T}, fsem::J, lwin::J, nwin::J, padv::T, ccerr::T, fqbands::Vector{Tuple{T,T}}) where T<:Real where J<:Int

    !isnothing(data) ||
        return nothing
    
    length(xStaUTM) == length(yStaUTM) ||
        throw(DimensionMismatch("xStaUTM and yStaUTM coordinate vectors not the same length"))

    # define baseparams
    nsta = length(xStaUTM)
    xysta = [xySta(xStaUTM[s], yStaUTM[s]) for s in 1:nsta]
    refxy = refsta(xysta)

    # array resp
    pmax_f = findmax(pmax)
    smax   = pmax_f[1]
    ds  = pinc[pmax_f[2]]
    f1  = findmin([x[1] for x in fqband])[1]
    f2  = findmax([x[2] for x in fqband])[1]
    arf = array_response(xStaUTM, yStaUTM, smax, ds, f1, f2, 0.1)

    cciter = Vector{Tuple{Int,Int}}()
    for ii in 1:nsta-1
        for jj in ii+1:nsta
            push!(cciter, (ii, jj))
        end
    end

    base = BaseParams(data, xysta, refxy[1], refxy[2], ccerr, fsem, lwin, nwin, nadv, pmax, pinc, cciter, fqbands)
    res = _core(base)

    return (res, arf)
end


function _idx(ip::J, nite::J, pxy0::Tuple{T,T}, base::BaseParams) where J<:Int where T<:Real 
    pxx = Array{Real}(undef, nite)
    pyy = Array{Real}(undef, nite)
    for i in 1:nite
        # for x
        px = pxy0[1] - base.pmax[ip] + base.pinc[ip]*(i-1)
        pxx[i] = base.pinc[ip] * floor(Int, px/base.pinc[ip])
        # for y
        py = pxy0[2] - base.pmax[ip] + base.pinc[ip]*(i-1)
        pyy[i] = base.pinc[ip] * floor(Int, py/base.pinc[ip])
    end

    return [[px, py] for px in pxx for py in pyy]
end


"""
   _pccorr(*args)

Calcula la correlacion cruzada para un slowness y banda de frequencia
"""
function _pccorr(nkk::J, pxy::Tuple{T,T}, nfq::J, base::BaseParams) where J<:Int where T<:Real
    # compute delat times for each station
    dtimes = [pxy[1]*(sta.x-base.xref) + pxy[2]*(sta.y-base.yref) for sta in base.stalist]
    
    # filter data
    fq_band = base.fqbands[nfq]

    # build cc matrix
    nsta = length(base.stalist)
    cc = zeros(Real, nsta, nsta)
    for ii in 1:nsta
        mii = nkk + floor(Int, base.fsem*dtimes[ii])
        dii = @view base.data[ii, 1+mii:base.lwin+mii]
        dii = _filt(dii, fq_band, base.fsem)
        
        for jj in ii:nsta
            mjj = nkk + floor(Int, base.fsem*dtimes[jj])
            djj = @view base.data[jj, 1+mjj:base.lwin+mjj]
            djj = _filt(djj, fq_band, base.fsem)
            cc[ii,jj] += dot(dii,djj)
        end
    end

    # computes crosscorr coefficient
    suma = sum([cc[ii,jj]/sqrt(cc[ii,ii]*cc[jj,jj]) for (ii, jj) in base.citer])
    return (2*suma + nsta) / nsta^2
end


"""
   _empty_dict(*args)

Genera un dict vacio para llenar durante el procesado.
"""
function _empty_dict(base::BaseParams)
    dict = Dict()
    
    nfq = length(base.fqbands)
    nip = length(base.pmax)

    for nfq in 1:nfq
        dict[nfq] = Dict()
        
        for ip in 1:nip
            dict[nfq][ip] = Dict()

            for attr in ("slow", "maac", "bazm")
                dict[nfq][ip][attr] = Array{Real}(undef, base.nwin)
            end

            dict[nfq][ip]["slomap"] = Array{Real}(undef, base.nwin, nite, nite)
        end

    end

    return dict
end


function _core(base::BaseParams)
    # initialize variables
    dict = _empty_dict(base)

    nro_bands = length(fqbands)
    for nfq in 1:nro_bands

        pxy0 = (0.,0.)
        for ip in 1:length(base.pmax)  
            nite = 1 + 2*floor(Int, base.pmax[ip]/base.pinc[ip])
            pxylist = _idx(ip, nite, pxy0, base)
            
            best_maac = -1.
            for nk in 1:base.nwin
                n0  = 1 + floor(Int, base.lwin*base.nadv*(nk-1))
                ccmap = map(pxyl->_pccorr(n0, pxyl, nfq, base), pxylist)

                # find max
                ccmax = findmax(ccmap)
                maac  = ccmax[1]
                px = pxylist[ccmax[2]][1]
                py = pxylist[ccmax[2]][2]
                
                slow, bazm = r2p(-px, -py)
                dict[nfq][ip]["slow"][nk] = slow
                dict[nfq][ip]["bazm"][nk] = bazm
                dict[nfq][ip]["maac"][nk] = ccmax
                dict[nfq][ip]["slomap"][nk,:,:] = reshape(ccmap, nite, nite)

                # get the best pxy
                if ccmax > best_maac
                    best_maac = ccmax
                    best_pxy = (px, py)
                end
            end

            # change px0 py0 and keep computing with next slowness domain only if
            # best macc > ccerr
            if best_maac > base.ccerr
                pxy0 = best_pxy
            else
                last_ip = ip
                break
            end
            
            last_ip = ip
        end

        if last_ip < length(fqbands)
            for ip in last_ip:length(fqbands)
                dict[nfq][ip] = nothing
            end
        end

    end

    return dict
end