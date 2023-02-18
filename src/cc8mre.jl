#!/usr/local/bin julia
# coding=utf-8

# GNU GPL v2 licenced to I. Melchor and J. Almendros 08/2022


function CC8(data::Array{T}, xStaUTM::Vector{T}, yStaUTM::Vector{T}, pmax::Vector{T}, pinc::Vector{T}, fsem::J, lwin::J, nwin::J, nadv::T, ccerr::T, toff::J) where T<:Real where J<:Integer
    
    # define delta time function
    dtime = _dtimefunc(xStaUTM, yStaUTM, fsem)

    # define base params
    nsta   = length(xStaUTM)
    cciter = _cciter(nsta)
    nites  = _nites(pmax, pinc)
    toff   = floor(Int64, toff*fsem)
    base   = BaseParams(nites, nwin, nsta, lwin, cciter)

    # initialize variables
    dict = _empty_dict(base)

    # init pxy
    pxy0::Vector{Float64} = [0.,0.]
    
    # loop over slowness domain
    for ip in 1:length(nites)
        # create slowness grid
        pxy_map  = _pxymap(pxy0, nites[ip], pinc[ip], pmax[ip])

        # create delta times grid
        time_map = _dtimemap(dtime, pxy_map, nsta)
        
        # init params
        best_maac::Float64 = -1.
        best_pxy::Vector{Float64} = [0., 0.]

        # iterate over time
        for nk in 1:nwin
            n0  = 1 + toff + floor(Int64, lwin*nadv*(nk-1)) 

            # get ccmap
            ccmap = _ccmap(data, n0, nites[ip], time_map, base)
            
            # find max value
            ccmax = findmax(ccmap)
            maac  = ccmax[1]
            (ii, jj) = ccmax[2].I
            dict[ip]["maac"][nk] = maac  # save max correlation coef
            
            best_pxy_n = pxy_map[ii, jj, :]
            slow, bazm = r2p(-1 .* best_pxy_n)
            dict[ip]["slow"][nk] = slow
            dict[ip]["bazm"][nk] = bazm
            dict[ip]["rms"][nk]  = _rms(data, n0, time_map[ii, jj, :], base)  # save rms of the best ccorr
            dict[ip]["slowmap"][nk,:,:] = ccmap
           
            # get bounds
            bounds = bm2(ccmap, pmax[ip], pinc[ip], maac, ccerr)
            dict[ip]["slowbnd"][nk,:] = [bounds.slomin, bounds.slomax]
            dict[ip]["bazmbnd"][nk,:] = [bounds.azimin, bounds.azimax]
            
            if maac > best_maac
                best_maac = maac
                best_pxy = best_pxy_n
            end
        end

        
        if best_maac > ccerr
            pxy0 = best_pxy
        end
    end

    return dict
end


function _pxymap(pxy0::Vector{T}, nite::J, pinc::T, pmax::T) where {T<:Real, J<:Integer}
    pxy_map = Array{Float64}(undef, nite, nite, 2)
    
    for ii in 1:nite # for x
        px = pxy0[1] - pmax + pinc*(ii-1)
        pxi = pinc * floor(Int64, px/pinc)
        for jj in 1:nite # for y
            py = pxy0[2] - pmax + pinc*(ii-1)
            pyj = pinc * floor(Int64, py/pinc)
            pxy_map[ii,jj,:] = [pxi,pyj]
        end
    end
    
    pxy_map[:,:,2] = adjoint(pxy_map[:,:,2])
    return pxy_map
end


function _dtimemap(dtime_func::Function, pxy_map::Array{T}, nsta::J) where {T<:Real, J<:Integer}
    nite = size(pxy_map, 1)
    time_map = Array{Int}(undef, nite, nite, nsta)
    for ii in 1:nite
        for jj in 1:nite
            time_map[ii,jj,:] = dtime_func(pxy_map[ii,jj,:])
        end
    end

    return time_map
end


function _pccorr(data::Array{T}, nkk::J, pxytime::Vector{J}, base::BaseParams) where {T<:Real, J<:Integer}
    cc = zeros(Float64, base.nsta, base.nsta)
    for ii in 1:base.nsta
        mii = nkk + pxytime[ii]
        dii = @view data[ii, mii:base.lwin+mii]
        for jj in ii:base.nsta
            mjj = nkk + pxytime[jj]
            djj = @view data[jj, mjj:base.lwin+mjj]
            cc[ii,jj] += dot(dii,djj)
        end
    end
    # computes crosscorr coefficient
    suma = 2*sum([cc[ii,jj]/sqrt(cc[ii,ii]*cc[jj,jj]) for (ii, jj) in base.citer])
    return (suma+base.nsta) / (base.nsta*base.nsta)
end


function _ccmap(data::Array{T}, n0::J, nite::J, time_map::Array{J}, base::BaseParams) where {T<:Real, J<:Integer}
    cc_map = zeros(Float64, nite, nite)
    
    for ii in 1:nite # for x
       for jj in 1:nite # for y
          cc_map[ii,jj] = _pccorr(data, n0, time_map[ii,jj,:], base)
       end
    end
    
    return cc_map
end


function _rms(data::Array{T}, nkk::J, pxytime::Vector{J}, base::BaseParams) where J<:Integer where T<:Real

    erg = 0.
    for ii in 1:base.nsta
        mii = nkk + pxytime[ii]
        dii = @view data[ii, 1+mii:base.lwin+mii]
        erg += sqrt(mean(dii.*dii))
    end
    
    return erg /= base.nsta
end