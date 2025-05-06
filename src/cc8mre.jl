#!/usr/local/bin julia
# coding=utf-8

# GNU GPL v2 licenced to I. Melchor and J. Almendros 08/2022

# Main ZLCC codes

function zlcc(data::Array{T}, xStaUTM::Array{T}, yStaUTM::Array{T}, slomax::T, sloint::T, fqband::Vector{T}, fsem::J, lwin::J, nwin::J, nadv::T, toff::J, slow0::Vector{T}=[0., 0.], ccerr::T=0.95, slow2::Bool=True, maac_thr::T=0.6, slomax2::T=0.3, sloint2::T=0.02) where {T<:Real, J<:Integer}
    
    # define base params
    nsta   = length(xStaUTM) # nro of stations
    cciter = _cciter(nsta)   # stations iterator
    toff   = toff*fsem       # off seconds in sp
    nite   = 1 + 2*round(Int64, slomax/sloint)   # slowness grid

    if slow2
        nite2  = 1 + 2*round(Int64, slomax2/sloint2) # slowness grid2
    end

    # base object for crosscorrelation
    base   = Base(nite, nwin, nsta, lwin, cciter, slow2)

    # compute base error of the array response
    # arresp = array_transfunc(xStaUTM, yStaUTM, slomax, sloint, fqband[1], fqband[2], 0.05)
    # powar = reshape(arresp.power,1,:)
    # base_error = 100*sum(powar .> ccerr)/(nite*nite)

    # init empty dictionary
    dict = _empty_dict(base)
    
    # create slowness main grid
    slow_grid  = _xygrid(slow0, sloint, slomax)

    # create deltatimes grid
    dtime     = _dtimefunc(xStaUTM, yStaUTM, fsem) # define delta time function
    time_grid = _dtimemap(dtime, slow_grid, nsta)

    # filter data
    _filter!(data, fsem, fqband)
    
    # iterate over time
    @inbounds for nk in 1:nwin
        n0  = 1 + toff + lwin*nadv*(nk-1)

        # get ccmap
        ccmap = _ccmap(data, n0, time_grid, nite, base)

        # find max value MAAC and position
        ccmax      = findmax(ccmap)
        maac       = ccmax[1]
        (ii, jj)   = ccmax[2].I
        best_slow  = slow_grid[ii, jj, :]

        # get slow, baz, and rms
        slow, bazm = r2p(-1 .* best_slow)
        bounds     = _bounds(ccmap, slow_grid, maac*ccerr)
        rms        = _rms(data, n0, time_grid[ii, jj, :], base)

        # compue error [%]
        npts  = sum(reshape(ccmap,1,:) .> maac*ccerr)
        error = 100*npts/(nite*nite)

        # compute slow2 onyl when maac is above threshold
        if slow2 && maac > maac_thr
            slow_grid2 = _xygrid(best_slow, sloint2, slomax2)
            time_grid2 = _dtimemap(dtime, slow_grid2, nsta)
            ccmap2 = _ccmap(data, n0, time_grid2, nite2, base)
            ccmax2 = findmax(ccmap2)
            # update maac
            maac = ccmax2[1]
            (ii, jj)    = ccmax2[2].I
            # update rms
            rms  = _rms(data, n0, time_grid2[ii, jj, :], base)
            best_slow2  = slow_grid2[ii, jj, :]
            slo2, baz2  = r2p(-1 .* best_slow2)
        else
            slo2 = NaN
            baz2 = NaN
        end

        # save values into dict
        dict["maac"][nk]   = maac
        dict["rms"][nk]    = rms
        dict["error"][nk]  = error
        dict["slowmap"][nk,:,:] = ccmap
        dict["slow"][nk]   = slow
        dict["baz"][nk]    = bazm
        dict["slowbnd"][nk,:] = [bounds.slomin, bounds.slomax]
        dict["bazbnd"][nk,:]  = [bounds.azimin, bounds.azimax]
        
        if slow2
            dict["slow2"][nk]  = slo2
            dict["baz2"][nk]   = baz2
        end
    
    end

    return dict
end


function _pccorr(data::Array{T}, nkk::T, pxytime::Vector{T}, base::Base) where T<:Real
    cc = zeros(T, base.nsta, base.nsta)
    
    for ii in 1:base.nsta
        mii = round(Int32, nkk + pxytime[ii])
        dii = @view data[ii, mii:base.lwin+mii]
        for jj in ii:base.nsta
            mjj = round(Int32, nkk + pxytime[jj])
            djj = @view data[jj, mjj:base.lwin+mjj]
            cc[ii,jj] += dot(dii,djj)
        end
    end

    # computes crosscorr coefficient
    cc_sum = 2*sum([cc[ii,jj]/sqrt(cc[ii,ii]*cc[jj,jj]) for (ii, jj) in base.citer])
    
    return (cc_sum+base.nsta) / (base.nsta*base.nsta)
end


function _ccmap(data::Array{T}, n0::T, time_map::Array{T}, nite::J, base::Base) where {T<:Real, J<:Integer}

    cc_map = zeros(T, nite, nite)
    
    @inbounds for ii in 1:nite, jj in 1:nite
        cc_map[ii,jj] = _pccorr(data, n0, time_map[ii,jj,:], base)
    end
    
    return cc_map
end


function _rms(data::Array{T}, nkk::T, pxytime::Vector{T}, base::Base) where T<:Real

    erg = 0.
    for ii in 1:base.nsta
        mii = round(Int32, nkk + pxytime[ii])
        dii = @view data[ii, 1+mii:base.lwin+mii]
        erg += sqrt(mean(dii.^2))
    end
    
    return erg /= base.nsta
end