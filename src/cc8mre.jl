#!/usr/local/bin julia
# coding=utf-8

# GNU GPL v2 licenced to I. Melchor and J. Almendros 08/2022

# Main ZLCC codes

function CC8(data::Array{T}, xStaUTM::Array{T}, yStaUTM::Array{T}, slomax::T, sloint::T, fqband::Vector{T}, fsem::J, lwin::J, nwin::J, nadv::T, ccerr::T, toff::J, slow0::Vector{T}) where {T<:Real, J<:Integer}
    
    # define base params
    nsta   = length(xStaUTM) # nro of stations
    cciter = _cciter(nsta)   # stations iterator
    toff   = toff*fsem       # off seconds in sp
    nite   = 1 + 2*round(Int64, slomax/sloint) # slowness grid
    base   = Base(nite, nwin, nsta, lwin, cciter) # base object for crosscorrelation

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
    for nk in 1:nwin
        n0  = 1 + toff + lwin*nadv*(nk-1)

        # get ccmap
        ccmap = _ccmap(data, n0, time_grid, base)

        # find max value MAAC and position
        ccmax      = findmax(ccmap)
        maac       = ccmax[1]
        (ii, jj)   = ccmax[2].I
        best_slow  = slow_grid[ii, jj, :]

        # get slow, baz, and rms
        slow, bazm = r2p(-1 .* best_slow)
        bounds = _bounds(ccmap, slow_grid, maac*0.95)
        rms = _rms(data, n0, time_grid[ii, jj, :], base)

        # compue error [%]
        npts  = sum(reshape(ccmap,1,:) .> maac*ccerr)
        error = 100*npts/(nite*nite)

        # save values into dict
        dict["maac"][nk]   = maac
        dict["rms"][nk]    = rms
        dict["error"][nk]  = error
        dict["slowmap"][nk,:,:] = ccmap
        dict["slow"][nk] = slow
        dict["baz"][nk]  = bazm
        dict["slowbnd"][nk,:] = [bounds.slomin, bounds.slomax]
        dict["bazbnd"][nk,:] = [bounds.azimin, bounds.azimax]
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


function _ccmap(data::Array{T}, n0::T, time_map::Array{T}, base::Base) where T<:Real
    cc_map = zeros(T, base.nite, base.nite)
    
    for ii in 1:base.nite, jj in 1:base.nite
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