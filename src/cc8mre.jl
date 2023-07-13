#!/usr/local/bin julia
# coding=utf-8

# GNU GPL v2 licenced to I. Melchor and J. Almendros 08/2022


function CC8(data::Array{T}, xStaUTM::Array{T}, yStaUTM::Array{T}, slomax::T, sloint::T, fqband::Vector{T}, fsem::J, lwin::J, nwin::J, nadv::T, ccerr::T, toff::J; slow0::Vector{T}=[0.,0.]) where {T<:Real, J<:Integer}
    
    # filter data
    _filter!(data, fsem, fqband)

    # define base params
    nsta   = length(xStaUTM) # nro of stations
    cciter = _cciter(nsta)   # stations iterator
    toff   = toff*fsem       # off seconds in sp
    nite   = 1 + 2*round(Int64, slomax/sloint) # slowness grid 
    base   = Base(nite, nwin, nsta, lwin, cciter) # base object

    # init empty dictionary
    dict = _empty_dict(base)
    
    # create slowness grid
    slow_grid  = _xygrid(slow0, sloint, slomax)

    # create deltatimes grid
    dtime     = _dtimefunc(xStaUTM, yStaUTM, fsem) # define delta time function
    time_grid = _dtimemap(dtime, slow_grid, nsta)
    
    # iterate over time
    for nk in 1:nwin
        n0  = 1 + toff + lwin*nadv*(nk-1)

        # get ccmap
        ccmap = _ccmap(data, n0, time_grid, base)
            
        # find max value
        ccmax      = findmax(ccmap)
        maac       = ccmax[1]
        (ii, jj)   = ccmax[2].I
        best_pxy_n = pxy_map[ii, jj, :]
        slow, bazm = r2p(-1 .* best_pxy_n)
        rms        = _rms(data, n0, time_grid[ii, jj, :], base)
        bounds     = bm2(ccmap, slomax, sloint, maac, ccerr)

        # save values into dict
        dict["maac"][nk] = maac 
        dict["slow"][nk] = slow
        dict["bazm"][nk] = bazm
        dict["rms"][nk]  = rms
        dict["slowmap"][nk,:,:] = ccmap
        dict["slowbnd"][nk,:] = [bounds.slomin, bounds.slomax]
        dict["bazmbnd"][nk,:] = [bounds.azimin, bounds.azimax]
    
    end

    return dict
end


"""
   get_dtimes(*args)

Devuelve los delta times correspondientes a un slowness y un azimuth
"""
function get_dtimes(slow::T, baz::T, slomax::T, sloint::T, xStaUTM::Array{T}, yStaUTM::Array{T}, etol::T; slow0::Vector{T}=[0.,0.]) where T<:Real

    # nro stations
    nsta   = length(xStaUTM)
    
    # create nites
    nite = 1 + 2*round(Int64, slomax/sloint)

    # create slowness grid
    slow_grid  = _xygrid(slow0, sloint, slomax)
    dtime      = _dtimefunc(xStaUTM, yStaUTM, fsem) # define delta time function
    time_grid  = _dtimemap(dtime, slow_grid, nsta)

    ijmin = [1, 1]
    tol    = [999., 999.]
  
    for ii in 1:nite, jj in 1:nite
        pxy = slow_grid[ii,jj,:]
        slow, baz = r2p(-1 .* pxy)
        slodif = abs(slow-slowness)
        bazdif = abs(baz-bazimuth)

        if slodif < etol && bazdif < etol
            tol[1] = slodif
            tol[2] = bazdif
            ijmin = [ii, jj]
            break
        end
    end

    ii, jj = ijmin
    deltas = time_grid[ii, jj, :]

    return deltas, tol
  
end


function _xygrid(slow0::Vector{T}, sloint::T, slomax::T) where T<:Real
    #
    # This function cretes the slownes grid
    #

    # define the size of the grid
    nite    = 1 + 2*round(Int64, slomax/sloint)
    
    # init the grid in memeory
    xy_grid = Array{T}(undef, nite, nite, 2)
    
    # fill the grid
    for ii in 1:nite, jj in 1:nite
        px = slow0[1] - slomax + sloint*(ii-1)
        # pxi = pinc * px/pinc
        py = slow0[2] - slomax + sloint*(ii-1)
        # pyj = pinc * py/pinc
        xy_grid[ii,jj,:] = [px,py]
    end
    
    xy_grid[:,:,2] = adjoint(xy_grid[:,:,2])
    
    return xy_grid
end


function _dtimemap(dtime_func::Function, pxy_map::Array{T}, nsta::J) where {T<:Real, J<:Integer}
    nite = size(pxy_map, 1)
    
    time_map = Array{T}(undef, nite, nite, nsta)
    for ii in 1:nite, jj in 1:nite 
        time_map[ii,jj,:] = dtime_func(pxy_map[ii,jj,:])
    end

    return time_map
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


function _ccmap(data::Array{T}, n0::T, time_map::Array{T}, base::Base) where {T<:Real, J<:Integer}
    cc_map = zeros(T, base.nite, base.nite)
    
    for ii in 1:nite, jj in 1:nite
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