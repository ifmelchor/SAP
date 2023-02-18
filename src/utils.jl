#!/usr/local/bin julia
# coding=utf-8

# Utility functions for cc8mre.jl

# GNU GPL v2 licenced to I. Melchor and J. Almendros 08/2022


"""
   _empty_dict(*args)

Genera un dict vacio para llenar durante el procesado.
"""
function _empty_dict(base::BaseParams)
    dict = Dict()
    
    for ip in 1:length(base.nites)
        dict[ip] = Dict()
        nite = base.nites[ip]

        for attr in ("slow", "maac", "bazm", "rms")
            dict[ip][attr] = Array{Float64}(undef, base.nwin)
        end

        dict[ip]["slowmap"] = Array{Float64}(undef, base.nwin, nite, nite)
        dict[ip]["slowbnd"] = Array{Float64}(undef, base.nwin, 2)
        dict[ip]["bazmbnd"] = Array{Float64}(undef, base.nwin, 2)
    end


    return dict
end

"""
   _dtimefunc(*args)

Genera la función que devuelve los delta times para un vector de lentidud aparente
"""
function _dtimefunc(stax::Vector{T}, stay::Vector{T}, fsem::J) where {T<:Real, J<:Integer}
    xref = mean(xStaUTM)
    yref = mean(yStaUTM)
    dtime(pxy) = floor.(Int, [pxy[1]*(stx-xref) + pxy[2]*(sty-yref) for (stx, sty) in zip(stax,stay)] .* fsem)

    return dtime
end


function _cciter(nsta::J) where J<:Integer
  
  cciter = Vector{Tuple{Int,Int}}()
  for ii in 1:nsta-1
      for jj in ii+1:nsta
          push!(cciter, (ii, jj))
      end
  end

  return cciter
end


"""
   _nites(*args)

Genera una lista con el número de intervalos de lentitud aparente
"""
function _nites(pmax::Vector{T}, pinc::Vector{T})
  
  nites = [1 + 2*floor(Int64, i/j) for (i, j) in zip(pmax, pinc)]

  return nites
end


"""
  r2p(x, y)
    
    Get slowness and back-azimuth angles
    
"""

function r2p(pxy::Vector{T}) where T<:Real
  x = pxy[1]
  y = pxy[2]

  slowness = sqrt(x^2 + y^2)
  azimuth = 180.
  
  if y < 0
    azimuth = 180. + 180*atan(x/y)/pi
  
  elseif y > 0
    azimuth = 180*atan(x/y)/pi
    if x < 0
      azimuth += 360.
    end
  
  else # y == 0

    if x > 0
      azimuth = 90.
    
    elseif x < 0
      azimuth = 270.
    
    else # x == 0
      azimuth = 666.
    
    end
  
  end
  
  return (slowness, azimuth)
end

"""
  r2p(x, y)
    
    Get slowness and back-azimuth bounds
    
"""
function bm2(msum::AbstractArray{T}, pmax::T, pinc::T, ccmax::T, ccerr::T) where T<:Real
  nite = size(msum, 1)
  bnd = Bounds(666., -1., 666., -1.)
  q = Array{Float64}(undef, nite, nite)

  ccth = ccmax - ccerr

  for i in 1:nite
    px = pinc * (i-1) - pmax  
    
    for j in 1:nite 
      py = pinc * (j-1) - pmax 
      
      ( (px == 0) && (py == 0) ) && continue # skip for

      if msum[i,j] >= ccth
        q[i,j] = 1
        
        for x in (-px+pinc, -px-pinc)
          for y in (-py+pinc, -py-pinc)
            s, a = r2p([x, y])
            if ( s > bnd.slomax ) ; bnd.slomax = s end
            if ( s < bnd.slomin ) ; bnd.slomin = s end
            if ( a > bnd.azimax ) ; bnd.azimax = a end
            if ( a < bnd.azimin ) ; bnd.azimin = a end
          end
        end
      
      else
        q[i,j] = 0
      
      end

    end

  end

  if (bnd.azimax > 355) && (bnd.azimin < 5)
    bnd.azimin = 666.
    bnd.azimax = 1.
    
    for i in 1:nite
      px = pinc * (i-1) - pmax 

      for j in 1:nite
        py = pinc * (j-1) - pmax

        ( (px == 0) && (py == 0) ) && continue # skip for
        
        if convert(Bool, q[i,j])
          for x in [-px+pinc, -px-pinc]
            for y in [-py+pinc, -py-pinc]
              s, a = r2p([x, y])
              if ( x > 0 ) && ( a > bnd.azimax ) ; bnd.azimax = a end
              if ( x < 0 ) && ( a < bnd.azimin ) ; bnd.azimin = a end
            end
          end
        end

      end
    end

  end

  return bnd
end



