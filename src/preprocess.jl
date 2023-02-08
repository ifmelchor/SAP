#!/usr/local/bin julia
# coding=utf-8

# Utility functions for .jl

# These functions were taken from SeisNoise julia package
# see more in github.com/tclements/SeisNoise.jl

# Ivan Melchor 2023

using FFTW

"""
    onebit!(A)

One-bit amplitude modification of Data A.
"""
function onebit!(A::AbstractArray)
  A .= sign.(A)
  return nothing
end

onebit(A::AbstractArray) = (U = deepcopy(A); onebit!(U);return U)


"""
   whiten!(A, freqmin, freqmax, fs, N; pad=50)

Whiten spectrum of rfft `A` between frequencies `freqmin` and `freqmax`.
Returns the whitened rfft of the time series.

# Arguments
- `A::AbstractArray`: Time series.
- `freqmin::Real`: Pass band low corner frequency.
- `freqmax::Real`: Pass band high corner frequency.
- `fs::Real`: Sampling rate of time series `A`.
- `N::Int`: Number of input time domain samples for each rfft.
- `pad::Int`: Number of tapering points outside whitening band.

"""
function whiten!(A::AbstractArray{Complex{T}}, freqmin::Real, freqmax::Real, fs::Real,N::Int;pad::Int=50) where T <: AbstractFloat
   Nrows,Ncols = size(A)

   # get whitening frequencies
   freqvec = FFTW.rfftfreq(N,fs)
   left = findfirst(x -> x >= freqmin, freqvec)
   right = findfirst(freqmax .<= freqvec)
   low, high = left - pad, right + pad

   if low <= 1
       low = 1
       left = low + pad
   end

   if high > length(freqvec)
       high = length(freqvec)- 1
       right = high - pad
   end

   compzero = complex(T(0))
   padarr = similar(A,T,pad)
   padarr .= T(0.):T(pad-1)
   # left zero cut-off
   A[1:low-1,:] .= compzero

   # left tapering
   A[low:left-1,:] .= cos.(T(pi) ./ T(2) .+ T(pi) ./ T(2) .* padarr ./ pad).^2 .* exp.(im .* angle.(A[low:left-1,:]))

   # pass band
   A[left:right-1,:] .= exp.(im .* angle.(A[left:right-1,:]))

   # right tapering
   A[right:high-1,:] .= cos.(T(pi) ./ T(2) .* padarr ./ pad).^2 .* exp.(im .* angle.(A[right:high-1,:]))

   # right zero cut-off
   A[high:end,:] .= compzero
   
   return nothing
end

whiten(A::AbstractArray, freqmin::Real, freqmax::Real, fs::Real, N::Int; pad::Int=50) = (U = deepcopy(A); whiten!(U,freqmin,freqmax,fs,N,pad=pad); return U)


"""
    mute(A,factor)

Set high amplitudes in array `A` to zero.
Uses median of envelope of `A` to find outliers.
"""
function mute!(A::AbstractArray,factor::Real=3)
    T = eltype(A)
    envelope = abs.(hilbert(A))
    levels = mean(envelope,dims=1)
    level = factor .* median(levels)
    A[envelope .> level] .= T(0)
    return nothing
end

mute(A::AbstractArray,factor::Real=3) = (U = deepcopy(A); mute!(U,factor); return U)


"""
    instant_phase(A::AbstractArray)

Extract instantaneous phase from signal A.

For time series `A`, its analytic representation ``S = A + H(A)``, where
``H(A)`` is the Hilbert transform of `A`. The instantaneous phase ``e^{iθ}``
of `A` is given by dividing ``S`` by its modulus: ``e^{iθ} = \\frac{S}{|S|}``
For more information on Phase Cross-Correlation, see:
[Ventosa et al., 2019](https://pubs.geoscienceworld.org/ssa/srl/article-standard/570273/towards-the-processing-of-large-data-volumes-with).
"""
function instant_phase(A::AbstractArray)
    # the analytic signal 
    s =  analytic(A)

    return s ./ abs.(s)
end

function analytic(A::AbstractArray)
    # the analytic signal 
    T = real(eltype(A))

    return hilberttransform(A) .* Complex(T(0),T(1)) .+ A
end

function hilberttransform(A::AbstractArray)
    Nrows = size(A,1)
    T = real(eltype(A))
    f = fft(A,1)
    f[1,:] .*= Complex(T(0),T(0))
    
    if iseven(Nrows)
        f[2:Nrows÷2 + Nrows % 2,:] .*= Complex(T(0),T(-1))
        f[Nrows÷2 + Nrows % 2 + 1,:] .*= Complex(T(0),T(0))
        f[Nrows÷2 + Nrows % 2 + 2: end,:] .*= Complex(T(0),T(1))
    else
        f[2:Nrows÷2 + Nrows % 2,:] .*= Complex(T(0),T(-1))
        f[Nrows÷2 + Nrows % 2 + 1 : end,:] .*= Complex(T(0),T(1))
    end

    return ifft(f,1)
end

"""
    detrend!(X::AbstractArray{<:AbstractFloat})

Remove linear trend from array `X` using least-squares regression.
"""
function detrend!(X::AbstractArray{<:AbstractFloat})
    T = eltype(X)
    N = size(X,1)

    # create linear trend matrix
    A = similar(X,T,N,2)
    A[:,2] .= T(1)
    A[:,1] .= range(T(0),T(1),length=N)
    # create linear trend matrix
    R = transpose(A) * A

    # do the matrix inverse for 2x2 matrix
    # this is really slow on GPU
    Rinv = inv(Array(R)) |> typeof(R)
    factor = Rinv * transpose(A)

    # remove trend
    X .-= A * (factor * X)

    return nothing
end

detrend(A::AbstractArray{<:AbstractFloat}) = (U = deepcopy(A); detrend!(U);return U)

"""
    taper!(A,fs; max_percentage=0.05, max_length=20.)

    Taper a time series `A` with sampling_rate `fs`.
Defaults to 'hann' window. Uses smallest of `max_percentage` * `fs`
or `max_length`.

# Arguments

- `A::AbstractArray`: Time series.
- `fs::AbstractFloat`: Sampling rate of time series `A`.
- `max_percentage::float`: Decimal percentage of taper at one end (ranging
from 0. to 0.5).
- `max_length::Real`: Length of taper at one end in seconds.

"""
function taper!(A::AbstractArray{<:AbstractFloat}, fs::Real;
                max_percentage::AbstractFloat=0.05, max_length::Real=20.)

    Nrows = size(A,1)
    wlen = min(Int(floor(Nrows * max_percentage)), Int(floor(max_length * fs)), Int(floor(Nrows/2)))
    taper_sides = hanningwindow(A, 2 * wlen)
    A[1:wlen,:] .*= taper_sides[1:wlen]
    A[end-wlen:end,:] .*= taper_sides[wlen:end]
    
    return nothing
end

taper(A::AbstractArray{<:AbstractFloat}, fs::Real; max_percentage::AbstractFloat=0.05,max_length::Real=20.) = (U = deepcopy(A); taper!(U,fs,max_percentage=max_percentage,max_length=max_length); return U)


"""
  hanningwindow(A,n)

Generate hanning window of length `n`.
Hanning window is sin(n*pi)^2.
"""
function hanningwindow(A::AbstractArray, n::Int)
   T = eltype(A)
   win = similar(A,T,n)
   win .= T(pi) .* range(T(0.),stop=T(n),length=n)
   win .= sin.(win).^2

   return win
end