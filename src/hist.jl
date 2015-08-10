using StatsBase

import Base.ndims
import Base.push!, Base.+

nbins(h::Histogram) = map(x->size(x, 1) - 1, edges(h))
edges(h::Histogram) = h.edges

contents(h::Histogram) = h.weights

integral(h::Histogram) = sum(h.weights)

type ErrorHistogram{T<:Real, N, E} <: AbstractHistogram{T, N, E}
    edges::E
    weights::Array{T,N}
    weights_sq::Array{T,N}
    closed::Symbol
    function ErrorHistogram(
        edges::NTuple{N,AbstractArray},
        weights::Array{T,N},
        weights_sq::Array{T,N},
        closed::Symbol)
        closed == :right || closed == :left || error("closed must :left or :right")
        map(x -> length(x)-1,edges) == size(weights) || error("Histogram edge vectors must be 1 longer than corresponding weight dimensions")
        new(edges, weights, weights_sq, closed)
    end 
end

Hist1D = ErrorHistogram{Float64, 1, (Vector{Float64}, )}
Hist2D = ErrorHistogram{Float64, 2, (Vector{Float64}, Vector{Float64}) }
Hist3D = ErrorHistogram{Float64, 2, (Vector{Float64}, Vector{Float64}, Vector{Float64})}

function ErrorHistogram(edges...)
    T = Float64
    N = length(edges)
    E = typeof(edges)

    return ErrorHistogram{T, N, E}(
        edges,
        zeros(T, map(x->size(x, 1) - 1, edges)...),
        zeros(T, map(x->size(x, 1) - 1, edges)...),
        :left
    )
end

#from StatsBase
function push!{T,E}(h::ErrorHistogram{T,1,E}, x::Real, w::Real)
    i = if h.closed == :right 
        searchsortedfirst(edges(h)[1], x) - 1 
    else
        searchsortedlast(edges(h)[1], x)
    end

    if 1 <= i <= length(h.weights)
        @inbounds h.weights[i] += w
        @inbounds h.weights_sq[i] += w^2
    end
    h
end

function push!{T,N}(h::ErrorHistogram{T,N},xs::NTuple{N,Real}, w::Real)
    is = if h.closed == :right
        map((edge, x) -> searchsortedfirst(edge,x) - 1, h.edges, xs)
    else
        map(searchsortedlast, h.edges, xs)
    end
    try
        h.weights[is...] += w
        h.weights_sq[is...] += w^2
    catch e
        isa(e,BoundsError) || rethrow(e)
    end
    h
end

nbins(h::ErrorHistogram) = map(x->size(x, 1) - 1, edges(h))
ndims{T, N, E}(h::ErrorHistogram{T, N, E}) = N

errors(h::ErrorHistogram) =  sqrt(h.weights_sq)
contents(h::ErrorHistogram) = h.weights
edges(h::ErrorHistogram) = h.edges

normalize(h::Histogram) = Histogram(h.edges, h.weights/sum(h.weights))
function normalize{T<:Real, N, E}(h::ErrorHistogram{T, N, E})
    sw = sum(h.weights)
    ErrorHistogram{T, N, E}(h.edges, h.weights/sw, h.weights_sq, h.closed)
end

function +{T1<:Real, T2<:Real, N, E}(h1::ErrorHistogram{T1, N, E}, h2::ErrorHistogram{T2, N, E})
    @assert ndims(h1) == ndims(h2)
    @assert nbins(h1) == nbins(h2)
    @assert h1.closed == h2.closed

    for i=1:ndims(h1)
        @assert all(edges(h1)[i] .== edges(h2)[i])
    end

    return ErrorHistogram{T1, N, E}(
        edges(h1),
        contents(h1) .+ contents(h2),
        h1.weights_sq .+ h2.weights_sq,
        h1.closed
    )
end

integral(h::ErrorHistogram) = (sum(contents(h)), sqrt(sum(h.weights_sq)))

import Base.mean
function mean(h::Hist1D)
    l = 0.0
    mids = midpoints(h.edges[1][2:end-1])
    for i=2:nbins(h)[1]-1
        l += mids[i-1] * h.weights[i]
    end
    return l / sum(h.weights)
end

import Base.std
function std(h::Hist1D)
    l = 0.0
    mids = midpoints(h.edges[1][2:end-1])
    m = mean(h)
    for i=2:nbins(h)[1]-1
        l += (mids[i-1] - m)^2 * h.weights[i] 
    end
    return sqrt(l / sum(h.weights))
end

export nbins, edges, contents, errors, integral
export push!
export normalize
export ErrorHistogram, Hist1D, Hist2D, Hist3D
export +
