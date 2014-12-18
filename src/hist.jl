using StatsBase

import Base.ndims
import Base.push!, Base.+

nbins(h::Histogram) = map(x->size(x, 1) - 1, edges(h))
edges(h::Histogram) = h.edges

contents(h::Histogram) = h.weights

integral(h::Histogram) = sum(h.weights)

type ErrorHistogram{T<:Real, N, E}
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
push!{T,E}(h::ErrorHistogram{T,1,E}, x::Real) = push!(h,x,one(T))

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
push!{T,N}(h::ErrorHistogram{T,N},xs::NTuple{N,Real}) = push!(h,xs, one(T))

nbins(h::ErrorHistogram) = map(x->size(x, 1) - 1, edges(h))
ndims{T, N, E}(h::ErrorHistogram{T, N, E}) = N

errors(h::ErrorHistogram) =  sqrt(h.weights_sq)
contents(h::ErrorHistogram) = h.weights
edges(h::ErrorHistogram) = h.edges

normalize(h::Histogram) = Histogram(h.edges, h.weights/sum(h.weights))

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

export nbins, edges, contents, errors, integral
export push!
export ErrorHistogram, Hist1D, Hist2D, Hist3D
export +