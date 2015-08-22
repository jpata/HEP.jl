module HEP

import Base.convert
import Base.eta

import Base.+
import Base.*

#A 4-tuple representing a 4 vector in the natural system of units
abstract AbstractFourVector

#Cartesian
immutable FourVector <: AbstractFourVector
    t::Float64
    x::Float64
    y::Float64
    z::Float64
end

#Spherical
immutable FourVectorSph <: AbstractFourVector
    perp::Float64 #perpendicular component (e.g. in case of 4-momentum -> pt or transverse momentum)
    eta::Float64
    phi::Float64
    l::Float64 #length
end

+(v1::FourVector, v2::FourVector) = FourVector(v1.t+v2.t, v1.x+v2.x, v1.y+v2.y, v1.z+v2.z)
*(v1::FourVector, a::FourVector) = FourVector(a * v1.t, a * v1.x, a * v1.y, a*v1.z)
-(v1::FourVector, v2::FourVector) = v1 + -1.0 * v2

#magnitude
l2(v::FourVector) = v.t^2 - v.x^2 - v.y^2 - v.z^2
function l(v::FourVector)
    sq = l2(v)
    return sq > 0 ? sqrt(sq) : -sqrt(-sq)
end

#polar
phi(v::FourVector) = atan2(v.y, v.x)
rho(v::FourVector) = sqrt(v.x^2 + v.y^2 + v.z^2)
theta(v::FourVector) = atan2(perp(v), v.z)

#pseudorapidity
function eta(v::FourVector)
    ct = cos(theta(v))
    return ct^2 < 1 ? -0.5 * log((1.0-ct)/(1.0+ct)) : sign(v.z)*Inf
end

#transverse component (perpendicular)
perp(v::FourVector) = sqrt(v.x^2 + v.y^2)

phi(v::FourVectorSph) = v.phi
eta(v::FourVectorSph) = v.eta
perp(v::FourVectorSph) = v.perp
l(v::FourVectorSph) = v.l
rho(v::FourVectorSph) = rho(convert(v, FourVector))
theta(v::FourVectorSph) = theta(convert(v, FourVector))

x(v::FourVector) = v.x
y(v::FourVector) = v.y
z(v::FourVector) = v.z
t(v::FourVector) = v.t

x(v::FourVectorSph) = perp(v) * cos(phi(v))
y(v::FourVectorSph) = perp(v) * sin(phi(v))
z(v::FourVectorSph) = perp(v) * sinh(eta(v))
t(v::FourVectorSph) = sqrt(sign(l(v)) * l(v)^2 + perp(v)^2 * (1 + sinh(eta(v))^2))

# #spherical components
# function FourVectorSph(perp, eta, phi, l)
#     perp = abs(perp)
#     x = perp * cos(phi)
#     y = perp * sin(phi)
#     z = perp * sinh(eta)
#     t = l>0 ? sqrt(l^2+x^2+y^2+z^2) : sqrt(-l^2+x^2+y^2+z^2)
#     return FourVector(t, x, y, z)
# end

function convert(a::FourVectorSph, b::Type{FourVector})
    FourVector(t(a), x(a), y(a), z(a))
end

function convert(a::FourVector, b::Type{FourVectorSph})
    FourVector(perp(a), eta(a), phi(a), l(a))
end

#calculates the delta R between 2 four momenta
function deltar{T1 <: AbstractFourVector, T2 <: AbstractFourVector}(v1::T1, v2::T2)
    deta = eta(v1) - eta(v2)
    dphi = phi(v1) - phi(v2)
   
    #bring dphi to -pi...+pi
    while dphi >= pi
        dphi -= 2 * pi
    end
    
    while dphi < -pi 
        dphi += 2 * pi
    end

    return sqrt(deta^2 + dphi^2)
end

export FourVector, FourVectorSph
export l2, l, phi, rho, theta, eta, perp
export x, y, z, t
export deltar
export +, -, *

include("hist.jl")

#if ROOT is available, create additional ROOT histogram conversion methods
if isfile(joinpath(Pkg.dir("ROOT"), "src", "ROOT.jl")) && haskey(ENV, "ROOTSYS")
using ROOT

function to_jl(o::TH1D; error_type=:errors)
    
    nb = int64(GetNbinsX(o))
    #+3 = underflow, overflow, superfluous extra bin
    conts = zeros(Float64, nb+2)
    errs = zeros(Float64, nb+2)
    ents = zeros(Float64, nb+2)

    #underflow low, overflow low, overflow high
    edges = zeros(Float64, nb+3)
    for n=0:nb+1
        conts[n+1] = GetBinContent(o, int32(n))|>float64
        errs[n+1] = GetBinError(o, int32(n))|>float64
        #entries[n+1] = GetEntries(h) * conts[n+1]
        edges[n+1] = GetBinLowEdge(o, int32(n))
    end
    
    if error_type == :errors
        ents = (conts ./ errs).^2
        ents[isnan(ents)] = 0.0
        ents[ents .== Inf] = 1.0
        #println(hcat(conts, errs, ents, conts./sqrt(ents)))
    end

    #this works for Poisson bin errors
    if error_type == :contents
        ents = conts ./ sum(conts) .* float64(GetEntries(o))
        ents[isnan(ents)] = 0.0
    end
    
    edges[1] = -Inf
    edges[nb+2] = edges[nb+1] + GetBinWidth(o, int32(nb))
    edges[nb+3] = Inf
    
    return (edges,), conts, errs
    
end

function to_jl(o::TH2D)
    
    nx = int64(GetNbinsX(o))
    ny = int64(GetNbinsY(o))
    
    arr = zeros(Float64, nx+2, ny+2)
    errs = zeros(Float64, nx+2, ny+2)
    
    edges_x = zeros(Float64, nx + 3)
    edges_y = zeros(Float64, ny + 3)
    
    for n=0:nx+1
        #edges_x[n+1] = GetBinLowEdge(o, int32(n))
        edges_x[n+1] = ROOT.GetBinLowEdgeX(o, int32(n))
    end
    for n=0:ny+1
        #edges_x[n+1] = GetBinLowEdge(o, int32(n))
        edges_y[n+1] = ROOT.GetBinLowEdgeY(o, int32(n))
    end
    for x=0:nx+1
        for y=0:ny+1
            arr[x+1, y+1] = GetBinContent(root_cast(TH1, o), int32(x), int32(y))
            errs[x+1, y+1] = GetBinError(root_cast(TH1, o), int32(x), int32(y))
        end
    end
    edges_x[1] = -Inf
    edges_x[end] = +Inf

    edges_y[1] = -Inf
    edges_y[end] = +Inf
    
    ents = arr .* GetEntries(o) ./ sum(arr)
    ents[isnan(ents)] = 0
    ents = round(ents)
    return (edges_x, edges_y), arr, errs
end

function to_jl(o::TH3D)
    
    nx = int64(GetNbinsX(o))
    ny = int64(GetNbinsY(o))
    nz = int64(GetNbinsZ(o))
    
    arr = zeros(Float64, nx+2, ny+2, nz+2)
    errs = zeros(Float64, nx+2, ny+2, nz+2)
    
    edges_x = zeros(Float64, nx + 3)
    edges_y = zeros(Float64, ny + 3)
    edges_z = zeros(Float64, nz + 3)
    
    for n=0:nx+1
        #edges_x[n+1] = GetBinLowEdge(o, int32(n))
        edges_x[n+1] = ROOT.GetBinLowEdgeX(o, int32(n))
    end
    for n=0:ny+1
        #edges_x[n+1] = GetBinLowEdge(o, int32(n))
        edges_y[n+1] = ROOT.GetBinLowEdgeY(o, int32(n))
    end
    for n=0:nz+1
        #edges_x[n+1] = GetBinLowEdge(o, int32(n))
        edges_z[n+1] = ROOT.GetBinLowEdgeZ(o, int32(n))
    end
    for x=0:nx+1
        for y=0:ny+1
            for z=0:nz+1
                arr[x+1, y+1, z+1] = GetBinContent(root_cast(TH1, o), int32(x), int32(y), int32(z))
                errs[x+1, y+1, z+1] = GetBinError(root_cast(TH1, o), int32(x), int32(y), int32(z))
            end
        end
    end
    edges_x[1] = -Inf
    edges_x[end] = +Inf

    edges_y[1] = -Inf
    edges_y[end] = +Inf

    edges_z[1] = -Inf
    edges_z[end] = +Inf
    
    ents = arr .* GetEntries(o) ./ sum(arr)
    ents[isnan(ents)] = 0
    ents = round(ents)
    return (edges_x, edges_y, edges_z), arr, errs
end

function read_file(fn)
    tf = TFile(fn)
    kl = GetListOfKeys(tf)
    
    hd = Dict()
    for i=1:length(kl)
        k = root_cast(TKey, kl[i])
        kn = k|>GetName|>bytestring
        o = ReadObj(k);
        clname = k|>GetClassName|>bytestring|>symbol
        
        o = root_cast(eval(clname), o)
        jo = to_jl(o)
         
        #println(jo)
        
        #h1 = Histogram(jo[1], jo[2])
        #h2 = Histogram(jo[1], jo[3])

    #    println(h2)
        hd[symbol(kn)] = ErrorHistogram{Float64, length(jo[1]), typeof(jo[1])}(jo[1], jo[2], jo[3].^2, :right)
    #    println(jo)
        #println(kn, " ", clname, " ", o,)
    end
    Close(tf)
    return hd
end

quadsum(arr::AbstractArray) = arr.^2 |> sum |> sqrt
quadsum(arr::AbstractArray, n::Int64) = sum(arr.^2, n) |> sqrt

import Base.sum
function sum{T<:Real, N, E}(h::ErrorHistogram{T,N,E}, n::Int64)
    inds = [1:N]
    splice!(inds, n)
    newedges = deepcopy(h.edges[inds])
    dims = map(x->size(x)[1]-1, newedges)
    #println(newedges, " ", dims)
    return ErrorHistogram{T, N-1, typeof(newedges)}(
        newedges, reshape(sum(h.weights, n), dims...), reshape(sum(sqrt(h.weights_sq), n).^2, dims...), :right
    )
end

import Base.getindex
function getindex{T<:Real, N, E}(h::ErrorHistogram{T, N, E}, inds...)
    newinds = Any[]
    for i in inds
        push!(newinds, i.start:i.stop+1)
    end
    edgs = tuple([h.edges[i][newinds[i]] for i=1:length(newinds)]...)
    return ErrorHistogram{T, N, E}(edgs, h.weights[inds...], h.weights_sq[inds...], :right)
end

import Base.size
size(h::ErrorHistogram, n) = size(h.weights, n)
import Base.ndims
ndims(h::ErrorHistogram) = length(h.edges)

export to_jl, read_file

else
warn("Could not import ROOT.jl, ROOT functionality not defined.")
end

end # module
