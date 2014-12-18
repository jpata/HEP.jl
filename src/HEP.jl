module HEP

import Base.convert

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
    perp::Float64 #perpendicular component (e.g. in case of 4-momentum -> pt or transverse mometum)
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
end # module
