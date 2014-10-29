module HEP

#A 4-tuple representing a 4 vector in the natural system of units
immutable FourVector
    t::Float64
    x::Float64
    y::Float64
    z::Float64
end

+(v1::FourVector, v2::FourVector) = FourVector(v1.t+v2.t, v1.x+v2.x, v1.y+v2.y, v1.z+v2.z)
*(v1::FourVector, a::Real) = FourVector(a * v1.t, a * v1.x, a * v1.y, a*v1.z)

#magnitude
l2(v::FourVector) = v.t^2-v.x^2-v.y^2-v.z^2
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

#spherical components
function FourVectorSph(perp, eta, phi, l)
    perp = abs(perp)
    x = perp * cos(phi)
    y = perp * sin(phi)
    z = perp * sinh(eta)
    t = l>0 ? sqrt(l^2+x^2+y^2+z^2) : sqrt(-l^2+x^2+y^2+z^2)
    return FourVector(t, x, y, z)
end

#calculates the delta R between 2 four momenta
function deltar(v1::FourVector, v2::FourVector)
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
export deltar
export +
end # module
