using HEP
using Base.Test

v = HEP.FourVector(4.0, 1.0, 2.0, 3.0)
@test_approx_eq perp(v) 2.23606797749979
@test_approx_eq eta(v) 1.1035868415601455
@test_approx_eq phi(v) 1.1071487177940904
@test_approx_eq rho(v) 3.7416573867739413
@test_approx_eq theta(v) 0.6405223126794246
@test_approx_eq perp(v) 2.23606797749979
@test_approx_eq l(v) 1.4142135623730951

v2 = HEP.FourVectorSph(10.0, -2.0, 3.14/2.0, 100.0)
@test_approx_eq perp(v2) 10.0
@test_approx_eq eta(v2) -2.0000000000000004
@test_approx_eq phi(v2) 1.57
@test_approx_eq rho(v2) 37.62195691083632
@test_approx_eq theta(v2) 2.872556662840912
@test_approx_eq perp(v2) 10.0
@test_approx_eq l(v2) 99.99999999999999

@test_approx_eq deltar(v, v2) deltar(v2, v)
@test_approx_eq deltar(v, v) 0
@test_approx_eq deltar(v, v2) 3.137910545656924

v2_c = convert(v2, FourVector)

@test_approx_eq l(v2_c) l(v2)
@test_approx_eq x(v2_c) x(v2)
@test_approx_eq y(v2_c) y(v2)
@test_approx_eq z(v2_c) z(v2)

vt = v + v2_c

@test_approx_eq x(vt) 1.0079632671073326
@test_approx_eq y(vt) 11.999996829318347
@test_approx_eq z(vt) -33.26860407847019
@test_approx_eq t(vt) 110.84292976983