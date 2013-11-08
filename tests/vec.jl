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

v = HEP.FourVectorSph(10.0, -2.0, 3.14/2.0, 100.0)
@test_approx_eq perp(v) 10.0
@test_approx_eq eta(v) -2.0000000000000004
@test_approx_eq phi(v) 1.57
@test_approx_eq rho(v) 37.62195691083632
@test_approx_eq theta(v) 2.872556662840912
@test_approx_eq perp(v) 10.0
@test_approx_eq l(v) 99.99999999999999