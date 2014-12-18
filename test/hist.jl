using HEP
using Base.Test

h1 = ErrorHistogram(linspace(0, 1, 11))
@test ndims(h1) == 1
@test nbins(h1) == (10, )

h2 = ErrorHistogram(linspace(0, 1, 11), linspace(0, 1, 3))
@test ndims(h2) == 2
@test nbins(h2) == (10, 2)

push!(h1, 0.0)
push!(h1, 0.0)
push!(h1, 0.2, 0.5)
push!(h1, 0.2, 0.5)
push!(h1, 0.2, 0.5)

@test all(contents(h1) .== Float64[2.0, 0.0, 1.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
@test_approx_eq sum(errors(h1) .- Float64[sqrt(2.0), 0.0, sqrt(0.75), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]) 0.0

push!(h2, (0.0, 0.8), 0.5)
@test contents(h2)[1,2] == 0.5

htot = h2 + h2

@test_approx_eq integral(htot)[1] 2*integral(h2)[1]
@test_approx_eq integral(htot)[2] sqrt(2)*integral(h2)[2]