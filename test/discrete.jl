# Numeric values have been obtained through the use of the `powerlaw` package in Python.

pkg_path = paths -> joinpath(dirname(pathof(PowerLaws)), "..", paths...)

cities = vec(readdlm(pkg_path(["data", "cities.txt"]), ' ', Int))
electrical_blackout = vec(readdlm(pkg_path(["data", "electrical_blackouts_US.txt"]), ' ', Int))
population = vec(readdlm(pkg_path(["data", "population.txt"]), ' ', '\n'))

@testitem "Discrete x_min estimation" begin
    using DelimitedFiles
    using Distributions

    data_dir = joinpath(dirname(pathof(PowerLaws)), "..", "data")
    moby_data = vec(readdlm(joinpath(data_dir, "moby_dick.txt"), Int))
    cities = vec(readdlm(joinpath(data_dir, "cities.txt"), Int))

    est = estimate_parameters(moby_data, DiscretePowerLaw)
    @test est[1].α ≈ 1.95015723
    @test est[1].θ == 7.0
    @test est[2] ≈ 0.00922886388

    dist = fit(DiscretePowerLaw, moby_data[moby_data .>= est[1].θ])
    @test est[1].α ≈ shape(dist)
    @test est[1].θ ≈ scale(dist)

    est1 = estimate_parameters(cities, DiscretePowerLaw)
    @test est1[1].α ≈ 1.61439261
    @test est1[1].θ == 1021.0
    @test est1[2] ≈ 0.0608858393
end

@testitem "Discrete x_min from options" begin
    using DelimitedFiles
    using Distributions

    data_dir = joinpath(dirname(pathof(PowerLaws)), "..", "data")
    moby_data = vec(readdlm(joinpath(data_dir, "moby_dick.txt"), Int))

    est = estimate_parameters(moby_data, DiscretePowerLaw, xmins=[2, 3, 4, 10, 20])
    @test est[1].α ≈ 1.95381938
    @test est[1].θ == 10.0
    @test est[2] ≈ 0.0122405536

    est = estimate_parameters(moby_data, DiscretePowerLaw, xmins=[1,])
    dist = fit(DiscretePowerLaw, moby_data)
    @test est[1].α ≈ shape(dist)
    @test est[1].θ ≈ scale(dist)
end

@testitem "Discrete bootstrap" begin
    using DelimitedFiles

    data_dir = joinpath(dirname(pathof(PowerLaws)), "..", "data")
    moby_data = vec(readdlm(joinpath(data_dir, "moby_dick.txt"), Int))
    electrical_blackout = vec(readdlm(joinpath(data_dir, "electrical_blackouts_US.txt"), Int))

    est = estimate_parameters(electrical_blackout, DiscretePowerLaw)
    @test est[1].α ≈ 1.22015235
    @test est[1].θ == 1000.0
    @test est[2] ≈ 0.362783061

    bootstr = bootstrap(moby_data, DiscretePowerLaw, no_of_sims=15)
    @test length(bootstr) == 15

    bootstr = bootstrap(electrical_blackout, est[1], no_of_sims=12)
    @test length(bootstr) == 12
end

@testitem "Compare discrete fitted-fitted" begin
    using Distributions

    data = collect(1:100)
    d1 = fit(DiscretePowerLaw, data)
    d2 = fit(Poisson, data)
    f_lpdf = (distribution, data) -> map(Base.Fix1(logpdf, distribution), data)
    ll_hoods_r = f_lpdf(d1, data) - f_lpdf(d2, data)
    cmpd = DistributionComparison(d1, d2, data)
    @test typeof(cmpd) == DistributionComparison
    @test cmpd.data == data
    @test cmpd.log_likelihood_ratio == ll_hoods_r
    @test cmpd.xmin == 1
    @test cmpd.sig_level == 0.05
    @test cmpd.V_test_stat ≈ 5.74635401
    @test cmpd.V_p_val ≈ 0.999999995
    @test cmpd.V_preff_distr == 1
    @test cmpd.C_b == 62
    @test cmpd.C_p_val ≈ 0.0209787356
    @test cmpd.C_preff_distr == 1
end

@testitem "Compare discrete fitted-type" begin
    using Distributions

    data = collect(1:100)
    d1 = fit(DiscretePowerLaw, data)
    d2 = fit(Poisson, data)
    f_lpdf = (distribution, data) -> map(Base.Fix1(logpdf, distribution), data)
    ll_hoods_r = f_lpdf(d1, data) - f_lpdf(d2, data)
    cmpd = DistributionComparison(d1, Poisson, data)
    @test typeof(cmpd) == DistributionComparison
    @test cmpd.data == data
    @test cmpd.log_likelihood_ratio == ll_hoods_r
    @test cmpd.xmin == 1
    @test cmpd.sig_level == 0.05
    @test cmpd.V_test_stat ≈ 5.74635401
    @test cmpd.V_p_val ≈ 0.999999995
    @test cmpd.V_preff_distr == 1
    @test cmpd.C_b == 62
    @test cmpd.C_p_val ≈ 0.0209787356
    @test cmpd.C_preff_distr == 1
end

@testitem "Compare discrete type-type" begin
    using Distributions

    data = collect(1:100)
    d1 = fit(DiscretePowerLaw, data)
    d2 = fit(Poisson, data)
    f_lpdf = (distribution, data) -> map(Base.Fix1(logpdf, distribution), data)
    ll_hoods_r = f_lpdf(d1, data) - f_lpdf(d2, data)
    cmpd = DistributionComparison(DiscretePowerLaw, Poisson, data)
    @test typeof(cmpd) == DistributionComparison
    @test cmpd.data == data
    @test cmpd.log_likelihood_ratio == ll_hoods_r
    @test cmpd.xmin == 1
    @test cmpd.sig_level == 0.05
    @test cmpd.V_test_stat ≈ 5.74635401
    @test cmpd.V_p_val ≈ 0.999999995
    @test cmpd.V_preff_distr == 1
    @test cmpd.C_b == 62
    @test cmpd.C_p_val ≈ 0.0209787356
    @test cmpd.C_preff_distr == 1
end

@testitem "Compare discrete estimated-fitted" begin
    using DelimitedFiles
    using Distributions

    data_dir = joinpath(dirname(pathof(PowerLaws)), "..", "data")
    moby_data = sort(vec(readdlm(joinpath(data_dir, "moby_dick.txt"), Int)))

    d1 = estimate_parameters(moby_data, DiscretePowerLaw)[1]
    d2 = fit(Poisson, moby_data[15898:end])
    cmpd = DistributionComparison(d1, d2, moby_data, 7.0)
    @test cmpd.xmin == 7.0
    @test cmpd.sig_level == 0.05
    @test cmpd.V_test_stat ≈ 4.448483568
    @test cmpd.V_p_val ≈ 0.999995676
    @test cmpd.V_preff_distr == 1
    @test cmpd.C_b == 2757
    @test cmpd.C_p_val ≈ 0.0
    @test cmpd.C_preff_distr == 1
end

@testitem "Discrete distribution functions" begin
    using Distributions
    using SpecialFunctions: zeta

    α, θ = 2.5, 2.0
    d = DiscretePowerLaw(α, θ)
    z = zeta(α, θ)

    @test params(d) == (α, θ)
    @test shape(d) == α
    @test scale(d) == θ
    @test params(DiscretePowerLaw(3.0)) == (3.0, 1.0)
    @test params(DiscretePowerLaw()) == (1.0, 1.0)

    # pdf / logpdf compared with the closed form (x below θ is outside the support)
    for x in (1.0, 2.0, 5.0, 10.0)
        expected = x < θ ? 0.0 : x^(-α) / z
        @test pdf(d, x) ≈ expected
        @test logpdf(d, x) ≈ (x < θ ? -Inf : -log(z) - α * log(x))
    end

    # ccdf / cdf / logccdf / logcdf via the Hurwitz zeta function
    for x in (2.0, 5.0, 10.0)
        @test ccdf(d, x) ≈ zeta(α, x) / z
        @test cdf(d, x) ≈ 1.0 - zeta(α, x) / z
        @test logccdf(d, x) ≈ log(zeta(α, x) / z)
        @test logcdf(d, x) ≈ log(1.0 - zeta(α, x) / z)
    end

    # array methods agree with the scalar ones
    xs = [1.0, 3.0, 7.0]
    @test pdf(d, xs) == [pdf(d, x) for x in xs]
    @test logpdf(d, xs) == [logpdf(d, x) for x in xs]
    @test ccdf(d, xs) ≈ [ccdf(d, x) for x in xs]
    @test cdf(d, xs) ≈ [cdf(d, x) for x in xs]
    @test logccdf(d, xs) ≈ [logccdf(d, x) for x in xs]
    @test logcdf(d, xs) ≈ [logcdf(d, x) for x in xs]
end

@testitem "Discrete sampling returns integers" begin
    using Distributions, Random

    d = DiscretePowerLaw(2.5, 3.0)
    @test rand(d) isa Integer
    @test rand(d) >= 3

    s = rand(d, 500)
    @test eltype(s) <: Integer
    @test length(s) == 500
    @test all(>=(3), s)

    @test rand(MersenneTwister(7), d, 10) == rand(MersenneTwister(7), d, 10)
end

@testitem "Discrete constructor validation" begin
    @test_throws ArgumentError DiscretePowerLaw(-1.0, 1.0)
    @test_throws ArgumentError DiscretePowerLaw(0.0, 1.0)
    @test_throws ArgumentError DiscretePowerLaw(1.0, -1.0)
    @test_throws ArgumentError DiscretePowerLaw(1.0, 0.0)
end

@testitem "Discrete estimate_parameters validates input" begin
    # Non-integer data is rejected
    @test_throws ArgumentError estimate_parameters([1.5, 2.5, 3.5], DiscretePowerLaw)
    # No supplied xmin matches the data -> clean ArgumentError, not a BoundsError
    @test_throws ArgumentError estimate_parameters(collect(1:10), DiscretePowerLaw, xmins=[9999])
end

@testitem "Discrete bootstrap_p" begin
    data = [1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 5, 6, 7, 8, 9,
        10, 12, 15, 20, 30, 40, 55, 70, 100]
    stats, p = bootstrap_p(data, DiscretePowerLaw, no_of_sims=5, seed=1)
    @test length(stats) == 5
    @test all(s -> s[1] isa DiscretePowerLaw, stats)
    @test 0.0 <= p <= 1.0
end
