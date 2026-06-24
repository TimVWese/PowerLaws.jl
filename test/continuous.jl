@testitem "Continuous x_min estimation" begin
    using DelimitedFiles, Distributions

    data_dir = joinpath(dirname(pathof(PowerLaws)), "..", "data")
    moby_data = vec(readdlm(joinpath(data_dir, "moby_dick.txt"), Int))
    population = vec(readdlm(joinpath(data_dir, "population.txt"), ' ', '\n'))

    est = estimate_parameters(moby_data, ContinuousPowerLaw)
    @test est[1].α ≈ 1.93352404
    @test est[1].θ == 26.0
    @test est[2] ≈ 0.0320477608

    dist = fit(ContinuousPowerLaw, moby_data[moby_data .>= est[1].θ])
    @test est[1].α ≈ shape(dist)
    @test est[1].θ ≈ scale(dist)

    est = estimate_parameters(collect(1:100), ContinuousPowerLaw)
    @test est[1].α ≈ 9.17882479
    @test est[1].θ == 79.0
    @test est[2] ≈ 0.182433050

    est1 = estimate_parameters(population, ContinuousPowerLaw)
    @test est1[1].α ≈ 2.03370715
    @test est1[1].θ == 96479.0
    @test est1[2] ≈ 0.0
end

@testitem "Continuous x_min given options" begin
    using DelimitedFiles
    using Distributions

    data_dir = joinpath(dirname(pathof(PowerLaws)), "..", "data")
    moby_data = vec(readdlm(joinpath(data_dir, "moby_dick.txt"), Int))

    est = estimate_parameters(moby_data, ContinuousPowerLaw, xmins=[2, 3, 4, 10, 20])
    @test est[1].α ≈ 1.951428245
    @test est[1].θ == 20.0
    @test est[2] ≈ 0.0441609421

    est = estimate_parameters(moby_data, ContinuousPowerLaw, xmins=[1,])
    dist = fit(ContinuousPowerLaw, moby_data)
    @test est[1].α ≈ shape(dist)
    @test est[1].θ ≈ scale(dist)
end

@testitem "Continuous bootstrap" begin
    using DelimitedFiles

    data_dir = joinpath(dirname(pathof(PowerLaws)), "..", "data")
    moby_data = vec(readdlm(joinpath(data_dir, "moby_dick.txt"), Int))

    bootstr = bootstrap(moby_data, ContinuousPowerLaw, no_of_sims=15)
    @test length(bootstr) == 15

    est = estimate_parameters(moby_data, ContinuousPowerLaw)
    bootstr = bootstrap(moby_data, est[1], no_of_sims=12)
    @test length(bootstr) == 12
end

@testitem "Compare continuous fitted-fitted" begin
    using DelimitedFiles
    using Distributions

    data_dir = joinpath(dirname(pathof(PowerLaws)), "..", "data")
    moby_data = sort(vec(readdlm(joinpath(data_dir, "moby_dick.txt"), Int)))

    d1 = estimate_parameters(moby_data, ContinuousPowerLaw)[1]
    d2 = fit(Exponential, moby_data[18072:end])
    cmpd = DistributionComparison(d1, d2, moby_data, 26)
    @test cmpd.xmin == 26
    @test cmpd.sig_level == 0.05
    @test cmpd.V_test_stat ≈ 8.71008966
    @test cmpd.V_p_val ≈ 1.0
    @test cmpd.V_preff_distr == 1
    @test cmpd.C_b == 566
    @test cmpd.C_p_val ≈ 0.0
    @test cmpd.C_preff_distr == 1
end

@testitem "Continuous distribution functions" begin
    using Distributions

    α, θ = 2.5, 2.0
    d = ContinuousPowerLaw(α, θ)

    @test params(d) == (α, θ)
    @test shape(d) == α
    @test scale(d) == θ
    @test mode(d) == θ
    @test params(ContinuousPowerLaw(3.0)) == (3.0, 1.0)
    @test params(ContinuousPowerLaw()) == (1.0, 1.0)

    # pdf / logpdf compared with the closed form (x below θ is outside the support)
    for x in (1.0, 2.0, 5.0, 10.0)
        expected = x < θ ? 0.0 : ((α - 1.0) / θ) * (x / θ)^(-α)
        @test pdf(d, x) ≈ expected
        @test logpdf(d, x) ≈ (x < θ ? -Inf : log(expected))
    end

    # ccdf / cdf / logccdf / logcdf
    for x in (2.0, 5.0, 10.0)
        @test ccdf(d, x) ≈ (x / θ)^(1.0 - α)
        @test cdf(d, x) ≈ 1.0 - (x / θ)^(1.0 - α)
        @test logccdf(d, x) ≈ log((x / θ)^(1.0 - α))
        @test logcdf(d, x) ≈ log(1.0 - (x / θ)^(1.0 - α))
    end

    # quantile / cquantile invert cdf / ccdf
    for p in (0.1, 0.5, 0.9)
        @test cdf(d, quantile(d, p)) ≈ p
        @test ccdf(d, cquantile(d, p)) ≈ p
    end

    # moments (evaluated where they are finite, NaN/Inf otherwise)
    @test mean(ContinuousPowerLaw(2.5, 2.0)) ≈ 2.0 * (1.5 / 0.5)
    @test mean(ContinuousPowerLaw(1.5, 2.0)) == Inf
    @test var(ContinuousPowerLaw(4.0, 2.0)) ≈ (2.0^2 * 3.0) / (2.0^2 * 1.0)
    @test var(ContinuousPowerLaw(2.5, 1.0)) == Inf
    @test median(ContinuousPowerLaw(2.5, 2.0)) ≈ 2.0 * 2.0^(1.0 / 1.5)
    @test isnan(median(ContinuousPowerLaw(0.5, 1.0)))
    @test skewness(ContinuousPowerLaw(5.0, 1.0)) ≈ (2.0 * 5.0 / 1.0) * sqrt(2.0 / 4.0)
    @test isnan(skewness(ContinuousPowerLaw(4.0, 1.0)))
    @test kurtosis(ContinuousPowerLaw(6.0, 1.0)) ≈
          (6.0 * (5.0^3 + 5.0^2 - 6.0 * 5.0 - 2.0)) / (5.0 * 2.0 * 1.0)
    @test isnan(kurtosis(ContinuousPowerLaw(5.0, 1.0)))
    @test entropy(d) ≈ log(θ / (α - 1.0)) + 1.0 / (α - 1.0) + 1.0
end

@testitem "Continuous pdf/logpdf on arrays match scalars" begin
    using Distributions

    # Regression: the array `pdf` method used to evaluate `x` (the whole array)
    # instead of the loop variable, so it never matched the scalar method.
    d = ContinuousPowerLaw(2.5, 2.0)
    xs = [0.5, 2.0, 3.0, 10.0]
    @test pdf(d, xs) == [pdf(d, x) for x in xs]
    @test logpdf(d, xs) == [logpdf(d, x) for x in xs]
end

@testitem "Continuous sampling" begin
    using Distributions, Random

    d = ContinuousPowerLaw(2.5, 3.0)
    @test rand(d) isa Float64
    @test rand(d) >= 3.0

    s = rand(d, 500)
    @test s isa Vector{Float64}
    @test length(s) == 500
    @test all(>=(3.0), s)

    # rand integrates with an explicit RNG, so it is reproducible
    @test rand(MersenneTwister(42), d, 10) == rand(MersenneTwister(42), d, 10)
end

@testitem "Continuous constructor validation" begin
    @test_throws ArgumentError ContinuousPowerLaw(-1.0, 1.0)
    @test_throws ArgumentError ContinuousPowerLaw(0.0, 1.0)
    @test_throws ArgumentError ContinuousPowerLaw(1.0, -1.0)
    @test_throws ArgumentError ContinuousPowerLaw(1.0, 0.0)
end

@testitem "Continuous bootstrap_p" begin
    data = Float64.([1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 5, 6, 7, 8, 9,
        10, 12, 15, 20, 30, 40, 55, 70, 100])
    stats, p = bootstrap_p(data, ContinuousPowerLaw, no_of_sims=5, seed=1)
    @test length(stats) == 5
    @test all(s -> s[1] isa ContinuousPowerLaw, stats)
    @test 0.0 <= p <= 1.0
end
