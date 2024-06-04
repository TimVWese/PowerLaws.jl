# Numeric values have been obtained through the use of the `powerlaw` package in Python.

pkg_path = paths -> joinpath(dirname(pathof(PowerLaws)), "..", paths...)

cities = vec(readdlm(pkg_path(["data", "cities.txt"]), ' ', Int))
electrical_blackout = vec(readdlm(pkg_path(["data", "electrical_blackouts_US.txt"]), ' ', Int))
population = vec(readdlm(pkg_path(["data", "population.txt"]), ' ', '\n'))

@testitem "Discrete x_min estimation" begin
    using DelimitedFiles

    data_dir = joinpath(dirname(pathof(PowerLaws)), "..", "data")
    moby_data = vec(readdlm(joinpath(data_dir, "moby_dick.txt"), Int))
    cities = vec(readdlm(joinpath(data_dir, "cities.txt"), Int))

    est = estimate_parameters(moby_data, DiscretePowerLaw)
    @test est[1].α ≈ 1.95015723
    @test est[1].θ == 7.0
    @test est[2] ≈ 0.00922886388

    est1 = estimate_parameters(cities, DiscretePowerLaw)
    @test est1[1].α ≈ 1.61439261
    @test est1[1].θ == 1021.0
    @test est1[2] ≈ 0.0608858393
end

@testitem "Discrete x_min from options" begin
    using DelimitedFiles

    data_dir = joinpath(dirname(pathof(PowerLaws)), "..", "data")
    moby_data = vec(readdlm(joinpath(data_dir, "moby_dick.txt"), Int))

    est = estimate_parameters(moby_data, DiscretePowerLaw, xmins = [2,3,4,10,20])
    @test est[1].α ≈ 1.95381938
    @test est[1].θ == 10.0
    @test est[2] ≈ 0.0122405536
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

    bootstr = bootstrap(moby_data,DiscretePowerLaw,no_of_sims = 15)
    @test length(bootstr) == 15

    bootstr = bootstrap(electrical_blackout,est[1],no_of_sims = 12)
    @test length(bootstr) == 12
end

@testitem "Compare discrete fitted-fitted" begin
    using Distributions

    data = collect(1:100)
    d1 = fit(DiscretePowerLaw,data)
    d2 = fit(Poisson,data)
    f_lpdf = (distribution, data) -> map(Base.Fix1(logpdf, distribution), data)
    ll_hoods_r = f_lpdf(d1,data) - f_lpdf(d2,data)
    cmpd = DistributionComparison(d1,d2,data)
    @test typeof(cmpd) == DistributionComparison
    @test cmpd.data == data
    @test cmpd.log_likehoods_ratio == ll_hoods_r
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
    d1 = fit(DiscretePowerLaw,data)
    d2 = fit(Poisson,data)
    f_lpdf = (distribution, data) -> map(Base.Fix1(logpdf, distribution), data)
    ll_hoods_r = f_lpdf(d1,data) - f_lpdf(d2,data)
    cmpd = DistributionComparison(d1,Poisson,data)
    @test typeof(cmpd) == DistributionComparison
    @test cmpd.data == data
    @test cmpd.log_likehoods_ratio == ll_hoods_r
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
    d1 = fit(DiscretePowerLaw,data)
    d2 = fit(Poisson,data)
    f_lpdf = (distribution, data) -> map(Base.Fix1(logpdf, distribution), data)
    ll_hoods_r = f_lpdf(d1,data) - f_lpdf(d2,data)
    cmpd = DistributionComparison(DiscretePowerLaw,Poisson,data)
    @test typeof(cmpd) == DistributionComparison
    @test cmpd.data == data
    @test cmpd.log_likehoods_ratio == ll_hoods_r
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
    d2 = fit(Poisson,moby_data[15898:end])
    cmpd = DistributionComparison(d1,d2,moby_data,7.0)
    @test cmpd.xmin == 7.0
    @test cmpd.sig_level == 0.05
    @test cmpd.V_test_stat ≈ 4.448483568
    @test cmpd.V_p_val ≈ 0.999995676
    @test cmpd.V_preff_distr == 1
    @test cmpd.C_b == 2757
    @test cmpd.C_p_val ≈ 0.0
    @test cmpd.C_preff_distr == 1
end
