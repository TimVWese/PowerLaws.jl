pkg_path = paths -> joinpath(dirname(pathof(PowerLaws)), "..", paths...)

cities = vec(readdlm(pkg_path(["data", "cities.txt"]), ' ', Int))
electrical_blackout = vec(readdlm(pkg_path(["data", "electrical_blackouts_US.txt"]), ' ', Int))
population = vec(readdlm(pkg_path(["data", "population.txt"]), ' ', '\n'))

@testitem "Discrete x_min estimation" begin
    using DelimitedFiles

    data_dir = joinpath(dirname(pathof(PowerLaws)), "..", "data")
    moby_data = vec(readdlm(joinpath(data_dir, "moby_dick.txt"), Int))
    cities = vec(readdlm(joinpath(data_dir, "cities.txt"), Int))

    est = estimate_xmin(moby_data, dis_powerlaw)
    @test est[1].α ≈ 1.95272751 broken = true
    @test est[1].θ == 7.0
    @test est[2] ≈ 0.00825295045 broken = true

    est1 = estimate_xmin(cities, dis_powerlaw)
    @test est1[1].α ≈ 1.61439264 broken = true
    @test est1[1].θ == 1021.0
    @test est1[2] ≈ 0.0608858298 broken = true
end

@testitem "Discrete x_min from options" begin
    using DelimitedFiles

    data_dir = joinpath(dirname(pathof(PowerLaws)), "..", "data")
    moby_data = vec(readdlm(joinpath(data_dir, "moby_dick.txt"), Int))

    est = estimate_xmin(moby_data, dis_powerlaw, xmins = [2,3,4,10,20])
    @test isapprox(est[1].α, 1.9550379794745847, atol=tolerance) broken = true
    @test est[1].θ == 10.0
    @test isapprox(est[2], 0.011867106465087929, atol=tolerance) broken = true
end

@testitem "Discrete bootstrap" begin
    using DelimitedFiles

    data_dir = joinpath(dirname(pathof(PowerLaws)), "..", "data")
    moby_data = vec(readdlm(joinpath(data_dir, "moby_dick.txt"), Int))
    electrical_blackout = vec(readdlm(joinpath(data_dir, "electrical_blackouts_US.txt"), Int))

    est = estimate_xmin(electrical_blackout, dis_powerlaw)
    @test est[1].α ≈ 1.22015235
    @test est[1].θ == 1000.0
    @test est[2] ≈ 0.362783061

    bootstr = bootstrap(moby_data,dis_powerlaw,no_of_sims = 15)
    @test length(bootstr) == 15

    bootstr = bootstrap(electrical_blackout,est[1],no_of_sims = 12)
    @test length(bootstr) == 12
end

@testitem "Compare discrete fitted-fitted" begin
    using Distributions

    data = collect(1:100)
    d1 = fit(dis_powerlaw,data)
    d2 = fit(Poisson,data)
    f_lpdf = (distribution, data) -> map(Base.Fix1(logpdf, distribution), data)
    ll_hoods_r = f_lpdf(d1,data) - f_lpdf(d2,data)
    cmpd = compare_distributions(d1,d2,data)
    @test typeof(cmpd) == compare_distributions
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
    d1 = fit(dis_powerlaw,data)
    d2 = fit(Poisson,data)
    f_lpdf = (distribution, data) -> map(Base.Fix1(logpdf, distribution), data)
    ll_hoods_r = f_lpdf(d1,data) - f_lpdf(d2,data)
    cmpd = compare_distributions(d1,Poisson,data)
    @test typeof(cmpd) == compare_distributions
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
    d1 = fit(dis_powerlaw,data)
    d2 = fit(Poisson,data)
    f_lpdf = (distribution, data) -> map(Base.Fix1(logpdf, distribution), data)
    ll_hoods_r = f_lpdf(d1,data) - f_lpdf(d2,data)
    cmpd = compare_distributions(dis_powerlaw,Poisson,data)
    @test typeof(cmpd) == compare_distributions
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