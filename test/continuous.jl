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
