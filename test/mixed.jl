@testitem "Impossible comparisons" begin
    using Distributions

    data = collect(1:100)
    d1 = fit(DiscretePowerLaw, data)

    @test_throws ArgumentError DistributionComparison(DiscretePowerLaw, ContinuousPowerLaw, data)
    @test_throws MethodError DistributionComparison(d1, ContinuousPowerLaw, data)
    @test_throws MethodError DistributionComparison(d1, ContinuousPowerLaw(), data)
end

@testitem "Vuong test compares against the critical value" begin
    using Distributions

    # This dataset yields a Vuong statistic of ≈1.25, which sits between the
    # buggy threshold cdf(Normal(), 0.975) ≈ 0.84 and the correct critical
    # value quantile(Normal(), 0.975) ≈ 1.96. Comparing against the cdf used
    # to (wrongly) declare distribution 1 preferred; against the critical
    # value the result is correctly inconclusive (0).
    data = [4, 1, 1, 1, 1, 2, 1, 1, 4, 1, 1, 1, 1, 1, 2, 1, 4, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 2, 1, 2, 1, 1, 1, 1, 2, 2, 1, 1, 2, 2, 67, 4, 2, 1, 1,
        1, 1, 5, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 2]
    d1 = DiscretePowerLaw(2.0, 1.0)
    cmpd = DistributionComparison(d1, Poisson, data)

    @test cmpd.V_test_stat ≈ 1.2470810909269956
    @test cdf(Normal(), 0.975) < cmpd.V_test_stat < quantile(Normal(), 0.975)
    @test cmpd.V_preff_distr == 0
end

@testitem "estimate_parameters errors when no xmin is valid" begin
    @test_throws ArgumentError estimate_parameters(collect(1:10), ContinuousPowerLaw, xmins=[9999])
    @test_throws ArgumentError estimate_parameters(collect(1:10), DiscretePowerLaw, xmins=[9999])
end