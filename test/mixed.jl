@testitem "Impossible comparisons" begin
    using Distributions

    data = collect(1:100)
    d1 = fit(DiscretePowerLawDistribution,data)

    @test_throws ArgumentError DistributionComparison(DiscretePowerLawDistribution, ContinuousPowerLawDistribution, data)
    @test_throws MethodError DistributionComparison(d1, ContinuousPowerLawDistribution, data)
    @test_throws MethodError DistributionComparison(d1, ContinuousPowerLawDistribution(), data)
end