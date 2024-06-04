@testitem "Impossible comparisons" begin
    using Distributions

    data = collect(1:100)
    d1 = fit(DiscretePowerLaw,data)

    @test_throws ArgumentError DistributionComparison(DiscretePowerLaw, ContinuousPowerLaw, data)
    @test_throws MethodError DistributionComparison(d1, ContinuousPowerLaw, data)
    @test_throws MethodError DistributionComparison(d1, ContinuousPowerLaw(), data)
end