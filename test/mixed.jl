@testitem "Impossible comparisons" begin
    using Distributions

    data = collect(1:100)
    d1 = fit(dis_powerlaw,data)

    @test_throws ArgumentError compare_distributions(dis_powerlaw,con_powerlaw,data)
    @test_throws MethodError compare_distributions(d1,con_powerlaw,data)
    @test_throws MethodError compare_distributions(d1,con_powerlaw(),data)
end