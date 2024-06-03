using PowerLaws
using DelimitedFiles
using Distributions
using Test

# tests are based moby dick word distribution
# and other data sets from http://tuvalu.santafe.edu/~aaronc/powerlaws/data.htm
# results are compared to result from R poweRlaw package
# test are done with tolerance to 1e-8
tolerance = 1e-8

pkg_path = paths -> joinpath(dirname(pathof(PowerLaws)), "..", paths...)
moby_data = vec(readdlm(pkg_path(["data", "moby_dick.txt"]), ' ', Int))
cities = vec(readdlm(pkg_path(["data", "cities.txt"]), ' ', Int))
electrical_blackout = vec(readdlm(pkg_path(["data", "electrical_blackouts_US.txt"]), ' ', Int))
population = vec(readdlm(pkg_path(["data", "population.txt"]), ' ', '\n'))

println("test: discrete powerlaw")
est = estimate_xmin(moby_data, dis_powerlaw)
@test isapprox(est[1].α, 1.9527275186009831, atol=tolerance) broken = true
@test est[1].θ == 7.0
@test isapprox(est[2], 0.008252950457938835,  atol=tolerance) broken = true

est = estimate_xmin(moby_data, dis_powerlaw, xmins = [2,3,4,10,20])
@test isapprox(est[1].α, 1.9550379794745847, atol=tolerance) broken = true
@test est[1].θ == 10.0
@test isapprox(est[2], 0.011867106465087929, atol=tolerance) broken = true

est1 = estimate_xmin(cities, dis_powerlaw)
@test isapprox(est1[1].α, 1.614392642791718, atol=tolerance) broken = true
@test est1[1].θ == 1021.0
@test isapprox(est1[2], 0.06088582981503132, atol=tolerance)

println("test: discrete powerlaw - bootstrap")
est1 = estimate_xmin(electrical_blackout, dis_powerlaw)
@test isapprox(est1[1].α, 1.2201523540878356, atol=tolerance)
@test est1[1].θ == 1000.0
@test isapprox(est1[2], 0.36278306128485405, atol=tolerance)

bootstr = bootstrap(moby_data,dis_powerlaw,no_of_sims = 15)
@test length(bootstr) == 15

bootstr = bootstrap(moby_data,est[1],no_of_sims = 15)
@test length(bootstr) == 15

# Unimplemented parallel tests:
# bootstr = bootstrap(moby_data,dis_powerlaw,3,no_of_sims = 15)
# @test length(bootstr) == 15

# bootstr = bootstrap(moby_data,est[1],2,no_of_sims = 15)
# @test length(bootstr) == 15

println("test: continuous powerlaw")
est = estimate_xmin(moby_data, con_powerlaw)
@test isapprox(est[1].α, 1.9335240463440206, atol=tolerance)
@test est[1].θ == 26.0
@test isapprox(est[2], 0.03204776085711486, atol=tolerance)

est = estimate_xmin(moby_data, con_powerlaw, xmins = [2,3,4,10,20])
@test isapprox(est[1].α, 1.9514282451732778, atol=tolerance)
@test est[1].θ == 20.0
@test isapprox(est[2], 0.04416094210009813, atol=tolerance)

est = estimate_xmin(collect(1:100),con_powerlaw)
@test isapprox(est[1].α, 9.178824791215408, atol=tolerance)
@test est[1].θ == 79.0
@test isapprox(est[2], 0.1824330509645657, atol=tolerance)

est1 = estimate_xmin(population,con_powerlaw)
@test isapprox(est1[1].α, 2.0337071519535046, atol=tolerance)
@test est1[1].θ == 96479.0
@test isapprox(est1[2], 0.0, atol=tolerance)

println("test: continuous powerlaw - bootstrap")
bootstr = bootstrap(moby_data,con_powerlaw,no_of_sims = 15)
@test length(bootstr) == 15

bootstr = bootstrap(moby_data,est[1],no_of_sims = 15)
@test length(bootstr) == 15

# Parallel tests:
# bootstr = bootstrap(moby_data,con_powerlaw,3,no_of_sims = 15)
# @test length(bootstr) == 15

# bootstr = bootstrap(moby_data,est[1],2,no_of_sims = 15)
# @test length(bootstr) == 15

println("test: compare distribution")
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
@test isapprox(cmpd.V_test_stat, 5.7463540130003254, atol=tolerance)
@test isapprox(cmpd.V_p_val, 0.9999999954405856, atol=tolerance)
@test cmpd.V_preff_distr == 1
@test cmpd.C_b == 62
@test isapprox(cmpd.C_p_val, 0.020978735677851468, atol=tolerance)
@test cmpd.C_preff_distr == 1

cmpd = compare_distributions(d1,Poisson,data)
@test typeof(cmpd) == compare_distributions
@test cmpd.data == data
@test cmpd.log_likehoods_ratio == ll_hoods_r
@test cmpd.xmin == 1
@test cmpd.sig_level == 0.05
@test isapprox(cmpd.V_test_stat, 5.7463540130003254, atol=tolerance)
@test isapprox(cmpd.V_p_val, 0.9999999954405856, atol=tolerance)
@test cmpd.V_preff_distr == 1
@test cmpd.C_b == 62
@test isapprox(cmpd.C_p_val, 0.020978735677851468, atol=tolerance)
@test cmpd.C_preff_distr == 1

cmpd = compare_distributions(dis_powerlaw,Poisson,data)
@test typeof(cmpd) == compare_distributions
@test cmpd.data == data
@test cmpd.log_likehoods_ratio == ll_hoods_r
@test cmpd.xmin == 1
@test cmpd.sig_level == 0.05
@test isapprox(cmpd.V_test_stat, 5.7463540130003254, atol=tolerance)
@test isapprox(cmpd.V_p_val, 0.9999999954405856, atol=tolerance)
@test cmpd.V_preff_distr == 1
@test cmpd.C_b == 62
@test isapprox(cmpd.C_p_val, 0.020978735677851468, atol=tolerance)
@test cmpd.C_preff_distr == 1

@test_throws ArgumentError compare_distributions(dis_powerlaw,con_powerlaw,data)
@test_throws MethodError compare_distributions(d1,con_powerlaw,data)
@test_throws MethodError compare_distributions(d1,con_powerlaw(),data)


moby_data = sort(moby_data)
d1 = estimate_xmin(moby_data,con_powerlaw)[1]
d2 = fit(Exponential,moby_data[18072:end])
cmpd = compare_distributions(d1,d2,moby_data,26)
@test cmpd.xmin == 26
@test cmpd.sig_level == 0.05
@test isapprox(cmpd.V_test_stat, 8.710089666437788, atol=tolerance)
@test isapprox(cmpd.V_p_val, 1.0, atol=tolerance)
@test cmpd.V_preff_distr == 1
@test cmpd.C_b == 566
@test isapprox(cmpd.C_p_val, 0.0, atol=tolerance)
@test cmpd.C_preff_distr == 1

d1 = estimate_xmin(moby_data,dis_powerlaw)[1]
d2 = fit(Poisson,moby_data[15898:end])
cmpd = compare_distributions(d1,d2,moby_data,7.0)
@test cmpd.xmin == 7.0
@test cmpd.sig_level == 0.05
@test isapprox(cmpd.V_test_stat, 4.448486373104684, atol=tolerance) broken = true
@test isapprox(cmpd.V_p_val, 0.9999956761232807, atol=tolerance)
@test cmpd.V_preff_distr == 1
@test cmpd.C_b == 2757
@test isapprox(cmpd.C_p_val, 0.0, atol=tolerance)
@test cmpd.C_preff_distr == 1
