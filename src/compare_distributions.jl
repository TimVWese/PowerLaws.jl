
"""
    DistributionComparison

# Fields

* data: array which according which should be distributions compared
* log_likelihood_ratio: log likelihood ratio of data
* sig_level: sigma level
* xmin: smallest element which was used for comparing distributions
* V_test_stat: Vuong test statistic
* V_p_val: p-value from Vuong test
* V_preff_distr: preffered distibution according to Vuong test
* C_b: Clarke test b value - sum of possitive values in log likehood ratio
* C_p_val: p-value from Clarke test
* C_preff_distr: preffered distibution according to Clarke test
"""
struct DistributionComparison
    data::AbstractArray
    log_likelihood_ratio::AbstractArray
    sig_level::Float64
    xmin::Number
    V_test_stat::Float64
    V_p_val::Float64
    V_preff_distr::Int64
    C_b::Int64
    C_p_val::Float64
    C_preff_distr::Int64
end

function show(io::IO, x::DistributionComparison)
    println(io, "Data:  $(x.data)")
    println(io, "Log likehood ratio: $(x.log_likelihood_ratio)")
    println(io, "Significance level: $(x.sig_level)")
    println(io, "Xmin: $(x.xmin)")
    println(io, "Vuong - test statictic: $(x.V_test_stat)")
    println(io, "Vuong - p-value: $(x.V_p_val)")
    println(io, "Vuong - preffered distribution: $(x.V_preff_distr)")
    println(io, "Clarke test - number of possitive values in log likehood ratio: $(x.C_b)")
    println(io, "Clarke test - p-value: $(x.C_p_val)")
    println(io, "Clarke test - preffered distribution: $(x.C_preff_distr)")
end


"""
    DistributionComparison

This function calculate Vuong test and Clarke test for non nested distributions.
This is necessary since it is possible to fit power law distribution to any data set.
Function was implemented according to this [Non nested model selection for spatial count regression models with application to health insurance](https://mediatum.ub.tum.de/doc/1083601/1083601.pdf).

# Arguments
- `d1`: First distribution to be compared.
- `d2`: Second distribution to be compared.
- `data`: Data to be compared.
- `sig_level`: Significance level. Default is `0.05`.

# Returns
- `DistributionComparison`: Struct containing all necessary information about comparison.
"""
function DistributionComparison(d1::ContinuousPowerLaw, d2::Type{<:ContinuousUnivariateDistribution}, data::AbstractArray; sig_level = 0.05)
    xmin = d1.θ
    data = sort(data)
    data = data[findfirst(x -> x >= xmin, data): end]
    d2 = fit(d2, data)

    _compare_distributions(d1,d2,data,xmin,sig_level)
end

function DistributionComparison(d1::DiscretePowerLaw, d2::Type{<:DiscreteUnivariateDistribution}, data::AbstractArray; sig_level = 0.05)
    xmin = d1.θ
    data = sort(data)
    data = data[findfirst(x -> x >= xmin, data): end]
    d2 = fit(d2, data)

    _compare_distributions(d1,d2,data,xmin,sig_level)
end

function DistributionComparison(d1::Type{<:UnivariateDistribution}, d2::Type{<:UnivariateDistribution}, data::AbstractArray, xmin::Number = 0; sig_level = 0.05)
    if !(((d1 <: ContinuousUnivariateDistribution) && (d2 <: ContinuousUnivariateDistribution)) ||
            ((d1 <: DiscreteUnivariateDistribution) && (d2 <: DiscreteUnivariateDistribution)))
        throw(ArgumentError("Both distributions should be either continuous or discrete."))
    end

    if (xmin == 0)
        if (d1 <: DiscreteUnivariateDistribution) xmin = 1 end
        xmin = minimum(data)
    else
        data = sort(data)
        data = data[findfirst(x -> x >= xmin, data): end]
    end
    d1 = fit(d1, data)
    d2 = fit(d2, data)

    _compare_distributions(d1,d2,data,xmin,sig_level)
end

function DistributionComparison(d1::DiscreteUnivariateDistribution, d2::DiscreteUnivariateDistribution, data::AbstractArray, xmin::Number = 1; sig_level = 0.05)
   _compare_distributions(d1,d2,data,xmin,sig_level) 
end

function DistributionComparison(d1::ContinuousUnivariateDistribution, d2::ContinuousUnivariateDistribution, data::AbstractArray, xmin::Number = 0; sig_level = 0.05)
   _compare_distributions(d1,d2,data,xmin,sig_level) 
end

function _compare_distributions(d1::UnivariateDistribution, d2::UnivariateDistribution, data::AbstractArray, xmin::Number, sig_level)
    data = sort(data)
    xmin == 0 ? xmin = minimum(data) : data = data[findfirst(x -> x >= xmin, data): end]
    # Vuong's test
    log_likelihood_ratio = map(Base.Fix1(logpdf, d1), data) - map(Base.Fix1(logpdf, d2), data)
    n = length(log_likelihood_ratio)
    m = mean(log_likelihood_ratio)
    standard_deviation = std(log_likelihood_ratio)
    test_stat = sqrt(n) * m / standard_deviation
    v_p_val = cdf(Normal(), test_stat)
    v_preff_distr = 0
    if (test_stat > cdf(Normal(), 1 - sig_level/2))
        v_preff_distr = 1
    elseif (test_stat < -cdf(Normal(), 1 - sig_level/2))
        v_preff_distr = 2
    end
    
    # Clarke's test
    # b = number of positive values in log_likelihood_ratio
    b = sum(x->x > 0 ? 1 : 0, log_likelihood_ratio)
    preff_distr = 0
    if (b >= n/2)
        pval = 2 * (1 - cdf(Binomial(n,0.5),b-1))
        if (pval <= sig_level) preff_distr = 1 end
    else
        pval = 2 * (cdf(Binomial(n,0.5),b))
        if (pval <= sig_level) preff_distr = 2 end
    end

    return DistributionComparison(data,log_likelihood_ratio,sig_level,xmin,test_stat,v_p_val,v_preff_distr,b,pval,preff_distr)
end