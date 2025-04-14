var documenterSearchIndex = {"docs":
[{"location":"reference/#PowerLaws.jl","page":"Reference","title":"PowerLaws.jl","text":"","category":"section"},{"location":"reference/","page":"Reference","title":"Reference","text":"Documentation for PowerLaws.jl","category":"page"},{"location":"reference/","page":"Reference","title":"Reference","text":"Modules = [PowerLaws]","category":"page"},{"location":"reference/#PowerLaws.DistributionComparison","page":"Reference","title":"PowerLaws.DistributionComparison","text":"DistributionComparison\n\nFields\n\ndata: array which according which should be distributions compared\nloglikelihoodratio: log likelihood ratio of data\nsig_level: sigma level\nxmin: smallest element which was used for comparing distributions\nVteststat: Vuong test statistic\nVpval: p-value from Vuong test\nVpreffdistr: preffered distibution according to Vuong test\nC_b: Clarke test b value - sum of possitive values in log likehood ratio\nCpval: p-value from Clarke test\nCpreffdistr: preffered distibution according to Clarke test\n\n\n\n\n\n","category":"type"},{"location":"reference/#PowerLaws.DistributionComparison-Tuple{ContinuousPowerLaw, Type{<:Distributions.Distribution{Distributions.Univariate, Distributions.Continuous}}, AbstractArray}","page":"Reference","title":"PowerLaws.DistributionComparison","text":"DistributionComparison\n\nThis function calculate Vuong test and Clarke test for non nested distributions. This is necessary since it is possible to fit power law distribution to any data set. Function was implemented according to this Non nested model selection for spatial count regression models with application to health insurance.\n\nArguments\n\nd1: First distribution to be compared.\nd2: Second distribution to be compared.\ndata: Data to be compared.\nsig_level: Significance level. Default is 0.05.\n\nReturns\n\nDistributionComparison: Struct containing all necessary information about comparison.\n\n\n\n\n\n","category":"method"},{"location":"reference/#Distributions.fit_mle-Tuple{Type{ContinuousPowerLaw}, AbstractArray{<:Real}}","page":"Reference","title":"Distributions.fit_mle","text":"fit_mle(::Type{ContinuousPowerLaw}, x::AbstractArray{<:Real})\n\nFit a ContinuousPowerLaw distribution to the data using maximum likelihood estimation (MLE). The x_min value is the minimum value of the data.\n\n\n\n\n\n","category":"method"},{"location":"reference/#Distributions.fit_mle-Tuple{Type{DiscretePowerLaw}, AbstractArray{<:Real}}","page":"Reference","title":"Distributions.fit_mle","text":"fit_mle(::Type{DiscretePowerLaw}, x::AbstractArray{<:Real})\n\nFits a discrete power law distribution to the data using an approximation to the maximum likelihood estimation (MLE). The x_min value is the minimum value of the data.\n\n\n\n\n\n","category":"method"},{"location":"reference/#PowerLaws.array_bins-Union{Tuple{AbstractArray{T}}, Tuple{T}} where T<:Real","page":"Reference","title":"PowerLaws.array_bins","text":"array_bins(arr::AbstractArray{T})::Dict{T,Int64} where {T <: Real}\n\nCreate a dictionary from a sorted array arr where the keys are the unique elements  and the values are the indices at which these elements first appear in the array.\n\n\n\n\n\n","category":"method"},{"location":"reference/#PowerLaws.bootstrap-Tuple{AbstractArray, Distributions.UnivariateDistribution}","page":"Reference","title":"PowerLaws.bootstrap","text":"bootstrap\n\nBootstrap method for estimating the parameters of a power law distribution. To quantify the uncertainty in our estimate for xmin you can use bootstrap method. More information can be found in this document Power-law distributions in empirical data.\n\nArguments\n\ndata::AbstractArray: Data to be tested.\nd::UnivariateDistribution: Distribution to be tested (ContinuousPowerLaw or DiscretePowerLaw).\nno_of_sims::Int64: Number of simulations. Default is 10.\nxmins::AbstractArray: Array of xmins to be tested, default is data.\nxmax::Int64: Maximum value of data to be considered. Default is 1e5.\nseed::Int64: Seed for random number generator. Default is 0.\n\nReturns\n\nstatistic::Array{Tuple{UnivariateDistribution, Float64}}: Array of tuples containing the distribution and the Kolmogorov-Smirnov distance between the data and the distribution.\n\n\n\n\n\n","category":"method"},{"location":"reference/#PowerLaws.bootstrap_p-Tuple{AbstractArray, Distributions.UnivariateDistribution}","page":"Reference","title":"PowerLaws.bootstrap_p","text":"bootstrap_p\n\nPerforms a bootstrapping hypothesis test to determine whether a power law distribution is plausible. Inspired by R poweRlaw documentation.\n\nArguments\n\ndata::AbstractArray: Data to be tested.\nd::UnivariateDistribution: Distribution to be tested (ContinuousPowerLaw or DiscretePowerLaw).\nno_of_sims::Int64: Number of simulations. Default is 10.\nxmins::AbstractArray: Array of xmins to be tested, default is data.\nxmax::Int64: Maximum value of data to be considered. Default is 1e5.\nseed::Int64: Seed for random number generator. Default is 0.\n\nReturns\n\nstatistic::Array{Tuple{UnivariateDistribution, Float64}}: Array of tuples containing the distribution and the Kolmogorov-Smirnov distance between the data and the distribution.\nP::Float64: p-value of the hypothesis test.\n\n\n\n\n\n","category":"method"},{"location":"reference/#PowerLaws.estimate_parameters-Tuple{AbstractArray, Type{ContinuousPowerLaw}}","page":"Reference","title":"PowerLaws.estimate_parameters","text":"estimate_parameters\n\nEstimate x_min and α for a given data set with respect to the Kolmogorov-Smirnov test.\n\nParameters\n\ndata::AbstractArray: Array of data which should be fit to a distribution.\ndistribution::Type: Distribution type, i.e. ContinuousPowerLaw or DiscretePowerLaw.\nxmins::AbstractArray: If not specified, all unique values in data are taken as possible x_mins. If specified, only values in xmins are considered when finding the best x_min.\nxmax::Int64: Maximum value considered in calculations. Values above xmax are not considered in, for example, calculating the Kolmogorov-Smirnov test.\n\nReturns\n\nbest_fit::distribution: A distribution of the type distribution fitted to the given parameters.\nKS::Float64: The Kolmogorov-Smirnov distance between the data and the fitted distribution.\n\n\n\n\n\n","category":"method"},{"location":"reference/#PowerLaws.kolmogorov_smirnov_test","page":"Reference","title":"PowerLaws.kolmogorov_smirnov_test","text":"kolmogorovsmirnovtest\n\nCalculate Kolmogorov Smirnov test on given data and distribution.\n\nArguments\n\ndat::AbstractArray: Data to be tested.\nd::UnivariateDistribution: Distribution to be tested (ContinuousPowerLaw or DiscretePowerLaw).\nxmin::Number: Minimum value of data to be considered.\nxmax::Int64: Maximum value of data to be considered. Default is 1e5.\n\nReturns\n\nKS::Float64: Kolmogorov-Smirnov distance between the data and the distribution. It is the maximum absolute difference between the empirical and theoretical cumulative distribution functions (`cdf\n\n\n\n\n\n","category":"function"},{"location":"#PowerLaws.jl","page":"Home","title":"PowerLaws.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Specification of discrete and continuous power-law distributions according to Distributions.jl. This package is implemented according to POWER-LAW DISTRIBUTIONS IN EMPIRICAL DATA","category":"page"},{"location":"","page":"Home","title":"Home","text":"This is a maintained fork of the archive PowerLaws.jl. Currently, the functionality replicates the original package, but with some compatibility upgrades, and modernised package structure. This will probably change in the future. Braking changes are however introduced by renaming structs and functions:","category":"page"},{"location":"","page":"Home","title":"Home","text":"powerlaw_con -> ContinuousPowerLaw\npowerlaw_dis -> DiscretePowerLaw\nestimate_xmin -> estimate_parameters\ncompare_distributions -> DistributionComparison","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This package is not (yet) registered, but can be installed with the following command:","category":"page"},{"location":"","page":"Home","title":"Home","text":"using Pkg\nPkg.add(url=\"https://www.github.com/TimVWese/PowerLaws.jl\")","category":"page"},{"location":"#Functionality","page":"Home","title":"Functionality","text":"","category":"section"},{"location":"#Discrete-Power-Law-Distribution","page":"Home","title":"Discrete Power-Law Distribution","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"DiscretePowerLaw(α, θ) represents a discrete power-law distribution with negative exponent α and minimum value θ. The probability mass function is given by","category":"page"},{"location":"","page":"Home","title":"Home","text":"p(x) = begincases fracx^-alphazeta(alpha theta)  x geq theta  0  x  theta endcases","category":"page"},{"location":"#Continuous-Power-Law-Distribution","page":"Home","title":"Continuous Power-Law Distribution","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"ContinuousPowerLaw(α, θ) represents a continuous power-law distribution with negative exponent α and minimum value θ. The probability density function is given by","category":"page"},{"location":"","page":"Home","title":"Home","text":"p(x) = begincases frac(alpha - 1)theta left(fracxthetaright)^-alpha  x geq theta  0  x  theta endcases","category":"page"},{"location":"","page":"Home","title":"Home","text":"Inpired by python powerlaw package and R poweRlaw package","category":"page"}]
}
