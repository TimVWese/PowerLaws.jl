module PowerLaws
  using StatsBase
  using Distributions
  using Optim
  using Compat
  using Random
  using SpecialFunctions: zeta
  import Distributions: rand, shape, pdf, ccdf, scale, params, cdf, cquantile, 
    quantile, fit_mle, mean, median, var, skewness, mode, kurtosis, logccdf, 
    logpdf, logcdf, entropy, @distr_support
  import Base: minimum, maximum, show

  export ContinuousPowerLawDistribution,DiscretePowerLawDistribution,Kolmogorov_smirnov_test,estimate_xmin,bootstrap,bootstrap_p,DistributionComparison

  include("powerlaw_discrete.jl")
  include("powerlaw_continuous.jl")
  include("estimate_xmin.jl")
  include("Kolmogorov-smirnov_test.jl")
  include("bootstrap.jl")
  include("compare_distributions.jl")

end
