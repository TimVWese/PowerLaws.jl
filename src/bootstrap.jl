"""
  bootstrap

Bootstrap method for estimating the parameters of a power law distribution.
To quantify the uncertainty in our estimate for xmin you can use bootstrap method. More information can be found in this document [Power-law distributions in empirical data](http://arxiv.org/pdf/0706.1062v2.pdf).

# Arguments

- `data::AbstractArray`: Data to be tested.
- `d::UnivariateDistribution`: Distribution to be tested (`ContinuousPowerLaw` or `DiscretePowerLaw`).
- `no_of_sims::Int64`: Number of simulations. Default is `10`.
- `xmins::AbstractArray`: Array of xmins to be tested, default is `data`.
- `xmax::Int64`: Maximum value of data to be considered. Default is `1e5`.
- `seed::Int64`: Seed for random number generator. Default is `0`.

# Returns
- `statistic::Array{Tuple{UnivariateDistribution, Float64}}`: Array of tuples containing the distribution and the Kolmogorov-Smirnov distance between the data and the distribution.
"""
function bootstrap(data::AbstractArray, d::UnivariateDistribution; no_of_sims::Int64=10, xmins::AbstractArray=[], xmax::Int64=round(Int, 1e5), seed::Int64=0)
  n = length(data)

  seed == 0 ? Random.seed!() : Random.seed!(seed)
  statistic = Array{Tuple{UnivariateDistribution,Float64}}(undef, no_of_sims)
  for i = 1:no_of_sims
    sim_data = sample(data, n, replace=true)
    statistic[i] = estimate_parameters(sim_data, typeof(d), xmins=xmins, xmax=xmax)
  end
  return statistic
end

function bootstrap(data::AbstractArray, distribution::Type{ContinuousPowerLaw}; no_of_sims::Int64=10, xmins::AbstractArray=[], xmax::Int64=round(Int, 1e5), seed::Int64=0)
  d, ks = estimate_parameters(data, distribution, xmins=xmins, xmax=xmax)
  bootstrap(data, d, no_of_sims=no_of_sims, xmins=xmins, xmax=xmax, seed=seed)
end

function bootstrap(data::AbstractArray, distribution::Type{DiscretePowerLaw}; no_of_sims::Int64=10, xmins::AbstractArray=[], xmax::Int64=round(Int, 1e5), seed::Int64=0)
  d, ks = estimate_parameters(data, distribution, xmins=xmins, xmax=xmax)
  bootstrap(data, d, no_of_sims=no_of_sims, xmins=xmins, xmax=xmax, seed=seed)
end

"""
  bootstrap_p

Performs a bootstrapping hypothesis test to determine whether a power law distribution is plausible.
Inspired by R [poweRlaw](http://arxiv.org/pdf/1407.3492v1.pdf) documentation.

# Arguments
- `data::AbstractArray`: Data to be tested.
- `d::UnivariateDistribution`: Distribution to be tested (`ContinuousPowerLaw` or `DiscretePowerLaw`).
- `no_of_sims::Int64`: Number of simulations. Default is `10`.
- `xmins::AbstractArray`: Array of xmins to be tested, default is `data`.
- `xmax::Int64`: Maximum value of data to be considered. Default is `1e5`.
- `seed::Int64`: Seed for random number generator. Default is `0`.

# Returns
- `statistic::Array{Tuple{UnivariateDistribution, Float64}}`: Array of tuples containing the distribution and the Kolmogorov-Smirnov distance between the data and the distribution.
- `P::Float64`: p-value of the hypothesis test.
"""
function bootstrap_p(data::AbstractArray, d::UnivariateDistribution; no_of_sims::Int64=10, xmins::AbstractArray=[], xmax::Int64=round(Int, 1e5), seed::Int64=0)
  sort_data = sort(data)
  α, θ = params(d)
  n = length(sort_data)
  tail_indx = findfirst(sort_data, θ)
  tail_p = length(sort_data[tail_indx:end]) / n
  KS_stat = kolmogorov_smirnov_test(sort_data[tail_indx:end], d)

  P = 0
  statistic = Array(Tuple{typeof(d),Float64}, no_of_sims)
  seed == 0 ? Random.seed!() : Random.seed!(seed)
  for i = 1:no_of_sims
    n1 = sum(map(x -> x > tail_p, rand(n)))
    n2 = n - n1
    sim_data = sample(sort_data[1:tail_indx-1], n1, replace=true)
    append!(sim_data, rand(d, n2))
    statistic[i] = estimate_parameters(sim_data, typeof(d), xmins=xmins, xmax=xmax)
    if (KS_stat <= statistic[i][2])
      P += 1
    end
  end
  return statistic, (P / no_of_sims)
end

function bootstrap_p(data::AbstractArray, distribution::Type{ContinuousPowerLaw}; no_of_sims::Int64=10, xmins::AbstractArray=[], xmax::Int64=round(Int, 1e5), seed::Int64=0)
  d, ks = estimate_parameters(data, distribution, xmins=xmins, xmax=xmax)
  bootstrap_p(data, d, no_of_sims=no_of_sims, xmins=xmins, xmax=xmax, seed=seed)
end

function bootstrap_p(data::AbstractArray, distribution::Type{DiscretePowerLaw}; no_of_sims::Int64=10, xmins::AbstractArray=[], xmax::Int64=round(Int, 1e5), seed::Int64=0)
  d, ks = estimate_parameters(data, distribution, xmins=xmins, xmax=xmax)
  bootstrap_p(data, d, no_of_sims=no_of_sims, xmins=xmins, xmax=xmax, seed=seed)
end
