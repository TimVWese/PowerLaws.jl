"""
  kolmogorov_smirnov_test

Calculate [Kolmogorov Smirnov test](https://www.encyclopediaofmath.org/index.php/Kolmogorov-Smirnov_test) on given data and distribution.

# Arguments
- `dat::AbstractArray`: Data to be tested.
- `d::UnivariateDistribution`: Distribution to be tested (`ContinuousPowerLaw` or `DiscretePowerLaw`).
- `xmin::Number`: Minimum value of data to be considered.
- `xmax::Int64`: Maximum value of data to be considered. Default is `1e5`.

# Returns
- `KS::Float64`: Kolmogorov-Smirnov distance between the data and the distribution. It is the maximum absolute difference between the empirical and theoretical cumulative distribution functions (`cdf
"""
function kolmogorov_smirnov_test(dat::AbstractArray, d::ContinuousUnivariateDistribution, xmin::Number, xmax::Int64=round(Int, 1e5))
  data = sort(dat)
  n = float(length(data))
  max_indx = findlast(x -> x <= xmax, data)
  min_indx = findfirst(x -> x >= xmin, data)
  data = data[min_indx:max_indx]

  act_cdf = collect(0:length(data)-1) / n
  thr_cdf = map(Base.Fix1(cdf, d), float(data))
  KS = maximum(abs.(act_cdf - thr_cdf))
  return KS
end

function kolmogorov_smirnov_test(dat::AbstractArray, d::DiscreteUnivariateDistribution, xmin::Int64, xmax::Int64=round(Int, 1e5))
  _, xmin = params(d)
  data = Int64.(round.(sort(dat)))
  max_indx = findlast(x -> x <= xmax, data)
  min_indx = findfirst(x -> x >= xmin, data)
  data = data[min_indx:max_indx]

  thr_cdf = cdf(d, collect(xmin:data[end]) .+ 1)
  occurence = counts(data)[Int(xmin)-(data[1]-1):end]
  act_cdf = occurence / sum(occurence)
  act_cdf = cumsum(act_cdf)
  KS = maximum(abs.(act_cdf - thr_cdf))
  return KS
end

#helper function
function create_histogram(x::AbstractArray)
  h = zeros(typeof(x[1]), 0)
  max = 0
  d = 0
  for i = 1:length(x)
    if max < x[i]
      for j = 1:(x[i]-max)
        push!(h, 0)
      end
      h[x[i]] = 1
      max = x[i]
    else
      h[x[i]] += 1
    end
  end
  return h
end
