function init_xmins(data::AbstractArray, xmins::AbstractArray, xmax::Int64)
  sorted_data = sort(data)
  bins_data = array_bins(sorted_data)
  if (xmins == [])
    xmins = sort(unique(sorted_data))
    if (xmins[end] > xmax)
      max_indx = findlast(x -> x <= xmax, xmins)
      xmins = xmins[1:max_indx]
    end
    xmins = xmins[1:end]
    if xmins[1] <= 0
      println("removing elements smaller than 0")
      xmins = xmins[findfirst(x -> x > 0, xmins):end]
    end
  else
    real_xmins = [haskey(bins_data, x) for x in xmins]
    xmins = xmins[real_xmins]
  end
  return sorted_data, bins_data, xmins
end

function _estimate_parameters(sorted_data::AbstractArray, bins_data::Dict{<:Real,Int64}, distribution::DataType, xmins::AbstractArray, xmax::Int64)
  if (length(xmins) == 0)
    throw(ArgumentError("No xmins"))
  end

  min_dist = Inf
  best_fit = Union{}

  for xmin in xmins
    fit_data = sorted_data[bins_data[xmin]:end]
    f = fit(distribution, fit_data)

    negloglike(alpha) = begin
      d = distribution(alpha[1], f.θ)
      r = -sum(logpdf(d, fit_data))
      if (Inf == r || -Inf == r)
        r = 1e12
      end
      return r
    end
    try
      opt_alfa = fminbox(DifferentiableFunction(negloglike), [f.α], [1.0], [Inf])
      f = distribution(opt_alfa.minimum[1], f.θ)
    catch
      #if fminbox throws error it means that function cannot be optimized
    end
    if (f.α == Inf)
      continue
    end
    d = kolmogorov_smirnov_test(fit_data, f, xmin, xmax)

    if ((min_dist > d))
      best_fit = f
      min_dist = d
    end
  end
  return best_fit, min_dist
end

"""
  estimate_parameters

Estimate `x_min` and `α` for a given data set with respect to the [Kolmogorov-Smirnov test](https://www.encyclopediaofmath.org/index.php/Kolmogorov-Smirnov_test).

# Parameters
- `data::AbstractArray`: Array of data which should be fit to a distribution.
- `distribution::Type`: Distribution type, i.e. `ContinuousPowerLaw` or `DiscretePowerLaw`.
- `xmins::AbstractArray`: If not specified, all unique values in `data` are taken as possible `x_min`s. If specified, only values in `xmins` are considered when finding the best `x_min`.
- `xmax::Int64`: Maximum value considered in calculations. Values above `xmax` are not considered in, for example, calculating the Kolmogorov-Smirnov test.

# Returns
- `best_fit::distribution`: A distribution of the type `distribution` fitted to the given parameters.
- `KS::Float64`: The Kolmogorov-Smirnov distance between the data and the fitted distribution.
"""
function estimate_parameters(data::AbstractArray, distribution::Type{ContinuousPowerLaw}; xmins::AbstractArray=[], xmax::Int64=round(Int, 1e5))

  sorted_data, bins_data, xmins = init_xmins(data, xmins, xmax)

  _estimate_parameters(sorted_data, bins_data, distribution, xmins, xmax)
end

function estimate_parameters(data::AbstractArray, distribution::Type{DiscretePowerLaw}; xmins::AbstractArray=[], xmax::Int64=round(Int, 1e5))
  if !all(data .== floor.(data))
    throw(ArgumentError("Data should be discreate. Use round or floor function."))
  end

  sorted_data, bins_data, xmins = init_xmins(data, xmins, xmax)
  if xmins[1] < 1
    xmins = xmins[findfirst(x -> x >= 1):end]
  end

  _estimate_parameters(sorted_data, bins_data, distribution, xmins, xmax)
end

"""
    array_bins(arr::AbstractArray{T})::Dict{T,Int64} where {T <: Real}

Create a dictionary from a sorted array `arr` where the keys are the unique elements 
and the values are the indices at which these elements first appear in the array.
"""
function array_bins(arr::AbstractArray{T})::Dict{T,Int64} where {T<:Real}
  bins = Dict{T,Int64}()
  for (i, num) in enumerate(arr)
    get!(bins, num, i)
  end
  return bins
end
