struct DiscretePowerLawDistribution <: DiscreteUnivariateDistribution
    α::Float64
    θ::Float64

    function DiscretePowerLawDistribution(α::Real, θ::Real)
        if !(α > zero(α) && θ > zero(θ))
            throw(ArgumentError("Both shape and scale parameters must be positive."))
        end
        new(α, θ)
    end
    DiscretePowerLawDistribution(α::Real) = DiscretePowerLawDistribution(α, 1.0)
    DiscretePowerLawDistribution() = new(1.0, 1.0)
end
@distr_support DiscretePowerLawDistribution d.θ Inf
#### Parameters

shape(d::DiscretePowerLawDistribution) = d.α
scale(d::DiscretePowerLawDistribution) = d.θ
params(d::DiscretePowerLawDistribution) = (d.α, d.θ)


#### Evaluation

function pdf(d::DiscretePowerLawDistribution, x::Real)
    (α, θ) = params(d)
    x >= θ ?  x ^ (-α) / zeta(α, θ) : 0.0
end

function pdf(d::DiscretePowerLawDistribution, x::AbstractArray{<:Real})
  (α, θ) = params(d)
  z = zeta(α, θ)
  pdfs = [num >= θ ?  num ^ (-α) / z : 0.0 for num in x]
  return pdfs
end

function logpdf(d::DiscretePowerLawDistribution, x::Real)
    (α, θ) = params(d)
    x >= θ ? -log(zeta(α, θ)) - α * log(x) : -Inf
end
function logpdf(d::DiscretePowerLawDistribution, x::AbstractArray{<:Real})
  (α, θ) = params(d)
  log_zeta = -log(zeta(α, θ))
  lpdfs = [num >= θ ? log_zeta - α * log(num) : -Inf for num in x]
  return lpdfs
end

function ccdf(d::DiscretePowerLawDistribution, x::Float64)
    (α, θ) = params(d)
    x >= θ ? zeta(α, x) / zeta(α, θ) : 1.0
end

function ccdf(d::DiscretePowerLawDistribution, x::AbstractArray{<:Real})
    (α, θ) = params(d)
    z = zeta(α, θ)
    ccdfs = [num >= θ ? zeta(α, num) / z : 1.0 for num in x]
    return ccdfs
end

cdf(d::DiscretePowerLawDistribution, x::Float64) = 1.0 - ccdf(d, x)
cdf(d::DiscretePowerLawDistribution, x::AbstractArray{<:Real}) = 1.0 .- ccdf(d, x)

logccdf(d::DiscretePowerLawDistribution, x::Float64) = log(ccfd(d,x))
logccdf(d::DiscretePowerLawDistribution, x::AbstractArray{<:Real}) = log.(ccfd(d,x))

logcdf(d::DiscretePowerLawDistribution, x::Float64) = log(cdf(d, x))
logcdf(d::DiscretePowerLawDistribution, x::AbstractArray{<:Real}) = log.(cdf(d, x))

cquantile(d::DiscretePowerLawDistribution, p::Float64) = (d.θ -0.5) *((p) ^ (-1.0 / (d.α - 1.0))) + 0.5
quantile(d::DiscretePowerLawDistribution, p::Float64) = cquantile(d, 1.0 - p)


#### Sampling

rand(d::DiscretePowerLawDistribution) = floor(quantile(d,rand()))



## Fitting
#The discrete MLE of α, is not available, instead we use the approximation
function fit_mle(::Type{DiscretePowerLawDistribution}, x::AbstractArray{<:Real})
    θ = minimum(x)
    n = length(x)
    lθ = log(θ-0.5)
    temp1 = sum(log.(x) .- lθ)
    α = 1.0 +n*(temp1)^(-1.0)

    return DiscretePowerLawDistribution(α, θ)
end
