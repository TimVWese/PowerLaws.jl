struct DiscretePowerLaw <: DiscreteUnivariateDistribution
    α::Float64
    θ::Float64

    function DiscretePowerLaw(α::Real, θ::Real)
        if !(α > zero(α) && θ > zero(θ))
            throw(ArgumentError("Both shape and scale parameters must be positive."))
        end
        new(α, θ)
    end
    DiscretePowerLaw(α::Real) = DiscretePowerLaw(α, 1.0)
    DiscretePowerLaw() = new(1.0, 1.0)
end
@distr_support DiscretePowerLaw d.θ Inf
#### Parameters

shape(d::DiscretePowerLaw) = d.α
scale(d::DiscretePowerLaw) = d.θ
params(d::DiscretePowerLaw) = (d.α, d.θ)


#### Evaluation

function pdf(d::DiscretePowerLaw, x::Real)
    (α, θ) = params(d)
    x >= θ ?  x ^ (-α) / zeta(α, θ) : 0.0
end

function pdf(d::DiscretePowerLaw, x::AbstractArray{<:Real})
  (α, θ) = params(d)
  z = zeta(α, θ)
  pdfs = [num >= θ ?  num ^ (-α) / z : 0.0 for num in x]
  return pdfs
end

function logpdf(d::DiscretePowerLaw, x::Real)
    (α, θ) = params(d)
    x >= θ ? -log(zeta(α, θ)) - α * log(x) : -Inf
end
function logpdf(d::DiscretePowerLaw, x::AbstractArray{<:Real})
  (α, θ) = params(d)
  log_zeta = -log(zeta(α, θ))
  lpdfs = [num >= θ ? log_zeta - α * log(num) : -Inf for num in x]
  return lpdfs
end

function ccdf(d::DiscretePowerLaw, x::Float64)
    (α, θ) = params(d)
    x >= θ ? zeta(α, x) / zeta(α, θ) : 1.0
end

function ccdf(d::DiscretePowerLaw, x::AbstractArray{<:Real})
    (α, θ) = params(d)
    z = zeta(α, θ)
    ccdfs = [num >= θ ? zeta(α, num) / z : 1.0 for num in x]
    return ccdfs
end

cdf(d::DiscretePowerLaw, x::Float64) = 1.0 - ccdf(d, x)
cdf(d::DiscretePowerLaw, x::AbstractArray{<:Real}) = 1.0 .- ccdf(d, x)

logccdf(d::DiscretePowerLaw, x::Float64) = log(ccfd(d,x))
logccdf(d::DiscretePowerLaw, x::AbstractArray{<:Real}) = log.(ccfd(d,x))

logcdf(d::DiscretePowerLaw, x::Float64) = log(cdf(d, x))
logcdf(d::DiscretePowerLaw, x::AbstractArray{<:Real}) = log.(cdf(d, x))

cquantile(d::DiscretePowerLaw, p::Float64) = (d.θ -0.5) *((p) ^ (-1.0 / (d.α - 1.0))) + 0.5
quantile(d::DiscretePowerLaw, p::Float64) = cquantile(d, 1.0 - p)


#### Sampling

rand(d::DiscretePowerLaw) = floor(quantile(d,rand()))



## Fitting
#The discrete MLE of α, is not available, instead we use the approximation
function fit_mle(::Type{DiscretePowerLaw}, x::AbstractArray{<:Real})
    θ = minimum(x)
    n = length(x)
    lθ = log(θ-0.5)
    temp1 = sum(log.(x) .- lθ)
    α = 1.0 +n*(temp1)^(-1.0)

    return DiscretePowerLaw(α, θ)
end
