struct dis_powerlaw <: DiscreteUnivariateDistribution
    α::Float64
    θ::Float64

    function dis_powerlaw(α::Real, θ::Real)
        @assert α > zero(α) && θ > zero(θ)
        new(α, θ)
    end
    dis_powerlaw(α::Real) = dis_powerlaw(α, 1.0)
    dis_powerlaw() = new(1.0, 1.0)
end
@distr_support dis_powerlaw d.θ Inf
#### Parameters

shape(d::dis_powerlaw) = d.α
scale(d::dis_powerlaw) = d.θ
params(d::dis_powerlaw) = (d.α, d.θ)


#### Evaluation

function pdf(d::dis_powerlaw, x::Real)
    (α, θ) = params(d)
    x >= θ ?  x ^ (-α) / zeta(α, θ) : 0.0
end

function pdf(d::dis_powerlaw, x::AbstractArray{<:Real})
  (α, θ) = params(d)
  z = zeta(α, θ)
  pdfs = [num >= θ ?  num ^ (-α) / z : 0.0 for num in x]
  return pdfs
end

function logpdf(d::dis_powerlaw, x::Real)
    (α, θ) = params(d)
    x >= θ ? -log(zeta(α, θ)) - α * log(x) : -Inf
end
function logpdf(d::dis_powerlaw, x::AbstractArray{<:Real})
  (α, θ) = params(d)
  log_zeta = -log(zeta(α, θ))
  lpdfs = [num >= θ ? log_zeta - α * log(num) : -Inf for num in x]
  return lpdfs
end

function ccdf(d::dis_powerlaw, x::Float64)
    (α, θ) = params(d)
    x >= θ ? zeta(α, x) / zeta(α, θ) : 1.0
end

function ccdf(d::dis_powerlaw, x::AbstractArray{<:Real})
    (α, θ) = params(d)
    z = zeta(α, θ)
    ccdfs = [num >= θ ? zeta(α, num) / z : 1.0 for num in x]
    return ccdfs
end

cdf(d::dis_powerlaw, x::Float64) = 1.0 - ccdf(d, x)
cdf(d::dis_powerlaw, x::AbstractArray{<:Real}) = 1.0 .- ccdf(d, x)

logccdf(d::dis_powerlaw, x::Float64) = log(ccfd(d,x))
logccdf(d::dis_powerlaw, x::AbstractArray{<:Real}) = log.(ccfd(d,x))

logcdf(d::dis_powerlaw, x::Float64) = log(cdf(d, x))
logcdf(d::dis_powerlaw, x::AbstractArray{<:Real}) = log.(cdf(d, x))

cquantile(d::dis_powerlaw, p::Float64) = (d.θ -0.5) *((p) ^ (-1.0 / (d.α - 1.0))) + 0.5
quantile(d::dis_powerlaw, p::Float64) = cquantile(d, 1.0 - p)


#### Sampling

rand(d::dis_powerlaw) = floor(quantile(d,rand()))



## Fitting
#The discrete MLE of α, is not available, instead we use the approximation
function fit_mle(::Type{dis_powerlaw}, x::AbstractArray{<:Real})
    θ = minimum(x)
    n = length(x)
    lθ = log(θ-0.5)
    temp1 = sum(log.(x) .- lθ)
    α = 1.0 +n*(temp1)^(-1.0)

    return dis_powerlaw(α, θ)
end
