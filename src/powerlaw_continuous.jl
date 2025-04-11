struct ContinuousPowerLaw <: ContinuousUnivariateDistribution
    α::Float64
    θ::Float64

    function ContinuousPowerLaw(α::Real, θ::Real)
        if !(α > zero(α) && θ > zero(θ))
            throw(ArgumentError("Both shape and scale parameters must be positive."))
        end
        new(α, θ)
    end
    ContinuousPowerLaw(α::Real) = ContinuousPowerLaw(α, 1.0)
    ContinuousPowerLaw() = new(1.0, 1.0)
end
@distr_support ContinuousPowerLaw d.θ Inf
#### Parameters

shape(d::ContinuousPowerLaw) = d.α
scale(d::ContinuousPowerLaw) = d.θ
params(d::ContinuousPowerLaw) = (d.α, d.θ)
#### Statistics

mean(d::ContinuousPowerLaw) = ((α, θ) = params(d);
α > 2.0 ? θ * ((α - 1.0) / (α - 2.0)) : Inf)

median(d::ContinuousPowerLaw) = ((α, θ) = params(d);
α > 1.0 ? 2.0^(1.0 / (α - 1.0)) * θ : NaN)
mode(d::ContinuousPowerLaw) = d.θ
function var(d::ContinuousPowerLaw)
    (α, θ) = params(d)
    α > 3.0 ? (θ^2 * (α - 1)) / ((α - 2.0)^2 * (α - 3.0)) : Inf
end

function skewness(d::ContinuousPowerLaw)
    α = shape(d)
    α > 4.0 ? ((2.0 * (α)) / (α - 4.0)) * sqrt((α - 3.0) / (α - 1)) : NaN
end

function kurtosis(d::ContinuousPowerLaw)
    α = shape(d)
    α > 5.0 ? (6.0 * ((α - 1)^3 + (α - 1)^2 - 6.0 * (α - 1) - 2.0)) / ((α - 1) * (α - 4.0) * (α - 5.0)) : NaN
end

entropy(d::ContinuousPowerLaw) = ((α, θ) = params(d);
log(θ / (α - 1)) + 1.0 / (α - 1) + 1.0)


#### Evaluation

function pdf(d::ContinuousPowerLaw, x::Real)
    (α, θ) = params(d)
    x >= θ ? ((α - 1.0) / θ) * ((x / θ)^(-α)) : 0.0
end

function pdf(d::ContinuousPowerLaw, x::AbstractArray)
    (α, θ) = params(d)
    cons = ((α - 1.0) / θ)
    pdfs = [num >= θ ? cons * ((x / θ)^(-α)) : 0.0 for num in x]
    return pdfs
end

function logpdf(d::ContinuousPowerLaw, x::Real)
    (α, θ) = params(d)
    x >= θ ? log(α - 1.0) - log(θ) - α * log(x / θ) : -Inf
end

function logpdf(d::ContinuousPowerLaw, x::AbstractArray{<:Real})
    (α, θ) = params(d)
    l_const = log(α - 1.0) - log(θ)
    lpdfs = [num >= θ ? l_const - α * log(num / θ) : -Inf for num in x]
    return lpdfs
end


function ccdf(d::ContinuousPowerLaw, x::Float64)
    (α, θ) = params(d)
    x >= θ ? (x / θ)^(-α + 1.0) : 1.0
end

cdf(d::ContinuousPowerLaw, x::Float64) = 1.0 - ccdf(d, x)

logccdf(d::ContinuousPowerLaw, x::Float64) = log(ccfd(d, x))

logcdf(d::ContinuousPowerLaw, x::Float64) = log(cdf(d, x))

cquantile(d::ContinuousPowerLaw, p::Float64) = d.θ * ((p)^(-1.0 / (d.α - 1.0)))
quantile(d::ContinuousPowerLaw, p::Float64) = cquantile(d, 1.0 - p)


#### Sampling

rand(d::ContinuousPowerLaw) = quantile(d, rand())


## Fitting
"""
    fit_mle(::Type{ContinuousPowerLaw}, x::AbstractArray{<:Real})

Fit a `ContinuousPowerLaw` distribution to the data using maximum likelihood
estimation (MLE). The `x_min` value is the minimum value of the data.
"""
function fit_mle(::Type{ContinuousPowerLaw}, x::AbstractArray{<:Real})
    θ = minimum(x)
    n = length(x)
    lθ = log(θ)
    temp1 = sum(log.(x) .- lθ)
    α = 1.0 + n * (temp1)^(-1.0)

    return ContinuousPowerLaw(α, θ)
end

