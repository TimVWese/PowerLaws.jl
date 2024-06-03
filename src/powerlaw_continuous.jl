struct ContinuousPowerLawDistribution <: ContinuousUnivariateDistribution
    α::Float64
    θ::Float64

    function ContinuousPowerLawDistribution(α::Real, θ::Real)
        if !(α > zero(α) && θ > zero(θ))
            throw(ArgumentError("Both shape and scale parameters must be positive."))
        end
        new(α, θ)
    end
    ContinuousPowerLawDistribution(α::Real) = ContinuousPowerLawDistribution(α, 1.0)
    ContinuousPowerLawDistribution() = new(1.0, 1.0)
end
@distr_support ContinuousPowerLawDistribution d.θ Inf
#### Parameters

shape(d::ContinuousPowerLawDistribution) = d.α
scale(d::ContinuousPowerLawDistribution) = d.θ
params(d::ContinuousPowerLawDistribution) = (d.α, d.θ)
#### Statistics

mean(d::ContinuousPowerLawDistribution) = ((α, θ) = params(d); α > 2.0 ? θ *((α - 1.0)/(α - 2.0)) : Inf)

median(d::ContinuousPowerLawDistribution) = ((α, θ) = params(d);α > 1.0 ? 2.0^(1.0/(α - 1.0)) * θ : NaN)
mode(d::ContinuousPowerLawDistribution) = d.θ
function var(d::ContinuousPowerLawDistribution)
    (α, θ) = params(d)
    α > 3.0 ? (θ^2 * (α-1)) / ((α - 2.0)^2 * (α - 3.0)) : Inf
end

function skewness(d::ContinuousPowerLawDistribution)
    α = shape(d)
    α > 4.0 ? ((2.0 * (α)) / (α - 4.0)) * sqrt((α - 3.0) / (α-1)) : NaN
end

function kurtosis(d::ContinuousPowerLawDistribution)
    α = shape(d)
    α > 5.0 ? (6.0 * ((α-1)^3 + (α-1)^2 - 6.0 * (α-1) - 2.0)) / ((α-1) * (α - 4.0) * (α - 5.0)) : NaN
end

entropy(d::ContinuousPowerLawDistribution) = ((α, θ) = params(d); log(θ / (α-1)) + 1.0 / (α-1) + 1.0)


#### Evaluation

function pdf(d::ContinuousPowerLawDistribution, x::Real)
    (α, θ) = params(d)
    x >= θ ? ((α-1.0)/θ) * ((x/θ)^(-α)) : 0.0
end

function pdf(d::ContinuousPowerLawDistribution, x::AbstractArray)
    (α, θ) = params(d)
    cons = ((α-1.0)/θ)
    pdfs = [num >= θ ? cons * ((x/θ)^(-α)) : 0.0 for num in x]
    return pdfs
end

function logpdf(d::ContinuousPowerLawDistribution, x::Real)
    (α, θ) = params(d)
    x >= θ ? log(α-1.0) - log(θ) - α * log(x/θ) : -Inf
end

function logpdf(d::ContinuousPowerLawDistribution, x::AbstractArray{<:Real})
    (α, θ) = params(d)
    l_const = log(α-1.0) - log(θ)
    lpdfs = [num >= θ ? l_const - α * log(num/θ) : -Inf for num in x]
    return lpdfs
end


function ccdf(d::ContinuousPowerLawDistribution, x::Float64)
    (α, θ) = params(d)
    x >= θ ? (x/θ)^(-α +1.0) : 1.0
end

cdf(d::ContinuousPowerLawDistribution, x::Float64) = 1.0 - ccdf(d, x)

logccdf(d::ContinuousPowerLawDistribution, x::Float64) = log(ccfd(d,x))

logcdf(d::ContinuousPowerLawDistribution, x::Float64) = log(cdf(d, x))

cquantile(d::ContinuousPowerLawDistribution, p::Float64) = d.θ *((p) ^ (-1.0 / (d.α - 1.0)))
quantile(d::ContinuousPowerLawDistribution, p::Float64) = cquantile(d, 1.0 - p)


#### Sampling

rand(d::ContinuousPowerLawDistribution) = quantile(d,rand())


## Fitting

function fit_mle(::Type{ContinuousPowerLawDistribution}, x::AbstractArray{<:Real})
    θ = minimum(x)

    n = length(x)
    lθ = log(θ)
    temp1 = sum(log.(x) .- lθ)
    α = 1.0+n*(temp1)^(-1.0)

    return ContinuousPowerLawDistribution(α, θ)
end
