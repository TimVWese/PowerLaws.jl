struct con_powerlaw <: ContinuousUnivariateDistribution
    α::Float64
    θ::Float64

    function con_powerlaw(α::Real, θ::Real)
        if !(α > zero(α) && θ > zero(θ))
            throw(ArgumentError("Both shape and scale parameters must be positive."))
        end
        new(α, θ)
    end
    con_powerlaw(α::Real) = con_powerlaw(α, 1.0)
    con_powerlaw() = new(1.0, 1.0)
end
@distr_support con_powerlaw d.θ Inf
#### Parameters

shape(d::con_powerlaw) = d.α
scale(d::con_powerlaw) = d.θ
params(d::con_powerlaw) = (d.α, d.θ)
#### Statistics

mean(d::con_powerlaw) = ((α, θ) = params(d); α > 2.0 ? θ *((α - 1.0)/(α - 2.0)) : Inf)

median(d::con_powerlaw) = ((α, θ) = params(d);α > 1.0 ? 2.0^(1.0/(α - 1.0)) * θ : NaN)
mode(d::con_powerlaw) = d.θ
function var(d::con_powerlaw)
    (α, θ) = params(d)
    α > 3.0 ? (θ^2 * (α-1)) / ((α - 2.0)^2 * (α - 3.0)) : Inf
end

function skewness(d::con_powerlaw)
    α = shape(d)
    α > 4.0 ? ((2.0 * (α)) / (α - 4.0)) * sqrt((α - 3.0) / (α-1)) : NaN
end

function kurtosis(d::con_powerlaw)
    α = shape(d)
    α > 5.0 ? (6.0 * ((α-1)^3 + (α-1)^2 - 6.0 * (α-1) - 2.0)) / ((α-1) * (α - 4.0) * (α - 5.0)) : NaN
end

entropy(d::con_powerlaw) = ((α, θ) = params(d); log(θ / (α-1)) + 1.0 / (α-1) + 1.0)


#### Evaluation

function pdf(d::con_powerlaw, x::Real)
    (α, θ) = params(d)
    x >= θ ? ((α-1.0)/θ) * ((x/θ)^(-α)) : 0.0
end

function pdf(d::con_powerlaw, x::AbstractArray)
    (α, θ) = params(d)
    cons = ((α-1.0)/θ)
    pdfs = [num >= θ ? cons * ((x/θ)^(-α)) : 0.0 for num in x]
    return pdfs
end

function logpdf(d::con_powerlaw, x::Real)
    (α, θ) = params(d)
    x >= θ ? log(α-1.0) - log(θ) - α * log(x/θ) : -Inf
end

function logpdf(d::con_powerlaw, x::AbstractArray{<:Real})
    (α, θ) = params(d)
    l_const = log(α-1.0) - log(θ)
    lpdfs = [num >= θ ? l_const - α * log(num/θ) : -Inf for num in x]
    return lpdfs
end


function ccdf(d::con_powerlaw, x::Float64)
    (α, θ) = params(d)
    x >= θ ? (x/θ)^(-α +1.0) : 1.0
end

cdf(d::con_powerlaw, x::Float64) = 1.0 - ccdf(d, x)

logccdf(d::con_powerlaw, x::Float64) = log(ccfd(d,x))

logcdf(d::con_powerlaw, x::Float64) = log(cdf(d, x))

cquantile(d::con_powerlaw, p::Float64) = d.θ *((p) ^ (-1.0 / (d.α - 1.0)))
quantile(d::con_powerlaw, p::Float64) = cquantile(d, 1.0 - p)


#### Sampling

rand(d::con_powerlaw) = quantile(d,rand())


## Fitting

function fit_mle(::Type{con_powerlaw}, x::AbstractArray{<:Real})
    θ = minimum(x)

    n = length(x)
    lθ = log(θ)
    temp1 = sum(log.(x) .- lθ)
    α = 1.0+n*(temp1)^(-1.0)

    return con_powerlaw(α, θ)
end
