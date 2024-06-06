# PowerLaws.jl

Specification of discrete and continuous power-law distributions according to [`Distributions.jl`](https://github.com/JuliaStats/Distributions.jl).
This package is implemented according to [POWER-LAW DISTRIBUTIONS IN EMPIRICAL DATA](http://arxiv.org/pdf/0706.1062v2.pdf)

This is a maintained fork of the archive [PowerLaws.jl](https://github.com/johnybx/PowerLaws.jl).
Currently, the functionality replicates the original package, but with some compatibility upgrades, and modernised package structure.
This will probably change in the future.
Braking changes are however introduced by renaming structs and functions:
- `powerlaw_con` -> `ContinuousPowerLaw`
- `powerlaw_dis` -> `DiscretePowerLaw`
- `estimate_xmin` -> `estimate_parameters`
- `compare_distributions` -> `DistributionComparison`

## Installation
This package is not (yet) registered, but can be installed with the following command:
```julia
using Pkg
Pkg.add(url="https://www.github.com/TimVWese/PowerLaws.jl")
```

## Functionality

### Discrete Power-Law Distribution

`DiscretePowerLaw(α, θ)` represents a discrete power-law distribution with negative exponent `α` and minimum value `θ`.
The probability mass function is given by

``p(x) = \\begin{cases} \\frac{x^{-\\alpha}}{\\zeta(\\alpha, \\theta)} & x \\geq \\theta \\\\ 0 & x < \\theta \\end{cases}``

### Continuous Power-Law Distribution

`ContinuousPowerLaw(α, θ)` represents a continuous power-law distribution with negative exponent `α` and minimum value `θ`.
The probability density function is given by

``p(x) = \\begin{cases} \\frac{(\\alpha - 1)}{\\theta} \\left(\\frac{x}{\\theta}\\right)^{-\\alpha} & x \\geq \\theta \\\\ 0 & x < \\theta \\end{cases}``



**Inpired by python** [powerlaw package](https://pypi.python.org/pypi/powerlaw)
**and R** [poweRlaw package](http://arxiv.org/pdf/1407.3492v1.pdf)

