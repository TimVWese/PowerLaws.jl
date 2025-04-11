# PowerLaws.jl

[![Build Status](https://github.com/TimVWese/PowerLaws.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/TimVWese/PowerLaws.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://timvwese.github.io/PowerLaws.jl/dev/)


This package is implemented according to [POWER-LAW DISTRIBUTIONS IN EMPIRICAL DATA](http://arxiv.org/pdf/0706.1062v2.pdf).


This is a maintained fork of the archive [PowerLaws.jl](https://github.com/johnybx/PowerLaws.jl).
Currently, the functionality replicates the original package, but with some compatibility upgrades, and modernised package structure.
This will probably change in the future.
Breaking changes are however introduced by renaming structs and functions:
- `powerlaw_con` -> `ContinuousPowerLaw`
- `powerlaw_dis` -> `DiscretePowerLaw`
- `estimate_xmin` -> `estimate_parameters`
- `compare_distributions` -> `DistributionComparison`

This package is not (yet) registered, but can be installed with the following command:
```julia
using Pkg
Pkg.add(url="https://www.github.com/TimVWese/PowerLaws.jl")
```

**Inpired by python** [powerlaw package](https://pypi.python.org/pypi/powerlaw)
**and R** [poweRlaw package](http://arxiv.org/pdf/1407.3492v1.pdf)
