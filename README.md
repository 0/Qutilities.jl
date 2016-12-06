# Qutilities.jl

Assorted utilities for quantum information.

Tested with Julia 0.5.


## Installation

1. `Pkg.clone("https://github.com/0/Qutilities.jl.git")`


## Examples

```julia
using Qutilities

rho = [[1,0,0,0] [0,3,3,0] [0,3,3,0] [0,0,0,1]]/8.
#-> 0.125  0.0    0.0    0.0
#   0.0    0.375  0.375  0.0
#   0.0    0.375  0.375  0.0
#   0.0    0.0    0.0    0.125

ptrace(rho)
#-> 0.5  0.0
#   0.0  0.5
ptranspose(rho)
#-> 0.125  0.0    0.0    0.375
#   0.0    0.375  0.0    0.0
#   0.0    0.0    0.375  0.0
#   0.375  0.0    0.0    0.125

@printf "%f" purity(rho)
#-> 0.593750

@printf "%f" S_renyi(rho, 0)
#-> 2.000000
@printf "%f" S_vn(rho)
#-> 1.061278
@printf "%f" S_renyi(rho)
#-> 0.752072
@printf "%f" S_renyi(rho, Inf)
#-> 0.415037

@printf "%f" mutinf(rho)
#-> 0.938722
@printf "%f" concurrence(rho)
#-> 0.500000
@printf "%f" formation(rho)
#-> 0.354579
@printf "%f" negativity(rho)
#-> 0.584963
```


## License

Provided under the terms of the MIT license.
See `LICENSE` for more information.
