# Breeding.jl

`Breeding.jl` provides mixed-model genetic evaluation, genomic
prediction, variance-component estimation and forward breeding
simulation in Julia. The methods follow Mrode and Pocrnic (2023) and the
original papers cited in each docstring. Relationship matrices come from
[`RelationshipMatrices.jl`](https://juliabng.github.io/RelationshipMatrices/)
and genotype containers from
[`BnGStructs.jl`](https://juliabng.github.io/BnGStructs/).

## Installation

```julia
using Pkg
Pkg.add(url = "https://github.com/JuliaBnG/Breeding.jl.git")
```

Or when registered:

```julia
using Pkg
Pkg.add("Breeding")
```

Then load the package:

```julia
using Breeding
using DataFrames
using RelationshipMatrices
```

## Quick start

Fit Henderson's animal model on weaning weight gain with fixed sex effect and random additive genetic effect:

```julia
using Breeding
using DataFrames
using RelationshipMatrices
using SparseArrays

# Pedigree and records (Mrode Example 4.1)
sire = [0, 0, 0, 1, 3, 1, 4, 3]
dam  = [0, 0, 0, 0, 2, 2, 5, 6]
calf = [4, 5, 6, 7, 8]
sex  = ["M", "F", "F", "M", "M"]
wwg  = [4.5, 2.9, 3.9, 3.5, 5.0]

# Setup incidence matrices
X, levels = incidence(sex; levels = ["M", "F"])
Z = incidence(calf, 8)

# Additive relationship matrix inverse A⁻¹ from RelationshipMatrices.jl
ped = DataFrame(sire = sire, dam = dam)
Ai = ainv(ped)

# Define random genetic effect (σ²a = 20.0, σ²e = 40.0, α = 2.0)
a = RandomEffect(Z, sparse(Ai), 20.0; name = "animal")

# Assemble Henderson's Mixed Model Equations
m = mme(wwg, X, [a]; σ²e = 40.0)

# Solve for BLUE (fixed) and BLUP (EBV)
sol = solve_mme(m)
b, (u,) = solutions(m, sol)

# Prediction error variance and reliabilities
P = pev(m, 1)
r2 = reliability(m, 1)
```

## Contents

```@contents
Pages = [
    "manual/mixed-models.md",
    "manual/iterative-solvers.md",
    "manual/genomic-prediction.md",
    "manual/variance-components.md",
    "manual/threshold-models.md",
    "manual/simulation-reproduction.md",
    "api.md",
]
Depth = 2
```
