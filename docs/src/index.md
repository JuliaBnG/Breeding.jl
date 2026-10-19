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

# Pedigree and records in one table (Mrode Example 4.1). Animal IDs are
# row numbers; founders 1–3 have no records.
ped = DataFrame(
    calf = 1:8,
    sex  = [missing, missing, missing, "M", "F", "F", "M", "M"],
    sire = [0, 0, 0, 1, 3, 1, 4, 3],
    dam  = [0, 0, 0, 0, 2, 2, 5, 6],
    wwg  = [missing, missing, missing, 4.5, 2.9, 3.9, 3.5, 5.0],
)
rec = dropmissing(ped, :wwg)

# Incidence matrices for fixed sex and random animal effects
X, levels = incidence(rec.sex; levels = ["M", "F"])
Z = incidence(rec.calf, nrow(ped))

# A⁻¹ from RelationshipMatrices.jl; only the sire and dam columns are used
Ai = ainv(ped)

# Random genetic effect (σ²a = 20.0, σ²e = 40.0, α = 2.0)
a = RandomEffect(Z, Ai, 20.0; name = "animal")

# Assemble Henderson's mixed model equations
m = mme(rec.wwg, X, [a]; σ²e = 40.0)

# Solve for BLUE (fixed) and BLUP (EBV)
sol = solve_mme(m)
b, (u,) = solutions(m, sol)

# Breeding values, prediction error variances and reliabilities
ped.ebv = vec(u)
ped.pev = vec(pev(m, 1))
ped.r2  = vec(reliability(m, 1))
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
