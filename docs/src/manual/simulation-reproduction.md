# Breeding Simulation & Gene Drop

Forward simulation lets selection strategies, mating schemes and
genomic-prediction designs be compared before they are used in the
field.

`Breeding.jl` works on genotypes in the formats of
[`BnGStructs.jl`](https://juliabng.github.io/BnGStructs/), e.g. founder
populations simulated with
[`FisherWright.jl`](https://juliabng.github.io/FisherWright/): it drops
recombined gametes, samples marker panels and QTL, and generates true
breeding values and phenotypes.

---

## Genome & Recombination Models

Recombination is modelled across multi-chromosome genomes using Poisson-distributed crossover
events with uniform physical placement within each chromosome:

- [`sum_map`](@ref): chromosome lengths in Morgans (1 Morgan per
  `M = 1e8` bp by default) and cumulative end positions in bp.
- [`sumMap`](@ref): per-chromosome summary of a linkage map (length in
  Morgans, number of loci, index of the first locus).
- `crossovers`: Sample the number of crossovers per chromosome $N_c \sim \operatorname{Poisson}(\lambda)$ and their locations.
- `gamete`: Construct a recombined haploid gamete from a diploid parent's two homologous haplotypes.

---

## Gene Dropping (`gene_drop`)

Given parental genotypes and a mating design $\mathbf{M}$ (pairs of sire and dam IDs), [`gene_drop`](@ref)
simulates the transmission of whole genomes to the progeny generation in place:

```julia
# Drop gametes into the offspring generation:
xo = gene_drop(parent_genotypes, offspring_genotypes, mating_matrix, map_summary, snp_positions)
```

The returned `DataFrame` (`generation`, `parent`, `ihap`, `crossovers`)
records, for every gamete, the parent, the haplotype it starts on and
its crossover positions, so identity-by-descent segments can be traced.
The `generation` column is a placeholder (`0`) for the caller to fill.

---

## Sampling Founders & Marker Panels

Starting from coalescent or forward-simulated populations (e.g. from `FisherWright.jl`):

- [`sampleID`](@ref): sample individuals, recompute allele frequencies
  (written to `lmp.frq` of the input map) and keep loci passing the MAF
  threshold.
- [`sampleLoci`](@ref): Partition available polymorphic loci into exclusive or overlapping SNP-chip panels and QTL subsets with defined MAF constraints.

```julia
# Sample founders:
xy_founders, lmp_founders = sampleID(mxy, lmp, 0.05, 500)

# Sample SNP panel and QTL panel:
chip = SNPSet("chip", 10_000, 0.05, true)
qtl  = SNPSet("qtl", 1_000, 0.01, true)
mxy_panel, lmp_panel = sampleLoci(xy_founders, lmp_founders, chip, qtl)
```

---

## QTL Effects & Phenotypic Realization

### Simulating QTL Effects (`eqtl`)

[`eqtl`](@ref) samples additive QTL effects normalized such that the resulting genetic values
have zero mean and target genetic standard deviation $\sigma_a$:

```julia
# Single trait:
eqtl(mxy_panel, lmp_panel, trait)

# Multiple correlated traits with genetic correlation matrix R:
eqtl(mxy_panel, lmp_panel, R_matrix, traits)
```

### True Breeding Values (`tbv!`)

Evaluates true breeding values for all animals from their QTL genotypes and assigned allele substitution effects:

```julia
tbv!(ped, mxy_panel, lmp_panel, traits)
```

### Phenotype Realization (`phenotype!`)

Adds environmental deviations $e_i \sim N(0, \sigma^2_e)$ based on trait heritability $h^2$,
and applies category thresholds for threshold traits:

```julia
phenotype!(ped, traits)
```
