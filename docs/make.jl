using Documenter
using Breeding

makedocs(
    sitename = "Breeding.jl",
    authors = "Xijiang Yu",
    modules = [Breeding],
    checkdocs = :exports,
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
    ),
    pages = [
        "Home" => "index.md",
        "Manual" => [
            "Linear Mixed Models & BLUP" => "manual/mixed-models.md",
            "Iterative Solvers & Iteration-on-Data" => "manual/iterative-solvers.md",
            "Genomic Prediction & Bayesian Alphabet" => "manual/genomic-prediction.md",
            "Variance Component Estimation" => "manual/variance-components.md",
            "Threshold & Categorical Models" => "manual/threshold-models.md",
            "Breeding Simulation & Gene Drop" => "manual/simulation-reproduction.md",
        ],
        "API reference" => "api.md",
    ],
)

if !isempty(get(ENV, "DOCUMENTER_KEY", ""))
    deploydocs(
        repo = "github.com/JuliaBnG/Breeding.jl.git",
        deploy_repo = "github.com/JuliaBnG/juliabng.github.io.git",
        dirname = "Breeding",
        forcepush = true,
    )
elseif get(ENV, "GITHUB_ACTIONS", "") == "true"
    @warn "Skipping documentation deployment because DOCUMENTER_KEY is not configured."
end
