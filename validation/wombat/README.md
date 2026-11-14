# REML against WOMBAT

`validate.jl` simulates one pedigree and two data sets, writes them in
WOMBAT's formats, runs `wombat` (it must be on the `PATH`) and
`Breeding.reml` from the same starting values, and prints estimates,
approximate standard errors and maximised log-likelihoods side by side.
The simulation model, sizes and seed are stated at the top of the script.

    julia -t 2 --project=<environment with Breeding> validation/wombat/validate.jl [workdir]

`results.txt` holds the output for the tagged version. Both programs
report the restricted log-likelihood without constants, so the values are
directly comparable. WOMBAT prints six significant digits, which bounds
the agreement that can be shown.

This compares two programs on the same data. It says nothing about the
sampling properties of the estimates: there is one replicate per model.
Random regression, heterogeneous residual classes and maternal effects
are covered by the package tests against dense likelihoods, not here.
