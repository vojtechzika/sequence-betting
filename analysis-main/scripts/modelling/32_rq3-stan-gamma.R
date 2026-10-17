# ============================================================
# 32_rq3-stan-gamma.R
#
# PURPOSE
#   Fits the Gamma-only model for RQ3 (certainty-equivalent welfare loss from
#   suboptimal betting) for every treatment and sample.
#   Gamma distribution for the strictly positive welfare losses, without a
#   hurdle. Preregistered for samples in which zero losses are rare (model-
#   implied probability of a zero loss below 0.05 in the hurdle-Gamma
#   model); trials with a zero loss are not used.
#   The script reports convergence and sample-size checks and takes no
#   action on them: which RQ3 model is used is decided by selecting its
#   folder under the key "rq3" in selected-models.csv.
#
# SAMPLES AND FITTING
#   Samples, welfare-loss construction, outputs and notes are described in
#   _rq3-fit.R and _rq3-welfare-loss.R, which this script calls.
#
# INPUT
#   model "mpl": r-draws-<tr>.rds, a-star/a-star-draws-<tr>-full.rds,
#                sample/analysis-trials.csv, sample/sample-definitions.csv,
#                drift/drift-decisions.rds
#   scripts/stan/rq3-gamma.stan
#
# OUTPUT (model folder path_mod/rq3-gamma/)
#   fit-, pid-levels-, seq-levels- and data-<tr>-<sample>.rds,
#   diagnostics.csv, inputs.csv, run-info.txt
# ============================================================

rq3_stan_gamma <- function(cfg) {
  source(here::here("scripts", "modelling", "_rq3-fit.R"))
  rq3_fit(cfg, list(
    folder     = "rq3-gamma",
    desc       = "Gamma-only",
    stan_file  = here::here("scripts", "stan", "rq3-gamma.stan"),
    settings   = cfg$model$stan$rq3_gamma,
    likelihood = "gamma"
  ))
}
