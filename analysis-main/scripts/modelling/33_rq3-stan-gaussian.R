# ============================================================
# 33_rq3-stan-gaussian.R
#
# PURPOSE
#   Fits the Gaussian model for RQ3 (certainty-equivalent welfare loss from
#   suboptimal betting) for every treatment and sample.
#   Gaussian diagnostic alternative with the same participant and sequence
#   structure. It ignores the boundary at zero and the right skew of the
#   losses; it is used to check the sign and ordering of the sequence
#   effects across likelihoods and as a convergence fallback.
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
#   scripts/stan/rq3-gaussian.stan
#
# OUTPUT (model folder path_mod/rq3-gaussian/)
#   fit-, pid-levels-, seq-levels- and data-<tr>-<sample>.rds,
#   diagnostics.csv, inputs.csv, run-info.txt
# ============================================================

rq3_stan_gaussian <- function(cfg) {
  source(here::here("scripts", "modelling", "_rq3-fit.R"))
  rq3_fit(cfg, list(
    folder     = "rq3-gaussian",
    desc       = "Gaussian",
    stan_file  = here::here("scripts", "stan", "rq3-gaussian.stan"),
    settings   = cfg$model$stan$rq3_gaussian,
    likelihood = "gaussian"
  ))
}
