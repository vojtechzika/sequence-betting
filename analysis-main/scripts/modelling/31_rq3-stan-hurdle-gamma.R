# ============================================================
# 31_rq3-stan-hurdle-gamma.R
#
# PURPOSE
#   Fits the hurdle-Gamma model for RQ3 (certainty-equivalent welfare loss from
#   suboptimal betting) for every treatment and sample.
#   The preregistered primary model, used in every sample: zero welfare losses
#   (expected-utility optimal bets) are modelled by a hurdle with participant
#   and sequence effects, positive losses by a Gamma distribution. The model-
#   implied mean probability of a zero loss is written to hurdle-summary.csv,
#   which shows whether the hurdle is negligible (below 0.05) or substantial
#   in the sample. The Gamma-only model (32_rq3-stan-gamma.R) is a robustness
#   check for samples in which zero losses are rare.
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
#   scripts/stan/rq3-hurdle-gamma.stan
#
# OUTPUT (model folder path_mod/rq3-hurdle-gamma/)
#   fit-, pid-levels-, seq-levels- and data-<tr>-<sample>.rds,
#   diagnostics.csv, hurdle-summary.csv, inputs.csv, run-info.txt
# ============================================================

rq3_stan_hurdle_gamma <- function(cfg) {
  source(here::here("scripts", "modelling", "_rq3-fit.R"))
  rq3_fit(cfg, list(
    folder     = "rq3-hurdle-gamma",
    desc       = "hurdle-Gamma",
    stan_file  = here::here("scripts", "stan", "rq3-hurdle-gamma.stan"),
    settings   = cfg$model$stan$rq3_hurdle_gamma,
    likelihood = "hurdle_gamma"
  ))
}
