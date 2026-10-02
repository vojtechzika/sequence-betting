# ============================================================
# 41_rq4_stan.R
#   RQ4 (side choice conditional on betting): Stan estimation
#
# SCOPE
#   Full treatment sample - no normative-better restriction and no
#   HL-consistency filter. Side choices do not depend on r_i, so
#   HL quality is not a measurement criterion here. Using the full
#   sample also ensures that the RQ4 posterior draws (mu_h[s],
#   hbar, mu_h_i[i]) are on the same population basis as EX1,
#   which constructs the GHI from these draws over all participants.
#
# PER TREATMENT:
#   1) Loads master_sequences.csv from path_src
#   2) Keeps betting trials only (stake > 0)
#   3) Uses cfg$design$seq$side_labels: heads = "H", tails = "T", nobet = "NB"
#   4) Constructs h = 1(side == "H"), 0(side == "T")
#   5) Reads drift_decisions.rds from path_mod; passes drift flags to Stan
#   6) Fits stan/rq4_side.stan (per treatment)
#   7) Saves Stan artifacts to path_mod:
#        rq4_fit_sequences_<tr>_full.rds
#        rq4_pid_levels_<tr>_full.rds
#        rq4_seq_levels_<tr>_full.rds
#
# NO CONFIRMATORY-SUBSET PASS
#   RQ4 has no normative-better restriction. A single "full" tag is
#   produced per treatment. Downstream scripts (42_, 43_, 54_) all
#   consume the full-sample fit.
#
# NO consistent_only FILTER
#   cfg$run$consistent_only applies to RQ1-3 only (those analyses
#   depend on r_i from the HL task). HL inconsistency is irrelevant
#   to the quality of side-choice observations.
#
# STAN MODEL: stan/rq4.stan
#   - Bernoulli-logit hierarchical model
#   - Participant RE u_i ~ N(0, sigma_u), non-centred
#   - Sequence RE beta_s ~ N(0, sigma_s), non-centred + sum-to-zero
#   - Optional drift (linear or categorical) via include_drift / drift_type
#   - Generated quantities: y_rep[T] (PPC / overdispersion check),
#     mu_h[S], hbar, mu_h_i[N], mu_h_sample[S], hbar_sample
#   - mc_M passed as data from cfg$model$ppc$rq4_mc_M (MC integration
#     resolution for population quantities mu_h and hbar)
#
# CONFIG KEYS USED
#   cfg$run$data_folder            -- resolves path_src, path_mod (via 00_setup.R)
#   cfg$run$treatment              -- c("m25", "m19")
#   cfg$run$seed
#   cfg$run$overwrite_models       -- governs should_skip()
#   cfg$design$seq$side_labels     -- heads/tails/nobet string labels
#   cfg$design$drift$params        -- prior_gamma_sd per drift type
#   cfg$design$drift$rq4$outcome   -- "side" (used to build drift_decisions key)
#   cfg$model$stan$rq4             -- iter/warmup/chains/adapt_delta/treedepth
#   cfg$model$ppc$rq4_mc_M         -- MC draw count for mu_h / hbar integration
#                                     (add to 02_models.R ppc list; default 2000)
#
# INPUTS
#   path_src/master_sequences.csv
#   path_mod/drift_decisions.rds
#
# OUTPUTS (consumed by 42_rq4_diagnostics.R, 43_rq4_tables.R, 54_ex1_2_stan_ghi.R)
#   path_mod/rq4_fit_sequences_<tr>_full.rds
#   path_mod/rq4_pid_levels_<tr>_full.rds
#   path_mod/rq4_seq_levels_<tr>_full.rds
#
# CALL ORDER IN PIPELINE (03_stage_analysis.R)
#   rq4_stan(cfg)         -- fits Bernoulli model (this script)
#   rq4_diagnostics(cfg)  -- PPC, overdispersion check, drift diagnostics;
#                            determines whether Beta-Binomial robustness
#                            fit is required (rq4_side_bb.stan);
#                            reads cfg$design$rq4$ppc_overdisp_cut
#   rq4_tables(cfg)       -- posterior summaries, sequence/participant tables;
#                            reads cfg$design$rq4$delta for tolerance bands
#   ex1_2_stan_ghi(cfg)   -- GHI construction from rq4 posterior draws
# ============================================================


library(data.table)
library(rstan)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

rq4_stan <- function(cfg) {
  
  # ----------------------------
  # Run-level settings
  # ----------------------------
  seed   <- as.integer(cfg$run$seed)
  design <- cfg$design
  model  <- cfg$model
  
  tr_vec <- unique(as.character(cfg$run$treatment))
  stopifnot(length(tr_vec) > 0L, all(nzchar(tr_vec)))
  
  # ----------------------------
  # Side labels  (design$seq$side_labels: heads="H", tails="T", nobet="NB")
  # ----------------------------
  stopifnot(!is.null(design$seq$side_labels),
            is.list(design$seq$side_labels))
  
  lab_heads <- as.character(design$seq$side_labels$heads)
  lab_tails <- as.character(design$seq$side_labels$tails)
  lab_nobet <- as.character(design$seq$side_labels$nobet)
  
  stopifnot(length(lab_heads) == 1L, nzchar(lab_heads))
  stopifnot(length(lab_tails) == 1L, nzchar(lab_tails))
  stopifnot(length(lab_nobet) == 1L, nzchar(lab_nobet))
  stopifnot(length(unique(c(lab_heads, lab_tails, lab_nobet))) == 3L)
  
  msg("RQ4 side labels: heads='", lab_heads,
      "' tails='", lab_tails,
      "' nobet='", lab_nobet, "'")
  
  # ----------------------------
  # Paths  (path_src and path_mod are global, set by 00_setup.R from
  #         cfg$run$data_folder; no per-dataset helper needed here)
  # ----------------------------
  infile    <- file.path(path_src, "master_sequences.csv")
  stan_file <- here::here("stan", "rq4.stan")
  
  stopifnot(file.exists(infile))
  stopifnot(file.exists(stan_file))
  
  dir.create(path_mod, showWarnings = FALSE, recursive = TRUE)
  
  # ----------------------------
  # Stan sampling settings
  # cfg$model$stan$rq4 is a flat list: iter / warmup / chains /
  # adapt_delta / treedepth  (see 02_models.R)
  # ----------------------------
  st <- model$stan$rq4
  stopifnot(!is.null(st))
  
  iter_val        <- as.integer(st$iter)
  warmup_val      <- as.integer(st$warmup)
  chains_val      <- as.integer(st$chains)
  adapt_delta_val <- as.numeric(st$adapt_delta)
  treedepth_val   <- as.integer(st$treedepth)
  
  stopifnot(iter_val > 0L, warmup_val >= 0L, chains_val > 0L)
  stopifnot(is.finite(adapt_delta_val),
            adapt_delta_val > 0, adapt_delta_val < 1)
  stopifnot(treedepth_val > 0L)
  
  # MC integration resolution for population quantities (mu_h, hbar).
  # Stored in cfg$model$ppc$rq4_mc_M alongside the other ppc keys
  # (rq1_k, rq3_k, rq4_k already present in 02_models.R).
  # Add  rq4_mc_M = 2000L  to model_cfg()$ppc to pin the value explicitly.
  mc_M <- as.integer(model$ppc$rq4_mc_M %||% 2000L)
  stopifnot(is.finite(mc_M), mc_M >= 100L)
  
  msg("RQ4 Stan settings:",
      " iter=", iter_val,
      " warmup=", warmup_val,
      " chains=", chains_val,
      " adapt_delta=", adapt_delta_val,
      " treedepth=", treedepth_val,
      " mc_M=", mc_M,
      " seed=", seed)
  
  # ----------------------------
  # Drift decisions
  # drift_decisions.rds written by 02_stage_indices.R.
  # Key convention mirrors RQ1-3: <tr>_full_<outcome>
  # outcome = design$drift$rq4$outcome = "side"  (from 01_design.R)
  # -> keys will be e.g. "m25_full_side", "m19_full_side"
  # ----------------------------
  f_drift <- file.path(path_mod, "drift_decisions.rds")
  stopifnot(file.exists(f_drift))
  drift_decisions <- readRDS(f_drift)
  
  outcome_key <- as.character(design$drift$rq4$outcome)   # "side"
  stopifnot(length(outcome_key) == 1L, nzchar(outcome_key))
  
  get_drift_cfg <- function(tr) {
    key <- paste(tr, "full", outcome_key, sep = "_")
    res <- drift_decisions[[key]]
    if (is.null(res)) {
      msg("RQ4 drift: no entry for key='", key, "' -- defaulting to no drift")
      return(list(drift = FALSE, type = "none"))
    }
    list(drift = isTRUE(res$drift), type = as.character(res$drift_type))
  }
  
  drift_type_to_int <- function(type_str) {
    switch(type_str,
           "none"        = 0L,
           "linear"      = 1L,
           "categorical" = 2L,
           {
             msg("RQ4 drift: unrecognised drift_type='", type_str,
                 "' -- defaulting to 0 (none)")
             0L
           })
  }
  
  # ----------------------------
  # Load master data
  # ----------------------------
  dt <- fread(infile, encoding = "UTF-8")
  
  required <- c("pid", "treat", "seq", "stake", "side", "block")
  missing  <- setdiff(required, names(dt))
  if (length(missing) > 0L) {
    stop(
      "master_sequences.csv missing columns: ",
      paste(missing, collapse = ", "),
      "\nExpected at least: ", paste(required, collapse = ", ")
    )
  }
  
  dt[, pid   := as.character(pid)]
  dt[, treat := as.character(treat)]
  dt[, seq   := as.character(seq)]
  dt[, stake := as.numeric(stake)]
  dt[, side  := as.character(side)]
  dt[, block := as.integer(block)]
  dt[is.na(stake), stake := 0]
  
  # Centred block index for drift.
  # block in {1,2,3,4} -> block_c in {-1.5,-0.5,0.5,1.5}
  # Matches design$drift$params$linear$block_centers.
  # Stan recovers block index as: blk = (int)(block_c + 2.5) -> {1,2,3,4}
  dt[, block_c := as.numeric(block) - 2.5]
  
  # ----------------------------
  # Main loop over treatments
  # Stan model compiled lazily on the first treatment that actually
  # needs fitting so no compilation occurs when all outputs are fresh.
  # ----------------------------
  sm <- NULL
  
  for (tr in tr_vec) {
    
    f_fit <- file.path(path_mod, paste0("rq4_fit_sequences_", tr, "_full.rds"))
    f_pid <- file.path(path_mod, paste0("rq4_pid_levels_",    tr, "_full.rds"))
    f_seq <- file.path(path_mod, paste0("rq4_seq_levels_",    tr, "_full.rds"))
    
    if (should_skip(
      paths = c(f_fit, f_pid, f_seq),
      cfg   = cfg,
      type  = "model",
      label = paste0("RQ4 Stan (", tr, ")")
    )) next
    
    if (is.null(sm)) {
      msg("RQ4 Stan: compiling ", stan_file)
      sm <- rstan::stan_model(stan_file)
    }
    
    # ---- Subset to treatment ----
    d <- dt[treat == tr]
    if (nrow(d) == 0L) {
      msg("RQ4 Stan: no data for tr='", tr, "' -- skipping")
      next
    }
    
    # ---- Betting trials only ----
    d <- d[is.finite(stake) & stake > 0]
    if (nrow(d) == 0L) {
      msg("RQ4 Stan: no betting trials for tr='", tr, "' -- skipping")
      next
    }
    
    # ---- Validate side labels on betting trials ----
    bad <- d[!(side %in% c(lab_heads, lab_tails))]
    if (nrow(bad) > 0L) {
      stop(
        "RQ4: side values outside {heads='", lab_heads,
        "', tails='", lab_tails, "'} found after stake > 0",
        " (tr='", tr, "').\n",
        "Examples: ", paste(unique(head(bad$side, 10L)), collapse = ", ")
      )
    }
    
    # ---- Outcome ----
    d[, h := as.integer(side == lab_heads)]   # 1 = Heads, 0 = Tails
    
    # ---- Factor levels ----
    pid_levels <- sort(unique(d$pid))
    seq_levels <- sort(unique(d$seq))
    
    d[, pid_i := match(pid, pid_levels)]
    d[, sid_s := match(seq, seq_levels)]
    stopifnot(!anyNA(d$pid_i), !anyNA(d$sid_s))
    
    # ---- Drift configuration ----
    drift_cfg      <- get_drift_cfg(tr)
    include_drift  <- as.integer(isTRUE(drift_cfg$drift))
    drift_type_int <- drift_type_to_int(drift_cfg$type)
    
    # prior_gamma_sd read from design$drift$params (01_design.R):
    #   linear:      design$drift$params$linear$prior_gamma_sd      = 0.3
    #   categorical: design$drift$params$categorical$prior_gamma_sd = 0.3
    # When include_drift == 0 the value is still required because the drift
    # parameters are always declared in the Stan model to keep sampler geometry
    # stable; 0.3 is used as the prior SD in that prior-only regime.
    prior_gamma_sd <- if (include_drift == 1L) {
      gsd <- design$drift$params[[drift_cfg$type]]$prior_gamma_sd
      if (is.null(gsd)) {
        msg("RQ4 drift: prior_gamma_sd not found in design$drift$params",
            " for type='", drift_cfg$type, "' -- using 0.3")
        0.3
      } else {
        as.numeric(gsd)
      }
    } else {
      0.3
    }
    
    msg("RQ4 Stan: fitting tr=", tr,
        " | N=", length(pid_levels),
        " | S=", length(seq_levels),
        " | T=", nrow(d),
        " | drift=", drift_cfg$type,
        " (include_drift=", include_drift,
        ", drift_type=", drift_type_int, ")",
        " | prior_gamma_sd=", prior_gamma_sd,
        " | mc_M=", mc_M)
    
    # ---- Stan data list ----
    data_list <- list(
      N              = length(pid_levels),
      S              = length(seq_levels),
      T              = nrow(d),
      pid            = as.integer(d$pid_i),
      sid            = as.integer(d$sid_s),
      h              = as.integer(d$h),
      include_drift  = include_drift,
      drift_type     = drift_type_int,
      block_c        = as.numeric(d$block_c),
      prior_gamma_sd = prior_gamma_sd,
      mc_M           = mc_M
    )
    
    # ---- Sample ----
    fit <- rstan::sampling(
      sm,
      data    = data_list,
      iter    = iter_val,
      warmup  = warmup_val,
      chains  = chains_val,
      seed    = seed,
      control = list(
        adapt_delta   = adapt_delta_val,
        max_treedepth = treedepth_val
      )
    )
    
    if (length(rstan::get_sampler_params(fit, inc_warmup = FALSE)) == 0L) {
      stop("RQ4: Stan produced no samples for tr='", tr,
           "'. Fit will NOT be saved.")
    }
    
    # ---- Save ----
    saveRDS(fit,        f_fit)
    saveRDS(pid_levels, f_pid)
    saveRDS(seq_levels, f_seq)
    
    msg("Saved: ", f_fit)
    msg("Saved: ", f_pid)
    msg("Saved: ", f_seq)
  }
  
  invisible(NULL)
}