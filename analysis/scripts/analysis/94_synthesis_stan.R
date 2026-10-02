# ============================================================
# 94_synthesis_stan.R
#
# PURPOSE
#   Fits the joint hierarchical bet/side-choice model used in the
#   Synthesis section: a two-equation model (extensive margin = bet
#   vs. no-bet; side choice | bet == 1) sharing one correlated
#   participant-level random-effects prior across both equations.
#   The cross-equation correlations this recovers are the formal
#   test of the "dissociation" hypothesis -- whether betting and
#   side choice are driven by a shared participant-level tendency
#   or by distinct ones (see Sec:Synthesis).
#
#   This file contains ONLY the fitting step. Diagnostics, tables,
#   participant-level heterogeneity classification, and figures all
#   read the saved .rds fit and live in 97_synthesis_stan_outputs.R.
#
# CUES (sequence-level, outcome-blind; built in synthesis_stan()):
#   -- betting equation: magnitude / presence only (unsigned) --
#   imb_mag = imbalance magnitude, |H count - T count| / 6, in [0,1]
#             (master_sequences.csv: imbalance_abs)
#   run_dum = 1 if seq ends in a run of >=2 identical symbols, else 0
#             (master_sequences.csv: has_terminal_run)
#   -- side-choice equation: signed --
#   run_h_s = signed terminal-run cue: +1 if seq ends in a run of >=2 H's,
#             -1 if >=2 T's, 0 if no terminal run
#             (master_sequences.csv: terminal_run)
#   imb_h_s = signed imbalance cue, (H count - T count) / 6
#             (master_sequences.csv: imbalance)
#   Additive specification (no run_h_s*imb_h_s interaction).
#
# INPUT
#   path_src/master_sequences.csv
#   stan/synthesis.stan
#
# OUTPUT
#   path_mod/synthesis_fit_<tr>.rds  -- fit + data_list + id maps
#   (everything downstream of the fit -- tables, diagnostics,
#   heterogeneity, figures -- is produced by 95_synthesis_stan_outputs.R)
#
# ============================================================
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
`%||%` <- function(a, b) if (is.null(a)) b else a
# ------------------------------------------------------------
# FIT
# ------------------------------------------------------------
synthesis_stan <- function(cfg) {
  
  seed   <- as.integer(cfg$run$seed)
  tr_vec <- unique(as.character(cfg$run$treatment))
  
  stan_file <- here::here("stan", "synthesis.stan")
  stopifnot(file.exists(stan_file))
  
  infile <- file.path(path_src, "master_sequences.csv")
  stopifnot(file.exists(infile))
  
  st <- cfg$model$stan$synthesis
  iter_val        <- as.integer(st$iter        %||% 8000)
  warmup_val      <- as.integer(st$warmup      %||% 1500)
  chains_val      <- as.integer(st$chains      %||% 4)
  adapt_delta_val <- as.numeric(st$adapt_delta %||% 0.97)
  treedepth_val   <- as.integer(st$treedepth   %||% 13)
  
  dt <- fread(infile, encoding = "UTF-8")
  stopifnot(all(c("pid", "treat", "seq", "side") %in% names(dt)))
  dt[, pid   := as.character(pid)]
  dt[, treat := as.character(treat)]
  dt[, seq   := as.character(seq)]
  dt[, side  := as.character(side)]
  
  # ---- sequence-level cue covariates (outcome-blind, purely combinatorial) ----
  # Read directly from master_sequences.csv:
  #   imbalance          = (H count - T count) / 6, in [-1, 1]  (SIGNED)
  #   imbalance_abs       = abs(imbalance), in [0, 1]
  #   terminal_run        = 1 if seq ends in a run of >=2 Heads,
  #                        -1 if >=2 Tails, 0 if no terminal run  (SIGNED)
  #   has_terminal_run    = 1 if terminal_run != 0, else 0  (UNSIGNED dummy)
  stopifnot(all(c("imbalance", "imbalance_abs", "terminal_run", "has_terminal_run") %in% names(dt)))
  dt[, run_h_s := terminal_run]
  dt[, imb_h_s := imbalance]
  dt[, run_dum := has_terminal_run]
  dt[, imb_mag := imbalance_abs]
  
  # ---- fit one treatment subset ----
  fit_one <- function(d0, tr) {
    
    f_fit <- file.path(path_mod, paste0("synthesis_fit_", tr, ".rds"))
    if (should_skip(f_fit, cfg, "model", paste0("Synthesis Stan (", tr, ")"))) return(invisible(NULL))
    if (nrow(d0) == 0) return(invisible(NULL))
    
    pid_levels <- sort(unique(d0$pid))
    seq_levels <- sort(unique(d0$seq))
    d0[, pid_i := match(pid, pid_levels)]
    d0[, sid_s := match(seq, seq_levels)]
    d0[, bet   := as.integer(side %in% c("H", "T"))]
    
    bet_d <- d0[side %in% c("H", "T")]
    bet_d[, side01 := as.integer(side == "H")]
    
    data_list <- list(
      P = length(pid_levels), S = length(seq_levels), N = nrow(d0),
      pid = as.integer(d0$pid_i), sid = as.integer(d0$sid_s),
      imb_mag = as.numeric(d0$imb_mag), run_dum = as.numeric(d0$run_dum),
      bet = as.integer(d0$bet),
      N_bet = nrow(bet_d),
      pid_s = as.integer(bet_d$pid_i), sid_s = as.integer(bet_d$sid_s),
      run_h_s = as.numeric(bet_d$run_h_s), imb_h_s = as.numeric(bet_d$imb_h_s),
      side = as.integer(bet_d$side01),
      prior_only = 0L
    )
    
    sm <- rstan::stan_model(stan_file)
    
    msg("Synthesis Stan: fitting tr=", tr,
        " | P=", data_list$P, " | S=", data_list$S,
        " | N=", data_list$N, " | N_bet=", data_list$N_bet)
    
    fit <- rstan::sampling(
      sm, data = data_list,
      iter = iter_val, warmup = warmup_val, chains = chains_val, seed = seed,
      control = list(adapt_delta = adapt_delta_val, max_treedepth = treedepth_val)
    )
    
    saveRDS(list(fit = fit, data_list = data_list,
                 pid_levels = pid_levels, seq_levels = seq_levels),
            f_fit)
    msg("Saved: ", f_fit)
    invisible(f_fit)
  }
  
  for (tr in tr_vec) {
    fit_one(dt[treat == tr], tr)
  }
  
  invisible(TRUE)
}