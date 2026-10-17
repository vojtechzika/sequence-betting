# ============================================================
# 91_synthesis-stan-joint.R
#
# PURPOSE
#   Fits the joint hierarchical model of the Synthesis section: a
#   two-equation model of the extensive margin (bet vs. no bet, all trials)
#   and of side choice conditional on betting (Heads vs. Tails), with one
#   correlated participant-level random-effects prior across both
#   equations (scripts/stan/synthesis-joint.stan). The cross-equation
#   correlations of the participant effects are the formal test of the
#   dissociation between betting and side choice. The script reports
#   convergence and posterior predictive checks and takes no action on
#   them: whether the model is used is decided by selecting its folder in
#   selected-models.csv.
#
#   Cues (sequence-level, outcome-blind, computed from the six symbols):
#     betting equation (unsigned)
#       imb_mag  = |#Heads - #Tails| / 6
#       run_dum  = 1 if the sequence ends in a run of at least two
#                  identical symbols, else 0
#     side equation (signed)
#       run_h_s  = +1 if it ends in a run of >= 2 Heads, -1 if >= 2 Tails,
#                  0 if there is no terminal run
#       imb_h_s  = (#Heads - #Tails) / 6
#   Additive specification (no interaction of the two cues).
#
# SAMPLE
#   synthesis-full for each treatment in cfg$run$treatment, from
#   analysis_trials() (09_analysis-sample.R): all trials of all
#   participants; the side equation uses the betting trials.
#
# POSTERIOR PREDICTIVE CHECKS
#   Replicates are generated in R for a subsample of cfg$model$ppc$synthesis_k
#   posterior draws, using the fitted participant and sequence effects, with
#   the observed betting trials fixed in the side equation. Statistics: the
#   overall betting rate, the Heads share among bets, the SD across
#   participants of the betting rate, the SD across sequences of the betting
#   rate, and the SD across sequences of the Heads share. p is two-sided,
#   2 * min(P(rep >= obs), P(rep <= obs)); a check passes if p >= 0.05.
#
# INPUT
#   model "mpl": sample/analysis-trials.csv, sample/sample-definitions.csv
#               (through analysis_trials())
#   scripts/stan/synthesis-joint.stan
#   cfg$model$stan$synthesis (sampler and prior), cfg$model$ppc$synthesis_k
#
# OUTPUT (model folder path_mod/synthesis-joint/)
#   fit-<tr>-full.rds, pid-levels-<tr>-full.rds, seq-levels-<tr>-full.rds,
#   data-<tr>-full.rds     list passed to Stan (cues and outcomes)
#   diagnostics.csv        convergence, sample sizes, PPC p-values
#   ppc-summary.csv        PPC statistics: observed, replicate quantiles, p
#   inputs.csv, run-info.txt
#
# NOTES
#   - The model folder is built as synthesis-joint.tmp and renamed when
#     complete; an existing folder is never overwritten
#   - Convergence criteria are cfg$model$convergence
#   - The side-equation slopes are signed toward Heads: a negative run slope
#     means betting against the run, a positive imbalance slope means betting
#     on the majority symbol
# ============================================================

options(mc.cores = parallel::detectCores())
rstan::rstan_options(auto_write = FALSE)

synthesis_stan_joint <- function(cfg) {

  run <- model_begin("synthesis-joint")
  if (is.null(run)) return(invisible(NULL))

  seed   <- as.integer(cfg$run$seed)
  tr_vec <- unique(as.character(cfg$run$treatment))
  design <- cfg$design

  stan_file <- here::here("scripts", "stan", "synthesis-joint.stan")
  stopifnot(file.exists(stan_file))

  lab_heads <- as.character(design$seq$side_labels$heads)
  lab_tails <- as.character(design$seq$side_labels$tails)
  stopifnot(length(lab_heads) == 1L, nzchar(lab_heads),
            length(lab_tails) == 1L, nzchar(lab_tails), lab_heads != lab_tails)

  st              <- cfg$model$stan$synthesis
  iter_val        <- as.integer(st$iter)
  warmup_val      <- as.integer(st$warmup)
  chains_val      <- as.integer(st$chains)
  adapt_delta_val <- as.numeric(st$adapt_delta)
  treedepth_val   <- as.integer(st$treedepth)

  pr <- st$prior
  pr_names <- c("beta_sd", "tau_scale", "lkj_eta", "sigma_v_scale")
  stopifnot(!is.null(pr), all(pr_names %in% names(pr)),
            all(vapply(pr[pr_names], function(x) length(x) == 1L && is.finite(x) && x > 0, logical(1))))

  K_ppc <- as.integer(cfg$model$ppc$synthesis_k)
  stopifnot(length(K_ppc) == 1L, is.finite(K_ppc), K_ppc >= 50L)

  sm <- rstan::stan_model(stan_file)

  stem <- function(prefix, tr, sample, ext = "rds") paste0(prefix, "-", tr, "-", sample, ".", ext)

  # Cues from the sequence string (outcome-blind)
  cues <- function(s) {
    ch <- strsplit(s, "", fixed = TRUE)[[1L]]
    stopifnot(length(ch) == 6L, all(ch %in% c(lab_heads, lab_tails)))
    n_h <- sum(ch == lab_heads); n_t <- 6L - n_h
    last <- ch[6L]
    run_len <- 1L
    while (run_len < 6L && ch[6L - run_len] == last) run_len <- run_len + 1L
    term <- if (run_len >= 2L) (if (last == lab_heads) 1 else -1) else 0
    list(imb_h_s = (n_h - n_t) / 6, imb_mag = abs(n_h - n_t) / 6,
         run_h_s = term, run_dum = as.numeric(term != 0))
  }

  # ---- Posterior predictive checks for one fit ----
  ppc_one <- function(fit, dl, tr, sample) {

    post  <- rstan::extract(fit, pars = c("beta_bet", "beta_side", "u", "v_bet", "v_side"))
    K_all <- nrow(post$beta_bet)
    K_use <- min(K_ppc, K_all)
    set.seed(seed + 1001L)
    k_idx <- sort(sample.int(K_all, K_use))

    n_pid <- tabulate(dl$pid, dl$P);  n_seq <- tabulate(dl$sid, dl$S)
    m_seq <- tabulate(dl$sid_s, dl$S)
    stats <- function(bet, side) {
      c(bet_rate         = mean(bet),
        heads_share      = mean(side),
        sd_pid_bet       = sd(as.numeric(rowsum(bet, dl$pid)) / n_pid),
        sd_seq_bet       = sd(as.numeric(rowsum(bet, dl$sid)) / n_seq),
        sd_seq_heads     = sd((as.numeric(rowsum(side, dl$sid_s)) / m_seq)[m_seq > 0]))
    }
    obs <- stats(dl$bet, dl$side)

    set.seed(seed + 2000L)
    rep_stats <- vapply(k_idx, function(k) {
      eta_b <- post$beta_bet[k, 1] + post$u[k, 1, dl$pid] + post$v_bet[k, dl$sid] +
        (post$beta_bet[k, 2] + post$u[k, 2, dl$pid]) * dl$run_dum +
        (post$beta_bet[k, 3] + post$u[k, 3, dl$pid]) * dl$imb_mag
      eta_s <- post$beta_side[k, 1] + post$u[k, 4, dl$pid_s] + post$v_side[k, dl$sid_s] +
        (post$beta_side[k, 2] + post$u[k, 5, dl$pid_s]) * dl$run_h_s +
        (post$beta_side[k, 3] + post$u[k, 6, dl$pid_s]) * dl$imb_h_s
      stats(rbinom(dl$N, 1L, plogis(eta_b)), rbinom(dl$N_bet, 1L, plogis(eta_s)))
    }, numeric(length(obs)))

    labels <- c(bet_rate     = "Betting rate",
                heads_share  = "Heads share among bets",
                sd_pid_bet   = "SD of betting rate across participants",
                sd_seq_bet   = "SD of betting rate across sequences",
                sd_seq_heads = "SD of Heads share across sequences")
    rbindlist(lapply(names(obs), function(s) {
      r <- rep_stats[s, ]
      p <- min(1, 2 * min(mean(r >= obs[[s]]), mean(r <= obs[[s]])))
      data.table(treatment = tr, sample = sample, statistic = s, label = labels[[s]],
                 observed = obs[[s]], rep_median = median(r),
                 rep_q025 = unname(quantile(r, 0.025)), rep_q975 = unname(quantile(r, 0.975)),
                 n_rep = K_use, p_value = p, threshold = 0.05, pass = p >= 0.05)
    }))
  }

  # ---- Fit one treatment ----
  fit_one <- function(tr, sample) {

    smp <- paste0("synthesis-", sample)
    d   <- analysis_trials(tr, smp)
    d[, `:=`(pid = as.character(pid), seq = as.character(seq))]
    d[, bet := as.integer(side %in% c(lab_heads, lab_tails) & is.finite(stake) & stake > 0)]

    bad <- d[bet == 1L & !(side %in% c(lab_heads, lab_tails))]
    stopifnot(nrow(bad) == 0L)

    pid_levels <- sort(unique(d$pid))
    seq_levels <- sort(unique(d$seq))
    cue_tbl <- rbindlist(lapply(seq_levels, function(s) as.data.table(cues(s))))
    cue_tbl[, seq := seq_levels]
    d <- merge(d, cue_tbl, by = "seq", sort = FALSE)
    d[, `:=`(pid_i = match(pid, pid_levels), sid_s = match(seq, seq_levels))]
    setorder(d, pid_i, block, pos)

    bd <- d[bet == 1L]
    data_list <- list(
      P = length(pid_levels), S = length(seq_levels), N = nrow(d),
      pid = as.integer(d$pid_i), sid = as.integer(d$sid_s),
      imb_mag = as.numeric(d$imb_mag), run_dum = as.numeric(d$run_dum),
      bet = as.integer(d$bet),
      N_bet = nrow(bd),
      pid_s = as.integer(bd$pid_i), sid_s = as.integer(bd$sid_s),
      run_h_s = as.numeric(bd$run_h_s), imb_h_s = as.numeric(bd$imb_h_s),
      side = as.integer(bd$side == lab_heads),
      prior_only = 0L,
      prior_beta_sd = as.numeric(pr$beta_sd), prior_tau_scale = as.numeric(pr$tau_scale),
      prior_lkj_eta = as.numeric(pr$lkj_eta), prior_sigma_v_scale = as.numeric(pr$sigma_v_scale)
    )

    msg("Synthesis joint: fitting ", tr, "/", sample, " | P=", data_list$P,
        " | S=", data_list$S, " | N=", data_list$N, " | N_bet=", data_list$N_bet)

    fit <- rstan::sampling(
      sm, data = data_list, iter = iter_val, warmup = warmup_val,
      chains = chains_val, seed = seed,
      control = list(adapt_delta = adapt_delta_val, max_treedepth = treedepth_val)
    )

    saveRDS(fit,        file.path(run, stem("fit",        tr, sample)))
    saveRDS(pid_levels, file.path(run, stem("pid-levels", tr, sample)))
    saveRDS(seq_levels, file.path(run, stem("seq-levels", tr, sample)))
    saveRDS(data_list,  file.path(run, stem("data",       tr, sample)))

    tag  <- paste(tr, sample, sep = "_")
    conv <- stan_diagnostics(fit, tag = tag, treatment = tr, sample = sample, cfg = cfg,
                             max_treedepth = treedepth_val)

    ppc <- ppc_one(fit, data_list, tr, sample)
    ppc_rows <- ppc[, .(tag = tag, treatment, sample,
                        check = paste0("ppc_", statistic, "_p"),
                        value = p_value, threshold, pass)]

    def  <- sample_definition(tr, smp)
    size <- data.table(
      tag = tag, treatment = tr, sample = sample,
      check = paste0(smp, c("_participants", "_trials")),
      value = c(length(pid_levels), nrow(d)),
      threshold = c(def$n_participants, def$n_trials)
    )
    size[, pass := value == threshold]

    list(diag = rbind(conv, size, ppc_rows, fill = TRUE), ppc = ppc)
  }

  results <- list()
  for (tr in tr_vec) results[[tr]] <- fit_one(tr, "full")

  fwrite(rbindlist(lapply(results, `[[`, "ppc")), file.path(run, "ppc-summary.csv"))

  # ---- Diagnostics (written last: marks the run as complete) ----
  write_diagnostics(run, rbindlist(lapply(results, `[[`, "diag"), fill = TRUE))

  model_end(run, cfg, settings = st, key = "synthesis")
}
