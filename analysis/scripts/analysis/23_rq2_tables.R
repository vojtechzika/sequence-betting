# ============================================================
# 23_rq2_tables.R
#
# PURPOSE
#   Produces sequence-level, participant-level, and model summary
#   tables for RQ2 (intensive margin of betting).
#   Reads rq2_diagnostics.csv to determine selected model.
#   Selected model tables -> path_out/
#   Other model tables    -> path_out/alternatives/
#   Sensitivity summary for confirmatory subset always produced.
#
# INPUT
#   path_src/master_sequences.csv
#   path_out/rq2_diagnostics.csv
#   path_mod/rq2_fit_sequences_<tr>_<tag>[_alt|_floor*|_mad].rds
#   path_mod/rq2_prepared_<tr>_<tag>[_alt|_floor*|_mad].rds
#
# OUTPUT
#   path_out/rq2_<tr>_<tag>_sequences.csv
#   path_out/rq2_<tr>_<tag>_sequences_summary.csv
#   path_out/rq2_<tr>_<tag>_participants.csv
#   path_out/rq2_<tr>_<tag>_participants_summary.csv
#   path_out/rq2_<tr>_<tag>_model_summary.csv
#   path_out/alternatives/alt_rq2_<tr>_<tag>_*.csv
#   path_out/rq2_<tr>_confirmatory_sensitivity_summary.csv
#
# TAGS
#   full          -- all participants
#   confirmatory  -- normative betters only
# ============================================================

rq2_tables <- function(cfg) {
  
  tr_vec   <- unique(as.character(cfg$run$treatment))
  design   <- cfg$design
  
  e <- as.numeric(design$seq$endowment)
  stopifnot(length(e) == 1L, is.finite(e), e > 0)
  
  rho_vec     <- as.numeric(design$rq2$rho)
  rho_main    <- rho_vec[1]
  rho_nm_main <- gsub("\\.", "", sprintf("%.2f", rho_main))
  
  sd_floor_sens_low  <- as.numeric(design$rq2$sd_floor_sens_low)
  sd_floor_sens_high <- as.numeric(design$rq2$sd_floor_sens_high)
  
  # ---- Read diagnostics to determine selected model ----
  f_diag <- file.path(path_out, "rq2_diagnostics.csv")
  diag   <- if (file.exists(f_diag)) fread(f_diag) else NULL
  
  get_selected_suffix <- function(tr, tg) {
    if (is.null(diag)) return("")
    row <- diag[treatment == tr & tag == tg]
    if (nrow(row) == 0L) return("")
    if (row$selected_model == "alternative") "_alt" else ""
  }
  
  # ---- Output dirs ----
  path_out_alt <- file.path(path_out, "alternatives")
  dir.create(path_out_alt, showWarnings = FALSE, recursive = TRUE)
  
  tags    <- c("full", "confirmatory")
  outputs <- list()
  
  # ---- Helper: extract mu_s and mu_i draws ----
  extract_draws <- function(fit, pid_levels, seq_levels, delta_bar_vec, sd_star_vec) {
    post  <- rstan::extract(fit)
    stopifnot(!is.null(post$alpha), !is.null(post$u), !is.null(post$b))
    
    alpha <- as.numeric(post$alpha)
    u     <- post$u
    b     <- post$b
    K     <- length(alpha)
    N     <- length(pid_levels)
    S     <- length(seq_levels)
    
    stopifnot(nrow(u) == K, ncol(u) == N, nrow(b) == K, ncol(b) == S)
    
    eta_base <- sweep(u, 1, alpha, "+")
    
    mu_s_draws <- matrix(NA_real_, nrow = K, ncol = S)
    for (s in seq_len(S)) {
      eta_mat         <- eta_base + b[, s]
      mapped          <- sweep(eta_mat, 2, sd_star_vec, "*")
      mapped          <- sweep(mapped,  2, delta_bar_vec, "+")
      mu_s_draws[, s] <- rowMeans(mapped) / e
    }
    
    mu_i_draws <- (sweep(eta_base, 2, sd_star_vec, "*") +
                     rep(delta_bar_vec, each = K)) / e
    
    list(mu_s = mu_s_draws, mu_i = mu_i_draws, K = K, N = N, S = S)
  }
  
  # ---- Helper: build sequence table ----
  build_seq_tbl <- function(mu_s_draws, seq_levels, n_trials_by_seq) {
    
    seq_tbl <- data.table(
      sequence    = seq_levels,
      n_trials    = as.integer(n_trials_by_seq[.(seq_levels), n_trials]),
      mu_a_median = apply(mu_s_draws, 2, median),
      mu_a_mean   = apply(mu_s_draws, 2, mean),
      mu_a_q025   = apply(mu_s_draws, 2, quantile, probs = 0.025),
      mu_a_q975   = apply(mu_s_draws, 2, quantile, probs = 0.975)
    )
    seq_tbl[is.na(n_trials), n_trials := 0L]
    
    for (rho in rho_vec) {
      nmU <- paste0("U_a_rho_", gsub("\\.", "", sprintf("%.2f", rho)))
      nmO <- paste0("O_a_rho_", gsub("\\.", "", sprintf("%.2f", rho)))
      seq_tbl[, (nmU) := apply(mu_s_draws, 2, function(x) mean(x < -rho))]
      seq_tbl[, (nmO) := apply(mu_s_draws, 2, function(x) mean(x >  rho))]
    }
    
    colU <- paste0("U_a_rho_", rho_nm_main)
    colO <- paste0("O_a_rho_", rho_nm_main)
    U    <- seq_tbl[[colU]]
    O    <- seq_tbl[[colO]]
    
    seq_tbl[, calib_label := {
      lab        <- rep("neutral", .N)
      pick_under <- (U >= 0.50 & O < 0.50) | (U >= 0.50 & O >= 0.50 & U >= O)
      pick_over  <- (O >= 0.50 & U < 0.50) | (U >= 0.50 & O >= 0.50 & O >  U)
      lab[pick_under & U >= 0.95] <- "under_strong"
      lab[pick_under & U >= 0.80 & U < 0.95] <- "under_moderate"
      lab[pick_under & U >= 0.50 & U < 0.80] <- "under_weak"
      lab[pick_over  & O >= 0.95] <- "over_strong"
      lab[pick_over  & O >= 0.80 & O < 0.95] <- "over_moderate"
      lab[pick_over  & O >= 0.50 & O < 0.80] <- "over_weak"
      lab
    }]
    
    setorder(seq_tbl, sequence)
    seq_tbl
  }
  
  # ---- Helper: build participant table ----
  build_pid_tbl <- function(mu_i_draws, pid_levels, n_bets_vec) {
    
    part_tbl <- data.table(
      pid         = pid_levels,
      n_bets      = n_bets_vec,
      mu_a_median = apply(mu_i_draws, 2, median),
      mu_a_mean   = apply(mu_i_draws, 2, mean),
      mu_a_q025   = apply(mu_i_draws, 2, quantile, probs = 0.025),
      mu_a_q975   = apply(mu_i_draws, 2, quantile, probs = 0.975)
    )
    
    for (rho in rho_vec) {
      nmU <- paste0("U_i_rho_", gsub("\\.", "", sprintf("%.2f", rho)))
      nmO <- paste0("O_i_rho_", gsub("\\.", "", sprintf("%.2f", rho)))
      part_tbl[, (nmU) := apply(mu_i_draws, 2, function(x) mean(x < -rho))]
      part_tbl[, (nmO) := apply(mu_i_draws, 2, function(x) mean(x >  rho))]
    }
    
    colU <- paste0("U_i_rho_", rho_nm_main)
    colO <- paste0("O_i_rho_", rho_nm_main)
    U    <- part_tbl[[colU]]
    O    <- part_tbl[[colO]]
    
    part_tbl[, calib_label := {
      lab        <- rep("neutral", .N)
      pick_under <- (U >= 0.75 & O < 0.75) | (U >= 0.75 & O >= 0.75 & U >= O)
      pick_over  <- (O >= 0.75 & U < 0.75) | (U >= 0.75 & O >= 0.75 & O >  U)
      lab[pick_under & U >= 0.95] <- "under_solid"
      lab[pick_under & U >= 0.90 & U < 0.95] <- "under_likely"
      lab[pick_under & U >= 0.75 & U < 0.90] <- "under_leaning"
      lab[pick_over  & O >= 0.95] <- "over_solid"
      lab[pick_over  & O >= 0.90 & O < 0.95] <- "over_likely"
      lab[pick_over  & O >= 0.75 & O < 0.90] <- "over_leaning"
      lab
    }]
    
    setorder(part_tbl, pid)
    part_tbl
  }
  
  # ---- Helper: load fit artifacts ----
  load_fit <- function(tr, tag, sfx = "") {
    f_fit  <- file.path(path_mod, paste0("rq2_fit_sequences_", tr, "_", tag, sfx, ".rds"))
    f_pid  <- file.path(path_mod, paste0("rq2_pid_levels_",    tr, "_", tag, sfx, ".rds"))
    f_seq  <- file.path(path_mod, paste0("rq2_seq_levels_",    tr, "_", tag, sfx, ".rds"))
    f_prep <- file.path(path_mod, paste0("rq2_prepared_",      tr, "_", tag, sfx, ".rds"))
    if (!all(file.exists(c(f_fit, f_pid, f_seq, f_prep)))) return(NULL)
    list(fit  = readRDS(f_fit),
         pids = as.character(readRDS(f_pid)),
         seqs = as.character(readRDS(f_seq)),
         prep = readRDS(f_prep))
  }
  
  # ---- Helper: get pid constants from prep ----
  get_pid_constants <- function(prep, pid_levels) {
    prep[, pid := as.character(pid)]
    pid_map <- prep[, .(delta_bar = delta_bar[1], sd_star = sd_star[1],
                        n_bets = n_bets[1]), by = pid]
    setkey(pid_map, pid)
    pid_map <- pid_map[.(pid_levels)]
    stopifnot(!anyNA(pid_map$pid))
    pid_map
  }
  
  # ---- Helper: produce tables for one fit ----
  produce_tables <- function(tr, tag, suffix, outdir, is_alt = FALSE) {
    
    arts <- load_fit(tr, tag, suffix)
    if (is.null(arts)) return(invisible(NULL))
    
    pid_levels <- arts$pids
    seq_levels <- arts$seqs
    prep       <- arts$prep
    
    pid_map       <- get_pid_constants(prep, pid_levels)
    delta_bar_vec <- as.numeric(pid_map$delta_bar)
    sd_star_vec   <- as.numeric(pid_map$sd_star)
    n_bets_vec    <- as.integer(pid_map$n_bets)
    
    n_trials_by_seq <- prep[pid %in% pid_levels, .(n_trials = .N), by = seq]
    setkey(n_trials_by_seq, seq)
    
    draws <- extract_draws(arts$fit, pid_levels, seq_levels,
                           delta_bar_vec, sd_star_vec)
    
    file_prefix <- if (is_alt) "alt_" else ""
    file_stem   <- paste0(file_prefix, "rq2_", tr, "_", tag)
    
    sum_draw <- function(x, nm) {
      data.table(parameter = nm, median = median(x), mean = mean(x),
                 q025 = quantile(x, 0.025), q975 = quantile(x, 0.975))
    }
    
    # ---- Sequence table ----
    f_seq_csv <- file.path(outdir, paste0(file_stem, "_sequences.csv"))
    if (!should_skip(f_seq_csv, cfg, "output",
                     paste0("RQ2 sequences (", tr, "/", tag,
                            if (is_alt) "/alt" else "", ")"))) {
      seq_tbl <- build_seq_tbl(draws$mu_s, seq_levels, n_trials_by_seq)
      fwrite(seq_tbl, f_seq_csv)
      msg("Saved: ", f_seq_csv)
    } else {
      seq_tbl <- fread(f_seq_csv)
    }
    
    # Sequence summary
    f_seq_sum <- file.path(outdir, paste0(file_stem, "_sequences_summary.csv"))
    if (!should_skip(f_seq_sum, cfg, "output",
                     paste0("RQ2 sequences summary (", tr, "/", tag,
                            if (is_alt) "/alt" else "", ")"))) {
      
      seq_summary <- seq_tbl[, .N, by = calib_label]
      seq_summary[, calib_label := factor(calib_label,
                                          levels = c("under_strong", "under_moderate", "under_weak",
                                                     "neutral",
                                                     "over_weak", "over_moderate", "over_strong"))]
      setorder(seq_summary, calib_label)
      seq_summary[, pct := round(100 * N / nrow(seq_tbl), 1)]
      
      seq_stats <- seq_tbl[, .(
        mu_mean   = round(mean(mu_a_mean),   3),
        mu_median = round(median(mu_a_mean), 3),
        mu_min    = round(min(mu_a_mean),    3),
        mu_max    = round(max(mu_a_mean),    3)
      ), by = calib_label]
      
      seq_ids <- seq_tbl[, .(
        sequences = paste(sort(sequence), collapse = ", ")
      ), by = calib_label]
      
      seq_summary <- merge(seq_summary, seq_stats, by = "calib_label")
      seq_summary <- merge(seq_summary, seq_ids,   by = "calib_label")
      fwrite(seq_summary, f_seq_sum)
      msg("Saved: ", f_seq_sum)
    }
    
    # ---- Participant table ----
    f_pid_csv <- file.path(outdir, paste0(file_stem, "_participants.csv"))
    if (!should_skip(f_pid_csv, cfg, "output",
                     paste0("RQ2 participants (", tr, "/", tag,
                            if (is_alt) "/alt" else "", ")"))) {
      part_tbl <- build_pid_tbl(draws$mu_i, pid_levels, n_bets_vec)
      fwrite(part_tbl, f_pid_csv)
      msg("Saved: ", f_pid_csv)
    } else {
      part_tbl <- fread(f_pid_csv)
    }
    
    # Participant summary
    f_pid_sum <- file.path(outdir, paste0(file_stem, "_participants_summary.csv"))
    if (!should_skip(f_pid_sum, cfg, "output",
                     paste0("RQ2 participants summary (", tr, "/", tag,
                            if (is_alt) "/alt" else "", ")"))) {
      
      pid_summary <- part_tbl[, .N, by = calib_label]
      pid_summary[, calib_label := factor(calib_label,
                                          levels = c("under_solid", "under_likely", "under_leaning",
                                                     "neutral",
                                                     "over_leaning", "over_likely", "over_solid"))]
      setorder(pid_summary, calib_label)
      pid_summary[, pct := round(100 * N / nrow(part_tbl), 1)]
      
      pid_stats <- part_tbl[, .(
        mu_mean   = round(mean(mu_a_mean),   3),
        mu_median = round(median(mu_a_mean), 3),
        mu_min    = round(min(mu_a_mean),    3),
        mu_max    = round(max(mu_a_mean),    3)
      ), by = calib_label]
      
      pid_ids <- part_tbl[, .(
        pids = paste(sort(pid), collapse = ", ")
      ), by = calib_label]
      
      pid_summary <- merge(pid_summary, pid_stats, by = "calib_label")
      pid_summary <- merge(pid_summary, pid_ids,   by = "calib_label")
      fwrite(pid_summary, f_pid_sum)
      msg("Saved: ", f_pid_sum)
    }
    
    # ---- Model summary ----
    # ---- Model summary ----
    f_mod_csv <- file.path(outdir, paste0(file_stem, "_model_summary.csv"))
    if (!should_skip(f_mod_csv, cfg, "output",
                     paste0("RQ2 model summary (", tr, "/", tag,
                            if (is_alt) "/alt" else "", ")"))) {
      
      post <- rstan::extract(arts$fit)
      
      safe_draw <- function(x, nm) {
        if (is.null(x)) return(NULL)
        data.table(parameter = nm, median = median(x), mean = mean(x),
                   q025 = quantile(x, 0.025), q975 = quantile(x, 0.975))
      }
      
      if (!is_alt) {
        rows <- list(
          safe_draw(post$alpha,       "alpha (grand mean z-deviation)"),
          safe_draw(post$sigma_u,     "sigma_u (between-participant SD)"),
          safe_draw(post$sigma_s,     "sigma_s (between-sequence SD)"),
          safe_draw(post$sigma,       "sigma (residual SD)"),
          safe_draw(post$gamma_drift, "gamma_drift (linear drift per block)")
        )
      } else {
        rows <- list(
          safe_draw(post$alpha,       "alpha (grand mean)"),
          safe_draw(post$sigma_u,     "sigma_u (between-participant SD)"),
          safe_draw(post$sigma_s,     "sigma_s (between-sequence SD)"),
          safe_draw(post$phi,         "phi (BB overdispersion)"),
          safe_draw(post$alpha_omega, "alpha_omega (all-in intercept)"),
          safe_draw(post$sigma_v,     "sigma_v (all-in between-participant SD)"),
          safe_draw(post$gamma_drift, "gamma_drift (linear drift per block)")
        )
      }
      
      mod_tbl <- rbindlist(Filter(Negate(is.null), rows), fill = TRUE)
      mod_tbl[, treatment  := tr]
      mod_tbl[, tag        := tag]
      mod_tbl[, likelihood := if (is_alt) "beta_binomial" else "gaussian"]
      
      grand_mu_a <- data.table(
        parameter  = "grand mean mu_a (absolute scale, proportion of endowment)",
        median     = median(apply(draws$mu_s, 1, mean)),
        mean       = mean(apply(draws$mu_s, 1, mean)),
        q025       = quantile(apply(draws$mu_s, 1, mean), 0.025),
        q975       = quantile(apply(draws$mu_s, 1, mean), 0.975),
        treatment  = tr,
        tag        = tag,
        likelihood = if (is_alt) "beta_binomial" else "gaussian"
      )
      
      mod_tbl <- rbindlist(list(mod_tbl, grand_mu_a), fill = TRUE)
      fwrite(mod_tbl, f_mod_csv)
      msg("Saved: ", f_mod_csv)
    }
    
    list(draws                    = draws,
         seq_tbl                  = seq_tbl,
         sequences_csv            = f_seq_csv,
         sequences_summary_csv    = f_seq_sum,
         participants_csv         = f_pid_csv,
         participants_summary_csv = f_pid_sum,
         model_summary_csv        = f_mod_csv)
  }
  
  # ============================================================
  # Main loop
  # ============================================================
  for (tr in tr_vec) {
    for (tag in tags) {
      
      if (tag == "confirmatory" && !isTRUE(design$a_flags$betting_normative[[tr]])) next
      
      selected_suffix <- get_selected_suffix(tr, tag)
      other_suffix    <- if (selected_suffix == "") "_alt" else ""
      
      # Selected model -> path_out
      res <- produce_tables(tr, tag,
                            suffix = selected_suffix,
                            outdir = path_out,
                            is_alt = (selected_suffix == "_alt"))
      
      # Other model -> path_out/alternatives
      produce_tables(tr, tag,
                     suffix = other_suffix,
                     outdir = path_out_alt,
                     is_alt = (other_suffix == "_alt"))
      
      # Sensitivity summary (confirmatory only, from selected model)
      if (tag == "confirmatory" && !is.null(res)) {
        
        f_sens <- file.path(path_out,
                            paste0("rq2_", tr, "_confirmatory_sensitivity_summary.csv"))
        if (!should_skip(f_sens, cfg, "output",
                         paste0("RQ2 sensitivity summary (", tr, ")"))) {
          
          sens_variants <- list(
            list(suffix = "",                                    label = "primary"),
            list(suffix = paste0("_floor", sd_floor_sens_low),  label = paste0("floor", sd_floor_sens_low)),
            list(suffix = paste0("_floor", sd_floor_sens_high), label = paste0("floor", sd_floor_sens_high)),
            list(suffix = "_mad",                                label = "mad")
          )
          
          mu_s_primary <- apply(res$draws$mu_s, 2, mean)
          mu_i_primary <- apply(res$draws$mu_i, 2, mean)
          sens_rows    <- list()
          
          for (sv in sens_variants) {
            arts_sv <- load_fit(tr, tag, sv$suffix)
            if (is.null(arts_sv)) next
            
            pid_map_sv <- get_pid_constants(arts_sv$prep, arts_sv$pids)
            draws_sv   <- extract_draws(arts_sv$fit, arts_sv$pids, arts_sv$seqs,
                                        as.numeric(pid_map_sv$delta_bar),
                                        as.numeric(pid_map_sv$sd_star))
            
            mu_s_sv <- apply(draws_sv$mu_s, 2, mean)
            mu_i_sv <- apply(draws_sv$mu_i, 2, mean)
            
            rho_seq <- cor(mu_s_primary, mu_s_sv, method = "spearman")
            rho_pid <- cor(mu_i_primary, mu_i_sv, method = "spearman")
            
            n_trials_sv <- arts_sv$prep[, .(n_trials = .N), by = seq]
            setkey(n_trials_sv, seq)
            seq_tbl_sv <- build_seq_tbl(draws_sv$mu_s, arts_sv$seqs, n_trials_sv)
            
            pct_stable <- mean(res$seq_tbl$calib_label == seq_tbl_sv$calib_label,
                               na.rm = TRUE)
            
            sens_rows[[sv$label]] <- data.table(
              treatment        = tr,
              variant          = sv$label,
              spearman_seq     = round(rho_seq, 4),
              spearman_pid     = round(rho_pid, 4),
              pct_stable_label = round(pct_stable, 4)
            )
          }
          
          if (length(sens_rows) > 0) {
            fwrite(rbindlist(sens_rows), f_sens)
            msg("Saved: ", f_sens)
          }
        }
      }
      
      if (!is.null(res)) {
        outputs[[paste(tr, tag, sep = "_")]] <- list(
          sequences_csv            = res$sequences_csv,
          sequences_summary_csv    = res$sequences_summary_csv,
          participants_csv         = res$participants_csv,
          participants_summary_csv = res$participants_summary_csv,
          model_summary_csv        = res$model_summary_csv
        )
      }
    }
  }
  
  invisible(outputs)
}