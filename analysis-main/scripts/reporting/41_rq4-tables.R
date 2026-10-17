# ============================================================
# 41_rq4-tables.R
#
# PURPOSE
#   Sequence-level, participant-level and model summary tables for RQ4
#   (side choice conditional on betting), from the Bernoulli model selected
#   under the key "rq4" in selected-models.csv. Labels and thresholds
#   (classification of sequences and participants by the posterior
#   probability of a Heads or Tails bias beyond the tolerance rho) are
#   defined in this script, not in the model folder.
#
# ESTIMANDS (generated quantities of rq4-bernoulli.stan)
#   mu_h[s]    population mean Heads probability of sequence s,
#              E_u[logistic(alpha + u + beta_s)], u ~ N(0, sigma_u)
#   hbar       population baseline, E_u[logistic(alpha + u)]
#   mu_h_i[i]  mean Heads probability of participant i, averaged over
#              sequences
#   H_rho = P(mu > hbar + rho) (Head-bias) and T_rho = P(mu < hbar - rho)
#   (Tail-bias), where rho is cfg$design$rq4$rho (first value = main, the
#   others sensitivity). P_absdev_rho = P(|mu - 0.5| > rho) is reported in
#   addition (deviation from parity).
#
# CLASSIFICATION
#   Sequences:     strong (H or T >= 0.95), moderate (>= 0.80), weak (>= 0.50)
#                  Head- or Tail-bias; neutral if both are below 0.50
#   Participants:  solid (>= 0.95), likely (>= 0.90), leaning (>= 0.75)
#                  Head- or Tail-bias; neutral if both are below 0.75
#
# SAMPLE
#   full           all participants of each treatment in cfg$run$treatment
#                  (betting trials; RQ4 has no confirmatory sample)
#
# INPUT
#   model "rq4": fit-<tr>-full.rds, pid-levels-<tr>-full.rds,
#                seq-levels-<tr>-full.rds
#   model "mpl": sample/analysis-trials.csv, sample/sample-definitions.csv
#                (through analysis_trials(); betting trials) and
#                drift/drift-decisions.rds
#
# OUTPUT (data/output/csv/)
#   rq4-<tr>-full-sequences.csv             one row per sequence
#   rq4-<tr>-full-sequences-summary.csv     sequences by label
#   rq4-<tr>-full-participants.csv          one row per participant
#   rq4-<tr>-full-participants-summary.csv  participants by label
#   rq4-<tr>-full-model-summary.csv         posterior summary of the
#                                           scalar model parameters
#
# NOTES
#   - The model summary lists every scalar parameter of the fit under its
#     Stan name; the drift parameters are listed only when drift was
#     active (gamma_drift for linear, gamma_block_1..4 for categorical)
#   - The estimands are checked against the levels files (dimensions) and
#     the probability scale (values in [0, 1]); the script stops otherwise
#   - Output files are overwritten at every run
# ============================================================

rq4_tables <- function(cfg) {

  tr_vec <- unique(as.character(cfg$run$treatment))
  design <- cfg$design
  sample <- "full"

  rho_vec  <- as.numeric(design$rq4$rho)
  rho_main <- rho_vec[1]
  stopifnot(length(rho_vec) >= 1L, all(is.finite(rho_vec)), all(rho_vec > 0))
  rho_nm   <- function(rho) gsub("\\.", "", sprintf("%.2f", rho))

  lab_heads <- as.character(design$seq$side_labels$heads)
  lab_tails <- as.character(design$seq$side_labels$tails)
  seq_n     <- as.integer(design$seq$seq_n)

  drift_decisions <- readRDS(model_path("mpl", file.path("drift", "drift-decisions.rds")))
  outcome_key     <- as.character(design$drift$rq4$outcome)

  summarise_draws <- function(x, nm) {
    data.table(parameter = nm, median = median(x), mean = mean(x),
               q025 = unname(quantile(x, 0.025)), q975 = unname(quantile(x, 0.975)))
  }

  # direction labels from the Head- and Tail-bias probabilities (vectorised)
  label_dir <- function(H, T, cuts, names, floor_p) {
    m   <- pmax(H, T)
    str <- fifelse(m >= cuts[1], names[1], fifelse(m >= cuts[2], names[2], names[3]))
    dir <- fifelse(H >= T, "head", "tail")
    fifelse(m < floor_p, "neutral", paste0(str, "_", dir))
  }

  produce_tables <- function(tr) {

    stem <- function(suffix) paste0("rq4-", tr, "-", sample, "-", suffix, ".csv")
    fit_file <- function(prefix)
      model_path("rq4", paste0(prefix, "-", tr, "-", sample, ".rds"))

    fit        <- readRDS(fit_file("fit"))
    pid_levels <- as.character(readRDS(fit_file("pid-levels")))
    seq_levels <- as.character(readRDS(fit_file("seq-levels")))
    stopifnot(length(seq_levels) == seq_n)

    d <- analysis_trials(tr, paste0("rq4-", sample))
    stopifnot(setequal(unique(d$pid), pid_levels),
              all(d$side %in% c(lab_heads, lab_tails)))

    drift <- drift_decisions[[paste(tr, sample, outcome_key, sep = "_")]]
    drift_type   <- if (!is.null(drift) && isTRUE(drift$drift)) as.character(drift$drift_type) else "none"

    # ---- Estimands and their checks ----
    post <- rstan::extract(fit)
    for (nm in c("mu_h", "mu_h_i", "hbar"))
      if (is.null(post[[nm]])) stop("RQ4 tables: '", nm, "' missing in the fit of ",
                                    tr, "/", sample, call. = FALSE)
    mu_s <- post$mu_h
    mu_i <- post$mu_h_i
    hbar <- as.numeric(post$hbar)
    stopifnot(is.matrix(mu_s), ncol(mu_s) == length(seq_levels),
              is.matrix(mu_i), ncol(mu_i) == length(pid_levels),
              all(mu_s >= 0 & mu_s <= 1), all(mu_i >= 0 & mu_i <= 1),
              all(hbar >= 0 & hbar <= 1))

    hbar_cols <- function(n) data.table(
      hbar_median = rep(median(hbar), n), hbar_mean = rep(mean(hbar), n),
      hbar_q025 = rep(unname(quantile(hbar, 0.025)), n),
      hbar_q975 = rep(unname(quantile(hbar, 0.975)), n))

    add_bias <- function(tbl, mu) {
      for (rho in rho_vec) {
        nm <- rho_nm(rho)
        tbl[, (paste0("H_rho_", nm))        := apply(mu, 2, function(x) mean(x > hbar + rho))]
        tbl[, (paste0("T_rho_", nm))        := apply(mu, 2, function(x) mean(x < hbar - rho))]
        tbl[, (paste0("P_absdev_rho_", nm)) := apply(mu, 2, function(x) mean(abs(x - 0.5) > rho))]
      }
      tbl
    }

    # ---- Sequence table ----
    n_seq <- d[, .(n_trials = .N), by = seq]
    setkey(n_seq, seq)
    seq_tbl <- cbind(
      data.table(sequence = seq_levels,
                 n_trials = as.integer(n_seq[.(seq_levels), n_trials]),
                 mu_h_median = apply(mu_s, 2, median), mu_h_mean = apply(mu_s, 2, mean),
                 mu_h_q025 = apply(mu_s, 2, quantile, probs = 0.025),
                 mu_h_q975 = apply(mu_s, 2, quantile, probs = 0.975)),
      hbar_cols(length(seq_levels)))
    seq_tbl <- add_bias(seq_tbl, mu_s)
    nm_main <- rho_nm(rho_main)
    seq_tbl[, direction_label := label_dir(get(paste0("H_rho_", nm_main)),
                                           get(paste0("T_rho_", nm_main)),
                                           c(0.95, 0.80), c("strong", "moderate", "weak"), 0.50)]
    setorder(seq_tbl, sequence)
    fwrite(seq_tbl, file.path(path_csv, stem("sequences")))
    msg("Saved: ", file.path(path_csv, stem("sequences")))

    # ---- Sequence summary ----
    dir_levs <- c("strong_head", "moderate_head", "weak_head", "neutral",
                  "weak_tail", "moderate_tail", "strong_tail")
    seq_sum <- merge(
      seq_tbl[, .N, by = direction_label],
      seq_tbl[, .(mu_mean = round(mean(mu_h_mean), 3), mu_median = round(median(mu_h_mean), 3),
                  mu_min = round(min(mu_h_mean), 3), mu_max = round(max(mu_h_mean), 3),
                  sequences = paste(sort(sequence), collapse = ", ")), by = direction_label],
      by = "direction_label")
    seq_sum[, direction_label := factor(direction_label, levels = dir_levs)]
    setorder(seq_sum, direction_label)
    seq_sum[, pct := round(100 * N / nrow(seq_tbl), 1)]
    fwrite(seq_sum, file.path(path_csv, stem("sequences-summary")))
    msg("Saved: ", file.path(path_csv, stem("sequences-summary")))

    # ---- Participant table ----
    n_pid <- d[, .(n_trials = .N), by = pid]
    setkey(n_pid, pid)
    pid_tbl <- cbind(
      data.table(pid = pid_levels,
                 n_trials = as.integer(n_pid[.(pid_levels), n_trials]),
                 mu_h_median = apply(mu_i, 2, median), mu_h_mean = apply(mu_i, 2, mean),
                 mu_h_q025 = apply(mu_i, 2, quantile, probs = 0.025),
                 mu_h_q975 = apply(mu_i, 2, quantile, probs = 0.975)),
      hbar_cols(length(pid_levels)))
    pid_tbl <- add_bias(pid_tbl, mu_i)
    pid_tbl[, side_label := label_dir(get(paste0("H_rho_", nm_main)),
                                      get(paste0("T_rho_", nm_main)),
                                      c(0.95, 0.90), c("solid", "likely", "leaning"), 0.75)]
    setorder(pid_tbl, pid)
    fwrite(pid_tbl, file.path(path_csv, stem("participants")))
    msg("Saved: ", file.path(path_csv, stem("participants")))

    # ---- Participant summary ----
    side_levs <- c("solid_head", "likely_head", "leaning_head", "neutral",
                   "leaning_tail", "likely_tail", "solid_tail")
    pid_sum <- merge(
      pid_tbl[, .N, by = side_label],
      pid_tbl[, .(mu_mean = round(mean(mu_h_mean), 3), mu_median = round(median(mu_h_mean), 3),
                  mu_min = round(min(mu_h_mean), 3), mu_max = round(max(mu_h_mean), 3),
                  pids = paste(sort(pid), collapse = ", ")), by = side_label],
      by = "side_label")
    pid_sum[, side_label := factor(side_label, levels = side_levs)]
    setorder(pid_sum, side_label)
    pid_sum[, pct := round(100 * N / nrow(pid_tbl), 1)]
    fwrite(pid_sum, file.path(path_csv, stem("participants-summary")))
    msg("Saved: ", file.path(path_csv, stem("participants-summary")))

    # ---- Model summary: scalar parameters of the fit ----
    scalar <- setdiff(names(post)[vapply(post, function(x) length(dim(x)) <= 1L, logical(1))], "lp__")
    scalar <- setdiff(scalar, "gamma_drift")
    rows <- lapply(scalar, function(nm) summarise_draws(as.numeric(post[[nm]]), nm))
    if (drift_type == "linear")
      rows <- c(rows, list(summarise_draws(as.numeric(post$gamma_drift), "gamma_drift")))
    if (drift_type == "categorical")
      rows <- c(rows, lapply(1:4, function(k)
        summarise_draws(as.numeric(post$gamma_block[, k]), paste0("gamma_block_", k))))
    mod_tbl <- rbindlist(rows)
    # (values are set outside data.table's scope, where column names such as
    # "sample" and "drift" would mask the objects of the same name)
    n_tr <- nrow(d); n_pid_total <- length(pid_levels); n_seq_total <- length(seq_levels)
    mod_tbl[, `:=`(treatment = tr, sample = sample, likelihood = "bernoulli",
                   drift = drift_type, n_trials = n_tr,
                   n_participants = n_pid_total, n_sequences = n_seq_total)]
    fwrite(mod_tbl, file.path(path_csv, stem("model-summary")))
    msg("Saved: ", file.path(path_csv, stem("model-summary")))

    invisible(NULL)
  }

  for (tr in tr_vec) produce_tables(tr)

  invisible(TRUE)
}
