# ============================================================
# 11_rq1-tables.R
#
# PURPOSE
#   Sequence-level, participant-level and model summary tables for RQ1
#   (probability of betting), from the Bernoulli model selected under the
#   key "rq1" in selected-models.csv. Labels and thresholds (the
#   classification of sequences and participants by the posterior
#   probability of under-betting) are defined in this script, not in the
#   model folder.
#
# ESTIMANDS (generated quantities of rq1-bernoulli.stan)
#   mu_b[s]          mean betting probability of sequence s, averaged over
#                    the participants of the sample
#   mu_b_i[i]        mean betting probability of participant i, averaged
#                    over sequences
#   grand_mean_bet   mean betting probability over participants and sequences
#   Under-betting of a sequence or participant is P(mu < 1 - rho), where
#   rho is cfg$design$rq1$rho (first value = main, the others sensitivity).
#
# CLASSIFICATION
#   Sequences:     strong (P >= 0.95), moderate (>= 0.80), weak (>= 0.50),
#                  neutral
#   Participants:  solid (P >= 0.95), likely (>= 0.90), leaning (>= 0.75),
#                  neutral
#   Grand-mean classes (sequences): above / likely_above (P(mu_s > grand
#   mean) >= 0.95 / 0.80), below / likely_below, neutral.
#
# SAMPLES
#   full           all participants of each treatment in cfg$run$treatment
#   confirmatory   treatments with cfg$design$a_flags$betting_normative =
#                  TRUE: participants kept by the confirmatory flags
#
# INPUT
#   model "rq1": fit-<tr>-<sample>.rds, pid-levels-<tr>-<sample>.rds,
#                seq-levels-<tr>-<sample>.rds
#   model "mpl": sample/analysis-trials.csv, sample/sample-definitions.csv
#                (through analysis_trials(); trial counts), drift/drift-decisions.rds
#
# OUTPUT (data/output/csv/)
#   rq1-<tr>-<sample>-sequences.csv             one row per sequence
#   rq1-<tr>-<sample>-sequences-summary.csv     sequences by label
#   rq1-<tr>-<sample>-participants.csv          one row per participant
#   rq1-<tr>-<sample>-participants-summary.csv  participants by label
#   rq1-<tr>-<sample>-model-summary.csv         posterior summary of the
#                                               scalar model parameters
#
# NOTES
#   - The model summary lists every scalar parameter of the fit under its
#     Stan name (see rq1-bernoulli.stan); gamma_drift is listed only when
#     linear drift was active for the treatment and sample
#   - The estimands are checked against the levels files (dimensions) and
#     the probability scale (values in [0, 1]); the script stops otherwise
#   - Output files are overwritten at every run
# ============================================================

rq1_tables <- function(cfg) {
  
  tr_vec <- unique(as.character(cfg$run$treatment))
  design <- cfg$design
  
  rho_vec  <- as.numeric(design$rq1$rho)
  rho_main <- rho_vec[1]
  rho_nm   <- function(rho) gsub("\\.", "", sprintf("%.2f", rho))
  seq_n    <- as.integer(design$seq$seq_n)
  
  drift_decisions <- readRDS(model_path("mpl", file.path("drift", "drift-decisions.rds")))
  
  summarise_draws <- function(x, nm) {
    data.table(parameter = nm, median = median(x), mean = mean(x),
               q025 = unname(quantile(x, 0.025)), q975 = unname(quantile(x, 0.975)))
  }
  
  # ---- Tables for one treatment and sample ----
  produce_tables <- function(tr, sample) {
    
    stem <- function(suffix) paste0("rq1-", tr, "-", sample, "-", suffix, ".csv")
    fit_file <- function(prefix)
      model_path("rq1", paste0(prefix, "-", tr, "-", sample, ".rds"))
    
    fit        <- readRDS(fit_file("fit"))
    pid_levels <- as.character(readRDS(fit_file("pid-levels")))
    seq_levels <- as.character(readRDS(fit_file("seq-levels")))
    stopifnot(length(seq_levels) == seq_n)
    
    # trials of the sample, for the trial counts
    d <- analysis_trials(tr, paste0("rq1-", sample))
    stopifnot(setequal(unique(d$pid), pid_levels))
    
    # drift decision used by the model
    drift <- drift_decisions[[paste(tr, sample, "bet", sep = "_")]]
    drift_active <- !is.null(drift) && isTRUE(drift$drift)
    
    # ---- Estimands and their checks ----
    post <- rstan::extract(fit)
    for (nm in c("mu_b", "mu_b_i", "grand_mean_bet"))
      if (is.null(post[[nm]])) stop("RQ1 tables: '", nm, "' missing in the fit of ",
                                    tr, "/", sample, call. = FALSE)
    mu_s  <- post$mu_b
    mu_i  <- post$mu_b_i
    grand <- as.numeric(post$grand_mean_bet)
    stopifnot(is.matrix(mu_s), ncol(mu_s) == length(seq_levels),
              is.matrix(mu_i), ncol(mu_i) == length(pid_levels),
              all(mu_s >= 0 & mu_s <= 1), all(mu_i >= 0 & mu_i <= 1))
    
    # ---- Sequence table ----
    n_seq <- d[, .(n_trials = .N), by = seq]
    setkey(n_seq, seq)
    
    seq_tbl <- data.table(
      sequence  = seq_levels,
      n_trials  = as.integer(n_seq[.(seq_levels), n_trials]),
      mu_b_median = apply(mu_s, 2, median),
      mu_b_mean   = apply(mu_s, 2, mean),
      mu_b_q025   = apply(mu_s, 2, quantile, probs = 0.025),
      mu_b_q975   = apply(mu_s, 2, quantile, probs = 0.975)
    )
    for (rho in rho_vec)
      seq_tbl[, (paste0("U_b_rho_", rho_nm(rho))) := apply(mu_s, 2, function(x) mean(x < 1 - rho))]
    
    U <- seq_tbl[[paste0("U_b_rho_", rho_nm(rho_main))]]
    seq_tbl[, underbet_label := fifelse(U >= 0.95, "strong",
                                fifelse(U >= 0.80, "moderate",
                                fifelse(U >= 0.50, "weak", "neutral")))]
    
    seq_tbl[, p_above_grand := apply(mu_s, 2, function(x) mean(x > grand))]
    seq_tbl[, p_below_grand := apply(mu_s, 2, function(x) mean(x < grand))]
    seq_tbl[, grand_label := fifelse(p_above_grand >= 0.95, "above",
                             fifelse(p_above_grand >= 0.80, "likely_above",
                             fifelse(p_below_grand >= 0.95, "below",
                             fifelse(p_below_grand >= 0.80, "likely_below", "neutral"))))]
    setorder(seq_tbl, sequence)
    fwrite(seq_tbl, file.path(path_csv, stem("sequences")))
    msg("Saved: ", file.path(path_csv, stem("sequences")))
    
    # ---- Sequence summary: counts, mean estimates and sequence ids by label ----
    summarise_by <- function(lab, levs, classification) {
      cnt <- seq_tbl[, .N, by = c(lab)]
      setnames(cnt, lab, "label")
      st <- seq_tbl[, .(mu_mean   = round(mean(mu_b_mean),   3),
                        mu_median = round(median(mu_b_mean), 3),
                        mu_min    = round(min(mu_b_mean),    3),
                        mu_max    = round(max(mu_b_mean),    3),
                        sequences = paste(sort(sequence), collapse = ", ")),
                    by = c(lab)]
      setnames(st, lab, "label")
      out <- merge(cnt, st, by = "label")
      out[, label := factor(label, levels = levs)]
      setorder(out, label)
      out[, pct := round(100 * N / nrow(seq_tbl), 1)]
      out[, classification := classification]
      out
    }
    seq_sum <- rbindlist(list(
      summarise_by("underbet_label", c("strong", "moderate", "weak", "neutral"), "preregistered"),
      summarise_by("grand_label", c("above", "likely_above", "neutral", "likely_below", "below"),
                   "grand_mean")
    ))
    fwrite(seq_sum, file.path(path_csv, stem("sequences-summary")))
    msg("Saved: ", file.path(path_csv, stem("sequences-summary")))
    
    # ---- Participant table ----
    n_pid <- d[, .(n_trials = .N), by = pid]
    setkey(n_pid, pid)
    
    pid_tbl <- data.table(
      pid         = pid_levels,
      n_trials    = as.integer(n_pid[.(pid_levels), n_trials]),
      mu_i_median = apply(mu_i, 2, median),
      mu_i_mean   = apply(mu_i, 2, mean),
      mu_i_q025   = apply(mu_i, 2, quantile, probs = 0.025),
      mu_i_q975   = apply(mu_i, 2, quantile, probs = 0.975)
    )
    for (rho in rho_vec)
      pid_tbl[, (paste0("U_i_rho_", rho_nm(rho))) := apply(mu_i, 2, function(x) mean(x < 1 - rho))]
    
    Ui <- pid_tbl[[paste0("U_i_rho_", rho_nm(rho_main))]]
    pid_tbl[, underbet_label := fifelse(Ui >= 0.95, "solid",
                                fifelse(Ui >= 0.90, "likely",
                                fifelse(Ui >= 0.75, "leaning", "neutral")))]
    setorder(pid_tbl, pid)
    fwrite(pid_tbl, file.path(path_csv, stem("participants")))
    msg("Saved: ", file.path(path_csv, stem("participants")))
    
    # ---- Participant summary ----
    pid_levs <- c("solid", "likely", "leaning", "neutral")
    pid_sum <- merge(
      pid_tbl[, .N, by = underbet_label],
      pid_tbl[, .(mu_mean   = round(mean(mu_i_mean),   3),
                  mu_median = round(median(mu_i_mean), 3),
                  mu_min    = round(min(mu_i_mean),    3),
                  mu_max    = round(max(mu_i_mean),    3),
                  pids      = paste(sort(pid), collapse = ", ")),
              by = underbet_label],
      by = "underbet_label")
    pid_sum[, underbet_label := factor(underbet_label, levels = pid_levs)]
    setorder(pid_sum, underbet_label)
    pid_sum[, pct := round(100 * N / nrow(pid_tbl), 1)]
    fwrite(pid_sum, file.path(path_csv, stem("participants-summary")))
    msg("Saved: ", file.path(path_csv, stem("participants-summary")))
    
    # ---- Model summary: scalar parameters of the fit ----
    scalar <- names(post)[vapply(post, function(x) length(dim(x)) <= 1L, logical(1))]
    scalar <- setdiff(scalar, "lp__")
    if (!drift_active) scalar <- setdiff(scalar, "gamma_drift")
    
    mod_tbl <- rbindlist(lapply(scalar, function(nm) summarise_draws(as.numeric(post[[nm]]), nm)))
    # (values are set outside data.table's scope, where column names such as
    # "sample" and "drift" would mask the objects of the same name)
    drift_lab <- if (drift_active) drift$drift_type else "none"
    n_tr <- nrow(d); n_pid_total <- length(pid_levels); n_seq_total <- length(seq_levels)
    mod_tbl[, `:=`(treatment = tr, sample = sample, likelihood = "bernoulli",
                   drift = drift_lab, n_trials = n_tr,
                   n_participants = n_pid_total, n_sequences = n_seq_total)]
    fwrite(mod_tbl, file.path(path_csv, stem("model-summary")))
    msg("Saved: ", file.path(path_csv, stem("model-summary")))
    
    invisible(NULL)
  }
  
  for (tr in tr_vec) produce_tables(tr, "full")
  for (tr in tr_vec) {
    if (!isTRUE(design$a_flags$betting_normative[[tr]])) next
    produce_tables(tr, "confirmatory")
  }
  
  invisible(TRUE)
}
