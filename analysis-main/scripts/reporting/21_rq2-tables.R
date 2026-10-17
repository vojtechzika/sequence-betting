# ============================================================
# 21_rq2-tables.R
#
# PURPOSE
#   Sequence-level, participant-level and model summary tables for RQ2
#   (stake calibration relative to the expected-utility-optimal stake),
#   from the model selected under the key "rq2" in selected-models.csv
#   (Gaussian or one-inflated Beta-Binomial). Labels and thresholds
#   (classification of sequences and participants by the posterior
#   probabilities of under- and over-staking) are defined in this script,
#   not in the model folder. If the selected folder also contains
#   sensitivity fits of the confirmatory sample (alternative SD floors,
#   MAD scaling), their agreement with the primary fit is summarised.
#
# ESTIMANDS (generated quantities of the RQ2 Stan models)
#   mu_a[s]     mean stake deviation of sequence s from the optimal stake,
#               averaged over the participants of the sample, as a
#               proportion of the endowment
#   mu_a_i[i]   mean stake deviation of participant i, averaged over
#               sequences, as a proportion of the endowment
#   Under-staking is P(mu < -rho) and over-staking P(mu > rho), where rho is
#   cfg$design$rq2$rho (first value = main, the others sensitivity).
#   grand_mean_a is the mean of mu_a over sequences, per posterior draw.
#
# CLASSIFICATION
#   Sequences:     under/over_strong (P >= 0.95), _moderate (>= 0.80),
#                  _weak (>= 0.50); the direction with the larger
#                  probability is used when both reach 0.50; neutral
#   Participants:  under/over_solid (P >= 0.95), _likely (>= 0.90),
#                  _leaning (>= 0.75); neutral
#   Grand-mean classes (sequences): above / likely_above (P(mu_s > grand
#   mean) >= 0.95 / 0.80), below / likely_below, neutral.
#
# SAMPLES
#   full           all participants of each treatment in cfg$run$treatment
#   confirmatory   treatments with cfg$design$a_flags$betting_normative =
#                  TRUE: participants kept by the confirmatory flags
#
# INPUT
#   model "rq2": fit-<tr>-<sample>.rds, pid-levels-<tr>-<sample>.rds,
#                seq-levels-<tr>-<sample>.rds, and, if present,
#                fit-<tr>-confirmatory-<variant>.rds with matching levels
#   model "mpl": sample/analysis-trials.csv, sample/sample-definitions.csv
#                (through analysis_trials(); betting trials), and
#                drift/drift-decisions.rds
#
# OUTPUT (data/output/csv/)
#   rq2-<tr>-<sample>-sequences.csv             one row per sequence
#   rq2-<tr>-<sample>-sequences-summary.csv     sequences by label
#   rq2-<tr>-<sample>-participants.csv          one row per participant
#   rq2-<tr>-<sample>-participants-summary.csv  participants by label
#   rq2-<tr>-<sample>-model-summary.csv         posterior summary of the
#                                               scalar model parameters
#   rq2-<tr>-confirmatory-sensitivity-summary.csv   (only if variant fits exist)
#                                               agreement of each variant with
#                                               the primary fit
#
# NOTES
#   - The model summary lists every scalar parameter of the fit under its
#     Stan name; gamma_drift is listed only when linear drift was active for
#     the treatment and sample, and the column "model" gives the folder
#   - The estimands are checked against the levels files (dimensions) and
#     the scale: absolute values above 2 stop the script; draws beyond 1,
#     possible under the unbounded Gaussian likelihood, are reported in a
#     warning and not clipped
#   - Output files are overwritten at every run
# ============================================================

rq2_tables <- function(cfg) {

  tr_vec <- unique(as.character(cfg$run$treatment))
  design <- cfg$design

  rho_vec  <- as.numeric(design$rq2$rho)
  rho_main <- rho_vec[1]
  rho_nm   <- function(rho) gsub("\\.", "", sprintf("%.2f", rho))
  seq_n    <- as.integer(design$seq$seq_n)

  model_name      <- basename(model_path("rq2"))
  drift_decisions <- readRDS(model_path("mpl", file.path("drift", "drift-decisions.rds")))

  summarise_draws <- function(x, nm) {
    data.table(parameter = nm, median = median(x), mean = mean(x),
               q025 = unname(quantile(x, 0.025)), q975 = unname(quantile(x, 0.975)))
  }

  # ---- Labels: the direction with the larger probability decides ----
  label_by <- function(U, O, cuts, names_under, names_over) {
    lab   <- rep("neutral", length(U))
    under <- U >= cuts[3] & U >= O
    over  <- O >= cuts[3] & O >  U
    for (k in 1:3) {
      lo <- cuts[k]; hi <- if (k == 1) Inf else cuts[k - 1]
      lab[under & U >= lo & U < hi] <- names_under[k]
      lab[over  & O >= lo & O < hi] <- names_over[k]
    }
    lab
  }
  seq_labels <- function(U, O)
    label_by(U, O, c(0.95, 0.80, 0.50),
             paste0("under_", c("strong", "moderate", "weak")),
             paste0("over_",  c("strong", "moderate", "weak")))
  pid_labels <- function(U, O)
    label_by(U, O, c(0.95, 0.90, 0.75),
             paste0("under_", c("solid", "likely", "leaning")),
             paste0("over_",  c("solid", "likely", "leaning")))

  # ---- Tables of the estimands ----
  build_seq_tbl <- function(mu_s, seq_levels, n_by_seq) {
    tbl <- data.table(
      sequence    = seq_levels,
      n_trials    = as.integer(n_by_seq[.(seq_levels), n_trials]),
      mu_a_median = apply(mu_s, 2, median),
      mu_a_mean   = apply(mu_s, 2, mean),
      mu_a_q025   = apply(mu_s, 2, quantile, probs = 0.025),
      mu_a_q975   = apply(mu_s, 2, quantile, probs = 0.975)
    )
    tbl[is.na(n_trials), n_trials := 0L]
    for (rho in rho_vec) {
      tbl[, (paste0("U_a_rho_", rho_nm(rho))) := apply(mu_s, 2, function(x) mean(x < -rho))]
      tbl[, (paste0("O_a_rho_", rho_nm(rho))) := apply(mu_s, 2, function(x) mean(x >  rho))]
    }
    tbl[, calib_label := seq_labels(tbl[[paste0("U_a_rho_", rho_nm(rho_main))]],
                                    tbl[[paste0("O_a_rho_", rho_nm(rho_main))]])]
    grand <- rowMeans(mu_s)
    tbl[, p_above_grand := apply(mu_s, 2, function(x) mean(x > grand))]
    tbl[, p_below_grand := apply(mu_s, 2, function(x) mean(x < grand))]
    tbl[, grand_label := fifelse(p_above_grand >= 0.95, "above",
                         fifelse(p_above_grand >= 0.80, "likely_above",
                         fifelse(p_below_grand >= 0.95, "below",
                         fifelse(p_below_grand >= 0.80, "likely_below", "neutral"))))]
    setorder(tbl, sequence)
    tbl
  }

  build_pid_tbl <- function(mu_i, pid_levels, n_by_pid) {
    tbl <- data.table(
      pid         = pid_levels,
      n_bets      = as.integer(n_by_pid[.(pid_levels), n_bets]),
      mu_a_median = apply(mu_i, 2, median),
      mu_a_mean   = apply(mu_i, 2, mean),
      mu_a_q025   = apply(mu_i, 2, quantile, probs = 0.025),
      mu_a_q975   = apply(mu_i, 2, quantile, probs = 0.975)
    )
    for (rho in rho_vec) {
      tbl[, (paste0("U_i_rho_", rho_nm(rho))) := apply(mu_i, 2, function(x) mean(x < -rho))]
      tbl[, (paste0("O_i_rho_", rho_nm(rho))) := apply(mu_i, 2, function(x) mean(x >  rho))]
    }
    tbl[, calib_label := pid_labels(tbl[[paste0("U_i_rho_", rho_nm(rho_main))]],
                                    tbl[[paste0("O_i_rho_", rho_nm(rho_main))]])]
    setorder(tbl, pid)
    tbl
  }

  # ---- Summary of a table by label ----
  summarise_labels <- function(tbl, levs, mean_col, id_col, id_name) {
    out <- merge(
      tbl[, .N, by = calib_label],
      tbl[, .(mu_mean   = round(mean(get(mean_col)),   3),
              mu_median = round(median(get(mean_col)), 3),
              mu_min    = round(min(get(mean_col)),    3),
              mu_max    = round(max(get(mean_col)),    3),
              ids       = paste(sort(get(id_col)), collapse = ", ")),
          by = calib_label],
      by = "calib_label")
    setnames(out, "ids", id_name)
    out[, calib_label := factor(calib_label, levels = levs)]
    setorder(out, calib_label)
    out[, pct := round(100 * N / nrow(tbl), 1)]
    out
  }

  # ---- Estimand draws of a fit, with their checks ----
  load_fit <- function(tr, sample, variant = NULL) {
    sfx <- paste0(tr, "-", sample, if (!is.null(variant)) paste0("-", variant))
    fit_file <- function(prefix) model_path("rq2", paste0(prefix, "-", sfx, ".rds"))
    fit        <- readRDS(fit_file("fit"))
    pid_levels <- as.character(readRDS(fit_file("pid-levels")))
    seq_levels <- as.character(readRDS(fit_file("seq-levels")))
    stopifnot(length(seq_levels) == seq_n)

    post <- rstan::extract(fit)
    for (nm in c("mu_a", "mu_a_i"))
      if (is.null(post[[nm]])) stop("RQ2 tables: '", nm, "' missing in the fit ", sfx, call. = FALSE)
    stopifnot(is.matrix(post$mu_a),   ncol(post$mu_a)   == length(seq_levels),
              is.matrix(post$mu_a_i), ncol(post$mu_a_i) == length(pid_levels))
    # Scale check: the estimands are proportions of the endowment. The
    # Gaussian likelihood is unbounded, so single draws can exceed 1 for
    # participants who are almost always all-in (reported, not clipped);
    # absolute values above 2 indicate a scale error.
    max_abs <- max(abs(post$mu_a), abs(post$mu_a_i))
    if (!is.finite(max_abs) || max_abs > 2)
      stop("RQ2 tables: estimands of ", sfx, " are not on the scale of the endowment (max |draw| = ",
           signif(max_abs, 3), ").", call. = FALSE)
    n_out <- sum(colMeans(abs(post$mu_a_i) > 1) > 0)
    if (n_out > 0L)
      warning("RQ2 tables: ", sfx, ": draws beyond +-1 for ", n_out, " participant(s) (max |draw| = ",
              round(max_abs, 3), "); not clipped.", call. = FALSE)
    list(post = post, pid_levels = pid_levels, seq_levels = seq_levels)
  }

  # ---- Tables for one treatment and sample ----
  produce_tables <- function(tr, sample) {

    stem <- function(suffix) paste0("rq2-", tr, "-", sample, "-", suffix, ".csv")
    out  <- function(suffix) file.path(path_csv, stem(suffix))

    fit <- load_fit(tr, sample)
    post <- fit$post; pid_levels <- fit$pid_levels; seq_levels <- fit$seq_levels

    # betting trials of the sample
    d <- analysis_trials(tr, paste0("rq2-", sample))
    stopifnot(setequal(unique(d$pid), pid_levels))
    n_by_seq <- d[, .(n_trials = .N), by = seq];  setkey(n_by_seq, seq)
    n_by_pid <- d[, .(n_bets   = .N), by = pid];  setkey(n_by_pid, pid)

    drift <- drift_decisions[[paste(tr, sample, "stake", sep = "_")]]
    drift_active <- !is.null(drift) && isTRUE(drift$drift)

    # ---- Sequences ----
    seq_tbl <- build_seq_tbl(post$mu_a, seq_levels, n_by_seq)
    fwrite(seq_tbl, out("sequences"));  msg("Saved: ", out("sequences"))

    seq_sum <- rbindlist(list(
      summarise_labels(seq_tbl[, .(calib_label = calib_label, mu_a_mean, sequence)],
                       c("under_strong", "under_moderate", "under_weak", "neutral",
                         "over_weak", "over_moderate", "over_strong"),
                       "mu_a_mean", "sequence", "sequences")[, classification := "preregistered"],
      summarise_labels(seq_tbl[, .(calib_label = grand_label, mu_a_mean, sequence)],
                       c("above", "likely_above", "neutral", "likely_below", "below"),
                       "mu_a_mean", "sequence", "sequences")[, classification := "grand_mean"]
    ))
    setnames(seq_sum, "calib_label", "label")
    fwrite(seq_sum, out("sequences-summary"));  msg("Saved: ", out("sequences-summary"))

    # ---- Participants ----
    pid_tbl <- build_pid_tbl(post$mu_a_i, pid_levels, n_by_pid)
    fwrite(pid_tbl, out("participants"));  msg("Saved: ", out("participants"))

    pid_sum <- summarise_labels(pid_tbl,
                                c("under_solid", "under_likely", "under_leaning", "neutral",
                                  "over_leaning", "over_likely", "over_solid"),
                                "mu_a_mean", "pid", "pids")
    setnames(pid_sum, "calib_label", "label")
    fwrite(pid_sum, out("participants-summary"));  msg("Saved: ", out("participants-summary"))

    # ---- Model summary: scalar parameters and the grand mean ----
    scalar <- names(post)[vapply(post, function(x) length(dim(x)) <= 1L, logical(1))]
    scalar <- setdiff(scalar, "lp__")
    if (!drift_active) scalar <- setdiff(scalar, "gamma_drift")

    mod_tbl <- rbindlist(c(
      lapply(scalar, function(nm) summarise_draws(as.numeric(post[[nm]]), nm)),
      list(summarise_draws(rowMeans(post$mu_a), "grand_mean_a"))
    ))
    # (values are set outside data.table's scope, where column names such as
    # "sample" and "drift" would mask the objects of the same name)
    drift_lab <- if (drift_active) drift$drift_type else "none"
    n_tr <- nrow(d); n_pid_total <- length(pid_levels); n_seq_total <- length(seq_levels)
    mod_tbl[, `:=`(treatment = tr, sample = sample, model = model_name,
                   drift = drift_lab, n_trials = n_tr,
                   n_participants = n_pid_total, n_sequences = n_seq_total)]
    fwrite(mod_tbl, out("model-summary"));  msg("Saved: ", out("model-summary"))

    invisible(list(seq_tbl = seq_tbl, pid_tbl = pid_tbl))
  }

  # ---- Agreement of sensitivity fits with the primary fit (confirmatory) ----
  produce_sensitivity <- function(tr, primary) {

    pattern  <- paste0("^fit-", tr, "-confirmatory-(.+)\\.rds$")
    files    <- list.files(model_path("rq2"), pattern = pattern)
    variants <- sub(pattern, "\\1", files)
    if (length(variants) == 0L) return(invisible(NULL))

    rows <- lapply(variants, function(v) {
      fit_v <- load_fit(tr, "confirmatory", v)
      stopifnot(setequal(fit_v$pid_levels, primary$pid_tbl$pid))

      mu_s_v <- apply(fit_v$post$mu_a, 2, mean)
      names(mu_s_v) <- fit_v$seq_levels
      mu_i_v <- apply(fit_v$post$mu_a_i, 2, mean)
      names(mu_i_v) <- fit_v$pid_levels

      n_by_seq <- primary$seq_tbl[, .(seq = sequence, n_trials)]
      setkey(n_by_seq, seq)
      seq_v <- build_seq_tbl(fit_v$post$mu_a, fit_v$seq_levels, n_by_seq)
      stopifnot(identical(seq_v$sequence, primary$seq_tbl$sequence))

      data.table(
        treatment        = tr,
        variant          = v,
        spearman_seq     = round(cor(primary$seq_tbl$mu_a_mean,
                                     mu_s_v[primary$seq_tbl$sequence], method = "spearman"), 4),
        spearman_pid     = round(cor(primary$pid_tbl$mu_a_mean,
                                     mu_i_v[primary$pid_tbl$pid], method = "spearman"), 4),
        pct_stable_label = round(mean(primary$seq_tbl$calib_label == seq_v$calib_label), 4)
      )
    })

    f <- file.path(path_csv, paste0("rq2-", tr, "-confirmatory-sensitivity-summary.csv"))
    fwrite(rbindlist(rows), f);  msg("Saved: ", f)
    invisible(NULL)
  }

  for (tr in tr_vec) produce_tables(tr, "full")
  for (tr in tr_vec) {
    if (!isTRUE(design$a_flags$betting_normative[[tr]])) next
    res <- produce_tables(tr, "confirmatory")
    produce_sensitivity(tr, res)
  }

  invisible(TRUE)
}
