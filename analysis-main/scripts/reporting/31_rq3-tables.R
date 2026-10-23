# ============================================================
# 31_rq3-tables.R
#
# PURPOSE
#   Sequence-level, participant-level and model summary tables for RQ3
#   (certainty-equivalent welfare loss from suboptimal betting), from the
#   model selected under the key "rq3" in selected-models.csv: hurdle-Gamma,
#   Gamma-only or Gaussian. Labels and thresholds (classification of sequences and
#   participants by the posterior probability of a welfare loss above rho)
#   are defined in this script, not in the model folder.
#
# ESTIMANDS (generated quantities of the RQ3 Stan models)
#   mu_c[s]     mean welfare loss of sequence s, averaged over the
#               participants of the sample, as a proportion of the endowment
#   mu_c_i[i]   mean welfare loss of participant i, averaged over sequences,
#               as a proportion of the endowment
#   L_rho = P(mu > rho), where rho is cfg$design$rq3$rho (first value =
#   main, the others sensitivity). grand_mean_c is the mean of mu_c over
#   sequences, per posterior draw.
#
# CLASSIFICATION
#   Sequences:     strong (L >= 0.95), moderate (>= 0.80), weak (>= 0.50),
#                  neutral
#   Participants:  solid (L >= 0.95), likely (>= 0.90), leaning (>= 0.75),
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
#   model "rq3": fit-<tr>-<sample>.rds,
#                pid-levels-<tr>-<sample>.rds, seq-levels-<tr>-<sample>.rds
#   model "mpl": sample/analysis-trials.csv, sample/sample-definitions.csv
#                (through analysis_trials(); all trials), and
#                drift/drift-decisions.rds
#   data/output/csv/rq1-<tr>-<sample>-participants.csv   (11_rq1-tables.R),
#                posterior median betting probability, merged if present
#
# OUTPUT (data/output/csv/)
#   rq3-<tr>-<sample>-sequences.csv             one row per sequence
#   rq3-<tr>-<sample>-sequences-summary.csv     sequences by label
#   rq3-<tr>-<sample>-participants.csv          one row per participant
#   rq3-<tr>-<sample>-participants-summary.csv  participants by label
#   rq3-<tr>-<sample>-model-summary.csv         posterior summary of the
#                                               scalar model parameters
#
# NOTES
#   - The model summary lists every scalar parameter of the fit under its
#     Stan name; gamma_drift is listed only when linear drift was active for
#     the treatment and sample. The column "model" gives the folder and
#     "likelihood" the likelihood (hurdle-gamma, gamma, gaussian)
#   - The number of trials of a sequence or participant counts all trials of
#     the sample, including those a likelihood does not use (zero losses in
#     the Gamma-only model)
#   - The estimands are checked against the levels files (dimensions) and
#     the scale: absolute values above 2 stop the script; draws outside
#     [0, 1], possible under the unbounded Gaussian likelihood, are reported
#     in a warning and not clipped
#   - Output files are overwritten at every run
# ============================================================

rq3_tables <- function(cfg) {

  tr_vec <- unique(as.character(cfg$run$treatment))
  design <- cfg$design

  rho_vec  <- as.numeric(design$rq3$rho)
  rho_main <- rho_vec[1]
  rho_nm   <- function(rho) gsub("\\.", "", sprintf("%.2f", rho))
  seq_n    <- as.integer(design$seq$seq_n)

  drift_decisions <- readRDS(model_path("mpl", file.path("drift", "drift-decisions.rds")))

  summarise_draws <- function(x, nm) {
    data.table(parameter = nm, median = median(x), mean = mean(x),
               q025 = unname(quantile(x, 0.025)), q975 = unname(quantile(x, 0.975)))
  }

  # ---- Labels from the probability of a loss above rho ----
  label_by <- function(L, cuts, names) {
    lab <- rep("neutral", length(L))
    # cuts ascending: later (higher) cut-offs overwrite earlier ones
    for (k in seq_along(cuts)) lab[L >= cuts[k]] <- names[k]
    lab
  }
  seq_labels <- function(L) label_by(L, c(0.50, 0.80, 0.95), c("weak", "moderate", "strong"))
  pid_labels <- function(L) label_by(L, c(0.75, 0.90, 0.95), c("leaning", "likely", "solid"))

  # ---- Summary of a table by label ----
  summarise_labels <- function(tbl, label_col, levs, id_col, id_name) {
    out <- merge(
      tbl[, .N, by = label_col],
      tbl[, .(mu_mean   = round(mean(mu_c_mean),   3),
              mu_median = round(median(mu_c_mean), 3),
              mu_min    = round(min(mu_c_mean),    3),
              mu_max    = round(max(mu_c_mean),    3),
              ids       = paste(sort(get(id_col)), collapse = ", ")),
          by = label_col],
      by = label_col)
    setnames(out, c(label_col, "ids"), c("label", id_name))
    out[, label := factor(label, levels = levs)]
    setorder(out, label)
    out[, pct := round(100 * N / nrow(tbl), 1)]
    out
  }

  # ---- Tables for one treatment and sample ----
  produce_tables <- function(tr, sample) {

    key <- "rq3"
    model_name <- basename(model_path(key))
    likelihood <- sub("^rq3-", "", model_name)

    stem <- function(suffix) paste0("rq3-", tr, "-", sample, "-", suffix, ".csv")
    out  <- function(suffix) file.path(path_csv, stem(suffix))

    fit_file   <- function(prefix) model_path(key, paste0(prefix, "-", tr, "-", sample, ".rds"))
    fit        <- readRDS(fit_file("fit"))
    pid_levels <- as.character(readRDS(fit_file("pid-levels")))
    seq_levels <- as.character(readRDS(fit_file("seq-levels")))
    stopifnot(length(seq_levels) == seq_n)

    post <- rstan::extract(fit)
    for (nm in c("mu_c", "mu_c_i"))
      if (is.null(post[[nm]])) stop("RQ3 tables: '", nm, "' missing in the fit ", tr, "/", sample, call. = FALSE)
    mu_s <- post$mu_c
    mu_i <- post$mu_c_i
    stopifnot(is.matrix(mu_s), ncol(mu_s) == length(seq_levels),
              is.matrix(mu_i), ncol(mu_i) == length(pid_levels))
    # Scale check: the estimands are proportions of the endowment. The
    # Gamma and Gaussian likelihoods have no upper bound at 1, so for
    # participants with a large participant effect single draws of the
    # participant-level mean can exceed 1 (reported, not clipped). The
    # sequence-level means average over participants and must stay on the
    # scale of the endowment: absolute values above 2 indicate a scale error.
    max_seq <- max(abs(mu_s))
    max_pid <- max(abs(mu_i))
    if (!is.finite(max_seq) || !is.finite(max_pid) || max_seq > 2)
      stop("RQ3 tables: sequence-level estimands of ", tr, "/", sample, " are not on the scale of the endowment (max |draw| = ",
           signif(max_seq, 3), ").", call. = FALSE)
    n_out <- sum(colMeans(mu_i < 0 | mu_i > 1) > 0)
    if (n_out > 0L)
      warning("RQ3 tables: ", tr, "/", sample, ": draws outside [0, 1] for ", n_out,
              " participant(s) (max |draw| = ", round(max_pid, 3), "); not clipped.", call. = FALSE)

    # trials of the sample
    d <- analysis_trials(tr, paste0("rq3-", sample))
    stopifnot(setequal(unique(d$pid), pid_levels))
    n_by_seq <- d[, .(n_trials = .N), by = seq];  setkey(n_by_seq, seq)
    n_by_pid <- d[, .(n_trials = .N), by = pid];  setkey(n_by_pid, pid)

    drift <- drift_decisions[[paste(tr, sample, "stake", sep = "_")]]
    drift_active <- !is.null(drift) && isTRUE(drift$drift)

    grand <- rowMeans(mu_s)

    # ---- Sequences ----
    seq_tbl <- data.table(
      sequence    = seq_levels,
      n_trials    = as.integer(n_by_seq[.(seq_levels), n_trials]),
      mu_c_median = apply(mu_s, 2, median),
      mu_c_mean   = apply(mu_s, 2, mean),
      mu_c_q025   = apply(mu_s, 2, quantile, probs = 0.025),
      mu_c_q975   = apply(mu_s, 2, quantile, probs = 0.975)
    )
    seq_tbl[is.na(n_trials), n_trials := 0L]
    for (rho in rho_vec)
      seq_tbl[, (paste0("L_rho_", rho_nm(rho))) := apply(mu_s, 2, function(x) mean(x > rho))]
    seq_tbl[, loss_label := seq_labels(seq_tbl[[paste0("L_rho_", rho_nm(rho_main))]])]
    seq_tbl[, p_above_grand := apply(mu_s, 2, function(x) mean(x > grand))]
    seq_tbl[, p_below_grand := apply(mu_s, 2, function(x) mean(x < grand))]
    seq_tbl[, grand_label := fifelse(p_above_grand >= 0.95, "above",
                             fifelse(p_above_grand >= 0.80, "likely_above",
                             fifelse(p_below_grand >= 0.95, "below",
                             fifelse(p_below_grand >= 0.80, "likely_below", "neutral"))))]
    setorder(seq_tbl, sequence)
    fwrite(seq_tbl, out("sequences"));  msg("Saved: ", out("sequences"))

    seq_sum <- rbindlist(list(
      summarise_labels(seq_tbl, "loss_label",
                       c("strong", "moderate", "weak", "neutral"),
                       "sequence", "sequences")[, classification := "preregistered"],
      summarise_labels(seq_tbl, "grand_label",
                       c("above", "likely_above", "neutral", "likely_below", "below"),
                       "sequence", "sequences")[, classification := "grand_mean"]
    ))
    fwrite(seq_sum, out("sequences-summary"));  msg("Saved: ", out("sequences-summary"))

    # ---- Participants ----
    pid_tbl <- data.table(
      pid         = pid_levels,
      n_trials    = as.integer(n_by_pid[.(pid_levels), n_trials]),
      mu_c_median = apply(mu_i, 2, median),
      mu_c_mean   = apply(mu_i, 2, mean),
      mu_c_q025   = apply(mu_i, 2, quantile, probs = 0.025),
      mu_c_q975   = apply(mu_i, 2, quantile, probs = 0.975)
    )
    for (rho in rho_vec)
      pid_tbl[, (paste0("L_rho_", rho_nm(rho))) := apply(mu_i, 2, function(x) mean(x > rho))]
    pid_tbl[, loss_label := pid_labels(pid_tbl[[paste0("L_rho_", rho_nm(rho_main))]])]
    setorder(pid_tbl, pid)

    # posterior median betting probability of RQ1, if the table exists
    f_rq1 <- file.path(path_csv, paste0("rq1-", tr, "-", sample, "-participants.csv"))
    if (file.exists(f_rq1)) {
      rq1 <- fread(f_rq1, encoding = "UTF-8", colClasses = c(pid = "character"))
      pid_tbl <- merge(pid_tbl, rq1[, .(pid, mu_b_median = mu_i_median)], by = "pid", all.x = TRUE)
      setorder(pid_tbl, pid)
    }
    fwrite(pid_tbl, out("participants"));  msg("Saved: ", out("participants"))

    pid_sum <- summarise_labels(pid_tbl, "loss_label",
                                c("solid", "likely", "leaning", "neutral"), "pid", "pids")
    fwrite(pid_sum, out("participants-summary"));  msg("Saved: ", out("participants-summary"))

    # ---- Model summary: scalar parameters and the grand mean ----
    scalar <- names(post)[vapply(post, function(x) length(dim(x)) <= 1L, logical(1))]
    scalar <- setdiff(scalar, "lp__")
    if (!drift_active) scalar <- setdiff(scalar, "gamma_drift")

    mod_tbl <- rbindlist(c(
      lapply(scalar, function(nm) summarise_draws(as.numeric(post[[nm]]), nm)),
      list(summarise_draws(grand, "grand_mean_c"))
    ))
    # (values are set outside data.table's scope, where column names such as
    # "sample" and "drift" would mask the objects of the same name)
    drift_lab <- if (drift_active) drift$drift_type else "none"
    n_tr <- nrow(d); n_pid_total <- length(pid_levels); n_seq_total <- length(seq_levels)
    mod_tbl[, `:=`(treatment = tr, sample = sample, model = model_name,
                   likelihood = likelihood, drift = drift_lab, n_trials = n_tr,
                   n_participants = n_pid_total, n_sequences = n_seq_total)]
    fwrite(mod_tbl, out("model-summary"));  msg("Saved: ", out("model-summary"))

    invisible(NULL)
  }

  for (tr in tr_vec) produce_tables(tr, "full")
  for (tr in tr_vec) {
    if (!isTRUE(design$a_flags$betting_normative[[tr]])) next
    produce_tables(tr, "confirmatory")
  }

  invisible(TRUE)
}
