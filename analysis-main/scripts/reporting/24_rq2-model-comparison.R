# ============================================================
# 24_rq2-model-comparison.R
#
# PURPOSE
#   Compares the RQ2 stake models on the same sample: the Gaussian model
#   (21_rq2-stan-gaussian.R), the preregistered one-inflated Beta-Binomial
#   model (22_rq2-stan-one-inflated.R) and its extension with a sequence
#   effect on the all-in probability (23_rq2-stan-one-inflated-seq.R). The
#   script reads the model folders directly, not the key "rq2"
#   in selected-models.csv, so the comparison does not depend on which model
#   is selected for the tables and figures. A model whose folder or fit for
#   a sample is missing is left out of that comparison.
#
#   For every sample the script reports, by model, (i) the posterior of the
#   sequence-level stake deviation mu_a against its raw counterpart, the mean
#   deviation observed in the data, (ii) the grand mean and the agreement of
#   the sequence ordering with the raw ordering and with the Gaussian model,
#   (iii) posterior predictive checks of the all-in share, and (iv) the
#   convergence status. Pure sequences (one symbol repeated) are shown
#   separately because they carry the largest effects.
#
# RAW COUNTERPART
#   For each sequence, the mean over participants of the participant's mean
#   stake deviation (stake minus the participant's optimal stake, as a
#   proportion of the endowment) on the betting trials of the sample. It is
#   not shrunk by partial pooling, so part of a gap to the model estimate is
#   expected; the gap relative to the posterior interval is reported.
#
# MODEL-FREE ANCHOR
#   Participant-clustered bootstrap (4000 resamples of the participants) of
#   the raw statistics: the mean stake deviation over the pure sequences and
#   over the other sequences, their difference, each pure sequence, and the
#   all-in share in pure and in other sequences and their difference. The
#   interval of the pure-sequence mean is added to the model summary, with a
#   flag whether each model's estimate lies inside it.
#
# POSTERIOR PREDICTIVE CHECKS
#   Taken from ppc-summary.csv of each model folder: overall all-in share
#   (preregistered), all-in share in pure minus other sequences, and the SD
#   across sequences of the all-in share. For one-inflated folders built
#   before the last two statistics were added, they are computed here from
#   the fit with rq2_ppc_stats() (_rq2-ppc.R), with the seeds of the model
#   scripts, so the values equal those the model script would write. The
#   Gaussian model has no all-in component and carries only the overall
#   check.
#
# SAMPLES
#   full           all participants of each treatment in cfg$run$treatment
#   confirmatory   treatments with cfg$design$a_flags$betting_normative =
#                  TRUE: participants kept by the confirmatory flags
#
# INPUT
#   data/models/rq2-gaussian, rq2-one-inflated, rq2-one-inflated-seq:
#     fit-, seq-levels-, pid-levels- and data-<tr>-<sample>.rds,
#     ppc-summary.csv, diagnostics.csv
#   model "mpl": sample/analysis-trials.csv, sample/sample-definitions.csv
#                (through analysis_trials()), drift/drift-decisions.rds
#
# OUTPUT (data/output/csv/, data/output/figures/)
#   rq2-model-comparison-<tr>-<sample>-sequences.csv   one row per sequence:
#       raw mean, and mean and 95% interval of mu_a for each model
#   rq2-model-comparison-<tr>-<sample>-raw-bootstrap.csv  raw statistics with
#       participant-clustered bootstrap intervals
#   rq2-model-comparison-<tr>-<sample>-summary.csv     one row per model
#   rq2-model-comparison-<tr>-<sample>-raw-vs-model.png
#       raw against model estimates, one panel per model
#
# NOTES
#   - Convergence: all of rhat_max, ess_min, divergences, treedepth_hits and
#     ebfmi_min pass in diagnostics.csv of the model
#   - sigma_c is reported for the model with a sequence effect on the all-in
#     probability only
#   - Output files are overwritten at every run
# ============================================================

rq2_model_comparison <- function(cfg) {

  tr_vec <- unique(as.character(cfg$run$treatment))
  design <- cfg$design
  e      <- as.numeric(design$seq$endowment)
  seed   <- as.integer(cfg$run$seed)
  K_ppc  <- as.integer(cfg$model$ppc$rq2_k)
  cut    <- as.numeric(design$rq2$ppc_cut)

  if (!exists("analysis_trials", mode = "function"))
    source(here::here("scripts", "modelling", "09_analysis-sample.R"))
  source(here::here("scripts", "modelling", "_rq2-ppc.R"))

  drift_decisions <- readRDS(model_path("mpl", file.path("drift", "drift-decisions.rds")))

  models <- data.table(
    folder = c("rq2-gaussian", "rq2-one-inflated", "rq2-one-inflated-seq"),
    name   = c("gaussian", "one_inflated", "one_inflated_seq"),
    label  = c("Gaussian", "One-inflated (preregistered)",
               "One-inflated, sequence effect on all-in")
  )
  models <- models[dir.exists(file.path(path_mod, folder))]
  stopifnot(nrow(models) > 0L)

  conv_checks <- c("rhat_max", "ess_min", "divergences", "treedepth_hits", "ebfmi_min")

  # ---- All-in checks of one model: from the folder, or computed here ----
  get_ppc <- function(folder, tr, sample, sfx) {
    f <- file.path(path_mod, folder, "ppc-summary.csv")
    smp <- sample   # (a column of the same name masks `sample` inside data.table)
    ppc <- if (file.exists(f)) fread(f)[treatment == tr & sample == smp] else NULL
    need <- c("allin", "pure_contrast", "seq_sd")
    if (!is.null(ppc) && all(need %in% ppc$statistic)) return(ppc)
    if (folder == "rq2-gaussian") return(ppc)         # overall check only

    fit        <- readRDS(model_path(folder, paste0("fit-", sfx, ".rds")))
    d          <- readRDS(model_path(folder, paste0("data-", sfx, ".rds")))
    pid_levels <- readRDS(model_path(folder, paste0("pid-levels-", sfx, ".rds")))
    seq_levels <- readRDS(model_path(folder, paste0("seq-levels-", sfx, ".rds")))
    drift <- drift_decisions[[paste(tr, sample, "stake", sep = "_")]]
    rq2_ppc_stats(fit, d, pid_levels, seq_levels, !is.null(drift) && isTRUE(drift$drift), e,
                  K_ppc, seed, cut, tr, sample, run = NULL,
                  model_label = paste("RQ2", folder),
                  seq_effect = grepl("-seq$", folder), figure = FALSE)
  }

  # ---- Model-free anchor: participant-clustered bootstrap of the raw statistics ----
  # Participants are resampled with replacement; for each resample the raw
  # sequence means, the means over pure and over other sequences, and the
  # all-in shares are recomputed (same definitions as for the data).
  B_boot <- 4000L
  raw_bootstrap <- function(d, seq_levels, pure) {

    pids <- sort(unique(d$pid)); N <- length(pids); S <- length(seq_levels)
    i <- match(d$pid, pids);     s <- match(d$seq, seq_levels)
    allin  <- as.numeric(d$stake >= e)
    tr_pure <- as.numeric(pure[s])

    # participant x sequence mean deviation, and which cells exist
    cell <- d[, .(dev = mean((stake - a_star) / e)), by = .(pid, seq)]
    ix <- cbind(match(cell$pid, pids), match(cell$seq, seq_levels))
    cm  <- matrix(0, N, S);  cm[ix]  <- cell$dev
    has <- matrix(0, N, S);  has[ix] <- 1

    # all-in counts per participant in pure and in other sequences
    cnt <- function(x) as.vector(rowsum(x, i, reorder = TRUE))
    a_p <- cnt(allin * tr_pure);  n_p <- cnt(tr_pure)
    a_o <- cnt(allin * (1 - tr_pure));  n_o <- cnt(1 - tr_pure)

    stats_of <- function(W) {            # W: resamples x participants (weights)
      sm <- (W %*% (cm * has)) / (W %*% has)                   # resamples x sequences
      out <- list(dev_pure_mean  = rowMeans(sm[, pure,  drop = FALSE]),
                  dev_other_mean = rowMeans(sm[, !pure, drop = FALSE]))
      out$dev_pure_minus_other <- out$dev_pure_mean - out$dev_other_mean
      out$allin_pure     <- as.vector(W %*% a_p) / as.vector(W %*% n_p)
      out$allin_other    <- as.vector(W %*% a_o) / as.vector(W %*% n_o)
      out$allin_contrast <- out$allin_pure - out$allin_other
      for (k in which(pure)) out[[paste0("dev_", seq_levels[k])]] <- sm[, k]
      out
    }

    est <- stats_of(matrix(1, 1, N))
    set.seed(seed + 4001L)
    W   <- t(rmultinom(B_boot, N, rep(1 / N, N)))
    bs  <- stats_of(W)
    rbindlist(lapply(names(est), function(nm) data.table(
      statistic = nm, estimate = est[[nm]],
      boot_q025 = unname(quantile(bs[[nm]], 0.025, na.rm = TRUE)),
      boot_q975 = unname(quantile(bs[[nm]], 0.975, na.rm = TRUE)),
      boot_sd   = sd(bs[[nm]], na.rm = TRUE), n_boot = B_boot)))
  }

  compare_set <- function(tr, sample) {

    sfx   <- paste0(tr, "-", sample)
    avail <- models[file.exists(file.path(path_mod, folder, paste0("fit-", sfx, ".rds")))]
    if (nrow(avail) == 0L) return(invisible(NULL))

    # ---- Raw sequence means on the betting trials of the sample ----
    d <- analysis_trials(tr, paste0("rq2-", sample))
    cell <- d[, .(dev = mean((stake - a_star) / e)), by = .(pid, seq)]
    raw  <- merge(cell[, .(raw_mean = mean(dev)), by = seq],
                  d[, .(n_trials = .N), by = seq], by = "seq")
    setnames(raw, "seq", "sequence")
    raw[, pure := grepl("^(.)\\1*$", sequence)]
    setorder(raw, sequence)

    # ---- Model-free anchor ----
    boot <- raw_bootstrap(d, raw$sequence, raw$pure)
    boot[, `:=`(treatment = tr, sample = sample)]
    f_boot <- file.path(path_csv, paste0("rq2-model-comparison-", sfx, "-raw-bootstrap.csv"))
    fwrite(boot, f_boot);  msg("Saved: ", f_boot)
    raw_pure_ci <- boot[statistic == "dev_pure_mean", .(boot_q025, boot_q975)]

    # ---- Posterior of mu_a by model ----
    fits <- lapply(seq_len(nrow(avail)), function(j) {
      fol  <- avail$folder[j]
      fit  <- readRDS(model_path(fol, paste0("fit-", sfx, ".rds")))
      lev  <- as.character(readRDS(model_path(fol, paste0("seq-levels-", sfx, ".rds"))))
      post <- rstan::extract(fit)
      stopifnot(is.matrix(post$mu_a), ncol(post$mu_a) == length(lev),
                setequal(lev, raw$sequence))
      mu <- post$mu_a[, match(raw$sequence, lev), drop = FALSE]
      list(mu = mu, sigma_c = post$sigma_c)
    })
    names(fits) <- avail$name

    seq_out <- copy(raw)
    for (nm in avail$name) {
      mu <- fits[[nm]]$mu
      seq_out[, (paste0(nm, "_mean"))  := apply(mu, 2, mean)]
      seq_out[, (paste0(nm, "_q025"))  := apply(mu, 2, quantile, probs = 0.025)]
      seq_out[, (paste0(nm, "_q975"))  := apply(mu, 2, quantile, probs = 0.975)]
    }
    f_seq <- file.path(path_csv, paste0("rq2-model-comparison-", sfx, "-sequences.csv"))
    fwrite(seq_out, f_seq);  msg("Saved: ", f_seq)

    # ---- One row per model ----
    rows <- lapply(seq_len(nrow(avail)), function(j) {
      nm  <- avail$name[j]
      fol <- avail$folder[j]
      m   <- seq_out[[paste0(nm, "_mean")]]
      lo  <- seq_out[[paste0(nm, "_q025")]]
      hi  <- seq_out[[paste0(nm, "_q975")]]
      pur <- seq_out$pure
      grand <- rowMeans(fits[[nm]]$mu)

      ppc <- get_ppc(fol, tr, sample, sfx)
      pv  <- function(s) if (is.null(ppc) || !(s %in% ppc$statistic)) NA_real_ else ppc[statistic == s, p_value]

      dg <- file.path(path_mod, fol, "diagnostics.csv")
      conv <- if (file.exists(dg)) {
        tg <- paste(tr, sample, sep = "_")
        x <- fread(dg)[tag == tg & check %in% conv_checks]
        nrow(x) == length(conv_checks) && all(as.logical(x$pass))
      } else NA

      sc <- fits[[nm]]$sigma_c
      data.table(
        treatment = tr, sample = sample, model = nm, folder = fol, label = avail$label[j],
        grand_mean = mean(grand), grand_q025 = unname(quantile(grand, 0.025)),
        grand_q975 = unname(quantile(grand, 0.975)),
        raw_pure_mean = mean(seq_out$raw_mean[pur]), pure_mean = mean(m[pur]),
        pure_gap = mean(seq_out$raw_mean[pur]) - mean(m[pur]),
        raw_pure_q025 = raw_pure_ci$boot_q025, raw_pure_q975 = raw_pure_ci$boot_q975,
        pure_mean_in_raw_ci = mean(m[pur]) >= raw_pure_ci$boot_q025 & mean(m[pur]) <= raw_pure_ci$boot_q975,
        raw_other_mean = mean(seq_out$raw_mean[!pur]), other_mean = mean(m[!pur]),
        cor_other = cor(m[!pur], seq_out$raw_mean[!pur]),
        rmse_other = sqrt(mean((m[!pur] - seq_out$raw_mean[!pur])^2)),
        spearman_raw = cor(m, seq_out$raw_mean, method = "spearman"),
        n_raw_outside = sum(seq_out$raw_mean < lo | seq_out$raw_mean > hi),
        n_pure_raw_outside = sum((seq_out$raw_mean < lo | seq_out$raw_mean > hi) & pur),
        sigma_c_median = if (is.null(sc)) NA_real_ else median(sc),
        sigma_c_q025   = if (is.null(sc)) NA_real_ else unname(quantile(sc, 0.025)),
        sigma_c_q975   = if (is.null(sc)) NA_real_ else unname(quantile(sc, 0.975)),
        ppc_allin_p = pv("allin"), ppc_pure_contrast_p = pv("pure_contrast"),
        ppc_seq_sd_p = pv("seq_sd"), convergence_pass = conv
      )
    })
    summ <- rbindlist(rows)
    if ("gaussian" %in% avail$name) {
      g <- seq_out$gaussian_mean
      summ[, spearman_vs_gaussian := vapply(model, function(nm)
        cor(seq_out[[paste0(nm, "_mean")]], g, method = "spearman"), numeric(1))]
    }
    f_sum <- file.path(path_csv, paste0("rq2-model-comparison-", sfx, "-summary.csv"))
    fwrite(summ, f_sum);  msg("Saved: ", f_sum)

    # ---- Figure: raw against model estimates ----
    long <- rbindlist(lapply(seq_len(nrow(avail)), function(j) {
      nm <- avail$name[j]
      data.table(model = factor(avail$label[j], levels = avail$label),
                 sequence = seq_out$sequence, pure = seq_out$pure, raw = seq_out$raw_mean,
                 est = seq_out[[paste0(nm, "_mean")]],
                 lo = seq_out[[paste0(nm, "_q025")]], hi = seq_out[[paste0(nm, "_q975")]])
    }))
    long[, kind := factor(ifelse(pure, "Pure sequences", "Other sequences"),
                          levels = c("Other sequences", "Pure sequences"))]
    p <- ggplot(long, aes(x = raw, y = est, colour = kind)) +
      geom_abline(slope = 1, intercept = 0, colour = "grey50", linetype = "dashed", linewidth = 0.4) +
      geom_linerange(aes(ymin = lo, ymax = hi), linewidth = 0.4, alpha = 0.5) +
      geom_point(size = 1.6) +
      facet_wrap(~model, nrow = 1) +
      scale_colour_manual(values = c("Other sequences" = "#2166AC", "Pure sequences" = "#C0392B"),
                          name = NULL) +
      scale_x_continuous(name = "Raw mean stake deviation of the sequence",
                         labels = scales::percent_format(accuracy = 1)) +
      scale_y_continuous(name = expression(hat(mu)[s]^a ~ "(posterior mean, 95% interval)"),
                         labels = scales::percent_format(accuracy = 1)) +
      theme_classic(base_size = 10) +
      theme(strip.background = element_blank(), strip.text = element_text(size = 8.5),
            legend.position = "bottom", panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
            plot.margin = margin(6, 10, 6, 6))
    f_fig <- file.path(path_fig, paste0("rq2-model-comparison-", sfx, "-raw-vs-model.png"))
    ggsave(f_fig, p, width = 3.4 * nrow(avail) + 0.6, height = 4, dpi = 300)
    msg("Saved: ", f_fig)

    invisible(NULL)
  }

  for (tr in tr_vec) compare_set(tr, "full")
  for (tr in tr_vec) {
    if (!isTRUE(design$a_flags$betting_normative[[tr]])) next
    compare_set(tr, "confirmatory")
  }

  invisible(TRUE)
}
