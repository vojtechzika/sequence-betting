# ============================================================
# 23_rq2-stake-vs-r.R
#
# PURPOSE
#   Checks whether the split in individual stake deviation between
#   under-stakers and over-stakers is related to the elicited CRRA risk
#   aversion r_i. The stake benchmark is constructed from r_i, so a link is
#   partly mechanical: a more risk-averse participant has a lower
#   expected-utility-optimal stake and therefore a larger deviation from
#   it. The script correlates r_i with the participant-level stake
#   deviation (Pearson and Spearman) and reports the mean r_i by
#   participant label.
#
# INPUT
#   data/output/csv/rq2-<tr>-<sample>-participants.csv   (21_rq2-tables.R)
#                    posterior median stake deviation and label per pid
#   model "mpl": scored-<tr>.csv    posterior mean of r_i per pid
#
# OUTPUT (data/output/csv/)
#   rq2-<tr>-<sample>-stake-vs-r-summary.csv    n, Pearson r, Spearman rho
#                                               with p-values
#   rq2-<tr>-<sample>-stake-vs-r-by-label.csv   n, mean and SD of r_i, and
#                                               mean stake deviation by label
#
# NOTES
#   - Participants are matched by pid; the participant set is that of the
#     RQ2 sample, not restricted further
#   - The r_i of the full-sample MPL fit are used (participants are not
#     pooled in that model, so these equal the estimates of any subset)
#   - Output files are overwritten at every run
# ============================================================

rq2_stake_vs_r <- function(cfg) {

  tr_vec <- unique(as.character(cfg$run$treatment))
  design <- cfg$design

  produce <- function(tr, sample) {

    stem <- paste0("rq2-", tr, "-", sample)
    f_pid <- file.path(path_csv, paste0(stem, "-participants.csv"))
    stopifnot(file.exists(f_pid))

    part <- fread(f_pid, encoding = "UTF-8")[, .(pid = as.character(pid), mu_a_median, calib_label)]
    mpl  <- fread(model_path("mpl", paste0("scored-", tr, ".csv")),
                  encoding = "UTF-8")[, .(pid = as.character(pid), r_mean)]

    dt <- merge(part, mpl, by = "pid")
    if (nrow(dt) != nrow(part))
      stop("RQ2 stake vs r: participants of ", stem, " missing in the MPL fit.", call. = FALSE)

    pear  <- cor.test(dt$r_mean, dt$mu_a_median, method = "pearson")
    spear <- suppressWarnings(cor.test(dt$r_mean, dt$mu_a_median,
                                       method = "spearman", exact = FALSE))
    summary_tbl <- data.table(
      treatment = tr, sample = sample, n = nrow(dt),
      pearson_r    = unname(pear$estimate),  pearson_p  = pear$p.value,
      spearman_rho = unname(spear$estimate), spearman_p = spear$p.value
    )
    f_sum <- file.path(path_csv, paste0(stem, "-stake-vs-r-summary.csv"))
    fwrite(summary_tbl, f_sum);  msg("Saved: ", f_sum)

    by_label <- dt[, .(n = .N, r_mean_avg = mean(r_mean), r_mean_sd = sd(r_mean),
                       mu_a_avg = mean(mu_a_median)), by = calib_label]
    setorder(by_label, calib_label)
    f_lab <- file.path(path_csv, paste0(stem, "-stake-vs-r-by-label.csv"))
    fwrite(by_label, f_lab);  msg("Saved: ", f_lab)

    invisible(NULL)
  }

  for (tr in tr_vec) produce(tr, "full")
  for (tr in tr_vec) {
    if (!isTRUE(design$a_flags$betting_normative[[tr]])) next
    produce(tr, "confirmatory")
  }

  invisible(TRUE)
}
