# ============================================================
# 09_r-descriptive-summary.R
#
# PURPOSE
#   Summarises the distribution of the posterior mean of the CRRA parameter
#   r_i (primary individual-level Holt-Laury model) across participants.
#   For each treatment the summary is given for all participants, for
#   participants with a consistent Holt-Laury choice pattern and for
#   participants with an inconsistent pattern, and within each of these
#   groups for all participants, women and men.
#
# INPUT
#   model "mpl": scored-<tr>.csv      posterior mean of r_i (primary model) and
#                                     the consistency flag, per treatment fit
#   path_src/participants.csv         sex
#
# OUTPUT (data/output/csv/)
#   r-descriptive-summary.csv   one row per treatment, consistency group and
#                               sex group: number of participants, minimum,
#                               first quartile, median, third quartile,
#                               maximum, mean and standard deviation of r_i
#
# NOTES
#   - The sex group "All" includes participants who do not report their sex
#     (four in the FN treatment), so the rows for women and men need not add
#     up to it
#   - Standard deviations are across participants (not posterior SDs) and are
#     reported only for groups with at least two participants
#   - Estimates of inconsistent participants are bounded above by the
#     indifference point of the last row of the list (r = 1.70), so their
#     distribution accumulates near this ceiling; the upper end of the range
#     is therefore informative about the bound and not about preferences
# ============================================================

r_descriptive_summary <- function(cfg) {
  
  style  <- treatment_style(cfg)
  tr_vec <- style$treat
  
  # ---- Posterior means of r, consistency flag and sex ----
  dt <- rbindlist(lapply(tr_vec, function(tr) {
    sc <- fread(model_path("mpl", paste0("scored-", tr, ".csv")), encoding = "UTF-8")
    sc[, .(pid = as.character(pid), treat = tr, r_mean, inconsistent = as.integer(inconsistent))]
  }))
  
  participants <- fread(file.path(path_src, "participants.csv"), encoding = "UTF-8")
  participants[, pid := as.character(pid)]
  dt <- merge(dt, participants[, .(pid, sex = toupper(trimws(sex)))], by = "pid", all.x = TRUE)
  dt[, sex_group := fcase(sex == "F", "Women", sex == "M", "Men", default = "Other")]
  
  # ---- Summary statistics of a vector of estimates ----
  summarise_r <- function(x) {
    data.table(
      n      = length(x),
      min    = if (length(x)) min(x) else NA_real_,
      q25    = if (length(x)) unname(quantile(x, 0.25)) else NA_real_,
      median = if (length(x)) median(x) else NA_real_,
      q75    = if (length(x)) unname(quantile(x, 0.75)) else NA_real_,
      max    = if (length(x)) max(x) else NA_real_,
      mean   = if (length(x)) mean(x) else NA_real_,
      sd     = if (length(x) >= 2L) sd(x) else NA_real_
    )
  }
  
  res <- rbindlist(lapply(tr_vec, function(tr) {
    rbindlist(lapply(c("all", "consistent", "inconsistent"), function(cs) {
      d <- dt[treat == tr]
      if (cs == "consistent")   d <- d[inconsistent == 0L]
      if (cs == "inconsistent") d <- d[inconsistent == 1L]
      rbindlist(lapply(c("All", "Women", "Men"), function(sx) {
        x <- if (sx == "All") d$r_mean else d[sex_group == sx, r_mean]
        cbind(data.table(treatment = tr, consistency = cs, sex = sx), summarise_r(x))
      }))
    }))
  }))
  
  f_out <- file.path(path_csv, "r-descriptive-summary.csv")
  fwrite(res, f_out)
  msg("Saved: ", f_out)
  
  invisible(res)
}
