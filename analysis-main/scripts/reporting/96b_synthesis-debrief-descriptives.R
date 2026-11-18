# ============================================================
# 96b_synthesis-debrief-descriptives.R
#
# PURPOSE
#   Descriptive statistics of the closed debriefing questions (Q1 to Q7) and the
#   distribution of the answers to Q8 (used a rule or strategy), by treatment
#   and pooled.
#
# SAMPLE
#   All participants who completed the debriefing questionnaire (FN = m25,
#   FP = m19).
#
# INPUT
#   data/source/debriefing.csv          questionnaire, one row per participant
#   data/derived/master-sequences.csv   pid, treat (treatment of each participant)
#
# OUTPUT (data/output/csv/)
#   synthesis-debrief-descriptives.csv
#     one row per item and sample (m25, m19, all): n, mean, sd, median, min, max
#   synthesis-debrief-used-strategy-shares.csv
#     one row per answer to Q8 and sample: n, share
#
# NOTES
#   - Q1 is a number between 0 and 100, Q2 to Q7 are answered on a 0 to 10 scale
#   - sd is the sample standard deviation
#   - Output files are overwritten at every run
# ============================================================

synthesis_debrief_descriptives <- function(cfg) {

  q <- fread(file.path(path_src, "debriefing.csv"), encoding = "UTF-8")
  m <- unique(fread(file.path(path_der, "master-sequences.csv"), encoding = "UTF-8",
                    select = c("pid", "treat")))
  stopifnot(!anyDuplicated(m$pid))
  q <- merge(q, m, by = "pid", all.x = TRUE)
  stopifnot(!anyNA(q$treat))

  items <- c("perceived_win_probability", "perceived_risk", "belief_independence",
             "reliance_on_sequence", "action_seeking", "fatigue", "self_risk_tolerance")

  samples <- list(m25 = q[treat == "m25"], m19 = q[treat == "m19"], all = q)

  # Moments of the closed items for each sample
  desc <- rbindlist(lapply(names(samples), function(s) {
    x <- samples[[s]]
    rbindlist(lapply(items, function(it) {
      v <- as.numeric(x[[it]])
      data.table(sample = s, item = it, n = sum(!is.na(v)), mean = mean(v, na.rm = TRUE),
                 sd = sd(v, na.rm = TRUE), median = median(v, na.rm = TRUE),
                 min = min(v, na.rm = TRUE), max = max(v, na.rm = TRUE))
    }))
  }))
  f1 <- file.path(path_csv, "synthesis-debrief-descriptives.csv")
  fwrite(desc, f1);  msg("Saved: ", f1)

  # Answers to Q8
  shares <- rbindlist(lapply(names(samples), function(s) {
    x <- samples[[s]]
    t <- x[, .(n = .N), by = .(answer = used_strategy)]
    t[, `:=`(sample = s, share = n / sum(n))]
    t
  }))
  f2 <- file.path(path_csv, "synthesis-debrief-used-strategy-shares.csv")
  fwrite(shares, f2);  msg("Saved: ", f2)
  invisible(list(desc = desc, shares = shares))
}
