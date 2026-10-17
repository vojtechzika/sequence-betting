# ============================================================
# 52_ex1-pure-sequences.R
#
# PURPOSE
#   Descriptive treatment of the pure sequences HHHHHH and TTTTTT, used in
#   place of the anchor-based index, which does not discriminate between
#   sequences (see 51_ex1-weights-scatter.R). Each participant is classified
#   from the side chosen on the two pure sequences:
#       hothand_consistent   Heads on HHHHHH and Tails on TTTTTT (with the streak)
#       gambler_consistent   Tails on HHHHHH and Heads on TTTTTT (against it)
#       inconsistent         the same side on both sequences
#       unclassified         no bet on at least one of the two sequences
#   The classification is the pure-sequence score chi_pure of the EX1
#   specification (+1 hot-hand, -1 gambler, 0 inconsistent; unclassified
#   missing). The script then
#     (1) counts the types by treatment,
#     (2) reports mean standardised optimism (LOT-R), log response time and
#         risk aversion r with 95% t-intervals by type,
#     (3) cross-tabulates the type against behaviour on four mirror pairs of
#         sequences ending in a terminal run, and reports the share of Heads
#         on the two pure sequences by type.
#
#   On a mirror sequence the response is gamblerish if the participant bets
#   against the terminal run and hot-handish otherwise. A participant is
#   gamblerish (hot-handish) on the mirrors if the mean score over their
#   bets on the eight sequences is positive (negative), and split if zero.
#
# SAMPLE
#   All participants of each treatment in cfg$run$treatment (sample rq1-full
#   defines the participants); side choices come from betting trials.
#
# INPUT
#   analysis_trials(tr, "rq1-full")     trials, side, stake, response time
#   model "mpl": scored-<tr>.csv        posterior mean of r_i
#   path_der/lotr-scored.csv            LOT-R score
#
# OUTPUT (data/output/csv/)
#   ex1-<tr>-pure-participants.csv      pid, chi_pure, type, mirror score and type
#   ex1-pure-counts.csv                 number and share of each type
#   ex1-pure-covariates.csv             covariate means, SD and 95% CI by type
#   ex1-pure-vs-mirror.csv              cross-tabulation, with shares within type
#   ex1-pure-heads-share.csv            share of Heads on each pure sequence by type
#
# NOTES
#   - Covariates are standardised within treatment over all participants with
#     a finite value; response time is the mean log screen time over the
#     participant's betting trials
#   - The mirror pairs are listed below; they are the four pairs with the
#     strongest side bias in RQ4 (main text, Sec. PsychMotives) and are
#     checked to be exact mirror images
#   - Output files are overwritten at every run
# ============================================================

ex1_pure_sequences <- function(cfg) {

  design <- cfg$design
  tr_vec <- unique(as.character(cfg$run$treatment))
  heads  <- as.character(design$seq$side_labels$heads)
  tails  <- as.character(design$seq$side_labels$tails)
  pure_H <- as.character(design$seq$anchor_labels$pure_heads)
  pure_T <- as.character(design$seq$anchor_labels$pure_tails)

  mirror_seq <- c("HHTHTT", "TTHTHH", "HTHHTT", "THTTHH",
                  "TTHHTT", "HHTTHH", "HHHTTT", "TTTHHH")
  flip <- function(x) chartr("HT", "TH", x)
  stopifnot(all(flip(mirror_seq[c(1, 3, 5, 7)]) == mirror_seq[c(2, 4, 6, 8)]))
  # gamblerish response = bet against the terminal run: Heads if the run is Tails
  mirror <- data.table(seq = mirror_seq,
                       gambler_h = as.integer(substring(mirror_seq, 6, 6) == tails))

  type_levels <- c("hothand_consistent", "gambler_consistent", "inconsistent", "unclassified")

  lotr <- fread(file.path(path_der, "lotr-scored.csv"), encoding = "UTF-8")
  lotr <- lotr[, .(pid = as.character(pid), lotr_score = as.numeric(lotr_score))]

  z_score <- function(x) {
    s <- sd(x, na.rm = TRUE)
    if (!is.finite(s) || s <= 0) return(rep(NA_real_, length(x)))
    (x - mean(x, na.rm = TRUE)) / s
  }
  ci_t <- function(x) {
    n <- length(x)
    if (n < 2L) return(c(NA_real_, NA_real_))
    mean(x) + c(-1, 1) * qt(0.975, n - 1) * sd(x) / sqrt(n)
  }

  counts_all <- list(); cov_all <- list(); xtab_all <- list(); ph_all <- list()

  for (tr in tr_vec) {

    d_all <- analysis_trials(tr, "rq1-full")
    pids  <- sort(unique(as.character(d_all$pid)))

    bets <- d_all[is.finite(stake) & stake > 0 & side %in% c(heads, tails)]
    bets[, `:=`(pid = as.character(pid), h = as.integer(side == heads))]

    # ---- Pure-sequence type (one trial per participant and sequence) ----
    hH <- bets[seq == pure_H, .(pid, h_H = h)]
    hT <- bets[seq == pure_T, .(pid, h_T = h)]
    part <- merge(merge(data.table(pid = pids), hH, by = "pid", all.x = TRUE),
                  hT, by = "pid", all.x = TRUE)
    part[, chi_pure := fcase(is.na(h_H) | is.na(h_T), NA_integer_,
                             h_H == 1L & h_T == 0L,  1L,
                             h_H == 0L & h_T == 1L, -1L,
                             default = 0L)]
    part[, type := fcase(is.na(chi_pure), "unclassified",
                         chi_pure ==  1L, "hothand_consistent",
                         chi_pure == -1L, "gambler_consistent",
                         default = "inconsistent")]
    part[, type := factor(type, levels = type_levels)]

    # ---- Behaviour on the mirror sequences ----
    m <- merge(bets[seq %in% mirror_seq, .(pid, seq, h)], mirror, by = "seq")
    m[, score := fifelse(h == gambler_h, 1L, -1L)]
    m_pid <- m[, .(n_mirror_bets = .N, mirror_score = mean(score)), by = pid]
    m_pid[, mirror_type := fcase(mirror_score > 0, "gambler", mirror_score < 0, "hothand",
                                 default = "split")]
    part <- merge(part, m_pid, by = "pid", all.x = TRUE)
    setorder(part, pid)

    f_part <- file.path(path_csv, paste0("ex1-", tr, "-pure-participants.csv"))
    fwrite(part, f_part);  msg("Saved: ", f_part)

    # ---- (1) counts ----
    cnt <- part[, .N, by = type][, treatment := tr]
    cnt[, share := N / sum(N)]
    counts_all[[tr]] <- cnt[order(type)]

    # ---- (2) covariates by type ----
    r_tbl <- fread(model_path("mpl", paste0("scored-", tr, ".csv")),
                   encoding = "UTF-8")[, .(pid = as.character(pid), r_mean)]
    rt <- bets[is.finite(screen_ms) & screen_ms > 0,
               .(rt_log_mean = mean(log(screen_ms))), by = pid]
    cov <- Reduce(function(a, b) merge(a, b, by = "pid", all.x = TRUE),
                  list(part[, .(pid, type)], lotr, r_tbl, rt))
    cov[, `:=`(Zopt = z_score(lotr_score), Zrt = z_score(rt_log_mean), Zr = z_score(r_mean))]
    cov_names <- c(Zopt = "Optimism (LOT-R)", Zrt = "Log response time", Zr = "Risk aversion (r)")
    cov_all[[tr]] <- rbindlist(lapply(names(cov_names), function(v) {
      rbindlist(lapply(type_levels, function(tp) {
        x <- cov[type == tp][[v]]
        x <- x[is.finite(x)]
        ci <- ci_t(x)
        data.table(treatment = tr, covariate = cov_names[[v]], type = tp, n = length(x),
                   mean = if (length(x)) mean(x) else NA_real_,
                   sd = if (length(x) > 1L) sd(x) else NA_real_,
                   ci_lo = ci[1], ci_hi = ci[2])
      }))
    }))

    # ---- (3) cross-tabulation and Heads share ----
    xt <- part[type %in% type_levels[1:3] & !is.na(mirror_type), .N, by = .(type, mirror_type)]
    xt[, share_within_type := N / sum(N), by = type]
    xt[, treatment := tr]
    xtab_all[[tr]] <- xt[order(type, mirror_type)]

    ph <- merge(bets[seq %in% c(pure_H, pure_T), .(pid, seq, h)],
                part[, .(pid, type)], by = "pid")
    ph_all[[tr]] <- ph[, .(n = .N, heads_share = mean(h),
                           ci_lo = ci_t(h)[1], ci_hi = ci_t(h)[2]),
                       by = .(seq, type)][, treatment := tr][order(seq, type)]
  }

  out <- list(`ex1-pure-counts.csv`     = counts_all,
              `ex1-pure-covariates.csv` = cov_all,
              `ex1-pure-vs-mirror.csv`  = xtab_all,
              `ex1-pure-heads-share.csv` = ph_all)
  for (nm in names(out)) {
    f <- file.path(path_csv, nm)
    fwrite(rbindlist(out[[nm]]), f);  msg("Saved: ", f)
  }

  invisible(TRUE)
}
