# ============================================================
# 99_synthesis_reaction_times.R
#
# PURPOSE
#   Part A -- RT BY SEQUENCE. Do the 12 balanced sequences without a
#   terminal run take longer to process, and do response times rise with
#   the number of alternations?
#   Part B -- RT AND PATTERN COMPLETION. Does longer deliberation go with
#   ACTING on the sequence's local pattern? On sequences that end in an
#   alternation (no terminal run), the pattern-completing bet is the side
#   OPPOSITE to the last symbol (THTHTH -> T, HTHTHT -> H, HHTTHT -> H,
#   TTHHTH -> T). If longer RTs reflect a search for structure, trials on
#   which a participant takes longer than usual should more often end in
#   the pattern-completing bet. Perceptual complexity alone (more switches
#   = longer encoding) does not make this prediction, which is what lets
#   Part B discriminate where Part A's alternation-RT correlation cannot.
#
# PART A -- method
#   Log RT is centred within participant (lrt_c); one difference per
#   participant: mean(lrt_c | balanced & no run) - mean(lrt_c | other 52).
#   exp(mean diff) - 1 = average proportional slowdown. Plus the
#   sequence-level Spearman correlation between alternations and mean lrt_c.
#
# PART B -- method
#   Trials: bets only (side H/T) on sequences WITHOUT a terminal run.
#     primary  : the 12 balanced no-run sequences
#     secondary: all 32 no-run sequences
#   Outcome : extend = 1 if the chosen side is opposite to the last symbol.
#   RT      : lrt_w = log RT centred within participant AND within
#             sequence (demeaned by participant, then by sequence), so
#             only a participant's unusual slowness on THAT sequence enters.
#   Models  :
#     (1) linear probability model extend ~ lrt_w + participant FE +
#         sequence FE, SEs clustered by participant (sandwich::vcovCL);
#         coefficient = change in P(extend) per 1 log-unit of RT;
#         also reported per +50% RT (log(1.5) units)
#     (2) logistic GLMM extend ~ lrt_w + (1|pid) + (1|seq) (lme4),
#         as a robustness check on functional form
#   Descriptive: P(extend) by within-participant RT tercile.
#   Prediction under pattern-seeking: coefficient on lrt_w > 0.
#
# DEFINITIONS (from the 6-symbol sequence string itself)
#   balanced  : 3 H and 3 T
#   run       : last two symbols identical
#   bal_norun : balanced AND no terminal run (12 sequences)
#   alt       : number of switches between adjacent symbols (0-5)
#   extend    : chosen side != last symbol (defined for no-run sequences)
#
# INPUT
#   path_src/master_sequences.csv  (pid, treat, seq, screen_ms, side)
#     side: "H", "T", or "NB" (no bet)
#
# OUTPUT
#   path_out/synthesis_reaction_times_summary.csv       -- Part A, one row per treatment
#   path_out/synthesis_reaction_times_sequences.csv     -- Part A, one row per treatment x sequence
#   path_out/synthesis_reaction_times_pattern_test.csv  -- Part B, model coefficients
#   path_out/synthesis_reaction_times_pattern_terciles.csv -- Part B, P(extend) by RT tercile
#
# NOTES
#   - Trials with missing, non-positive, or non-finite screen_ms are
#     dropped. RTs are analysed on the log scale (no trimming rule).
#   - Part A uses all trials (bet and no-bet); Part B uses bets only.
#   - Requires sandwich and lmtest (Part B, model 1) and lme4 (model 2);
#     model 2 is skipped with an NA row if lme4 fails to fit.
#   - Only console output: "Saved: <path>".
# ============================================================
synthesis_reaction_times <- function(cfg) {
  
  for (pkg in c("sandwich", "lmtest", "lme4")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("synthesis_reaction_times: package '", pkg, "' is required.")
    }
  }
  
  tr_vec <- unique(as.character(cfg$run$treatment))
  
  f_in <- file.path(path_src, "master_sequences.csv")
  stopifnot(file.exists(f_in))
  m <- fread(f_in, encoding = "UTF-8")
  stopifnot(all(c("pid", "treat", "seq", "screen_ms", "side") %in% names(m)))
  
  m[, `:=`(pid = as.character(pid), treat = as.character(treat), seq = as.character(seq),
           side = as.character(side), rt_s = as.numeric(screen_ms) / 1000)]
  m <- m[is.finite(rt_s) & rt_s > 0]
  m[, lrt := log(rt_s)]
  
  m[, nH        := nchar(gsub("T", "", seq))]
  m[, last      := substr(seq, 6, 6)]
  m[, run       := substr(seq, 5, 5) == last]
  m[, alt       := vapply(strsplit(seq, ""), function(x) sum(x[-1] != x[-length(x)]), integer(1))]
  m[, bal_norun := nH == 3L & !run]
  m[, lrt_c     := lrt - mean(lrt), by = .(treat, pid)]    # within-participant centring
  
  summ <- list(); seqs <- list(); pat <- list(); terc <- list()
  
  for (tr in tr_vec) {
    d <- m[treat == tr]
    if (nrow(d) == 0L) next
    
    # ======================================================
    # PART A -- RT by sequence
    # ======================================================
    pd <- d[, .(diff = mean(lrt_c[bal_norun]) - mean(lrt_c[!bal_norun])), by = pid]
    pd <- pd[is.finite(diff)]
    
    sq <- d[, .(lrt_c = mean(lrt_c), raw_mean_rt_s = mean(rt_s), n_trials = .N,
                alt = alt[1], balanced = nH[1] == 3L, run = run[1], bal_norun = bal_norun[1]),
            by = seq]
    sq[, treatment := tr]
    setorder(sq, -lrt_c)
    sq[, rt_rank := seq_len(.N)]                            # 1 = slowest sequence
    
    ct <- suppressWarnings(cor.test(sq$alt, sq$lrt_c, method = "spearman", exact = FALSE))
    tt <- t.test(pd$diff)
    wt <- suppressWarnings(wilcox.test(pd$diff))
    
    summ[[tr]] <- data.table(
      treatment                 = tr,
      n_participants            = nrow(pd),
      pct_slower_bal_norun      = 100 * (exp(mean(pd$diff)) - 1),
      pct_slower_ci_low         = 100 * (exp(tt$conf.int[1]) - 1),
      pct_slower_ci_high        = 100 * (exp(tt$conf.int[2]) - 1),
      t_p                       = tt$p.value,
      wilcox_p                  = wt$p.value,
      share_participants_slower = mean(pd$diff > 0),
      mean_rank_bal_norun       = mean(sq[bal_norun == TRUE, rt_rank]),   # out of 64; 1 = slowest
      mean_alt_bal_norun        = mean(sq[bal_norun == TRUE,  alt]),
      mean_alt_other            = mean(sq[bal_norun == FALSE, alt]),
      share_alt_ge4_that_are_bal_norun = mean(sq[alt >= 4, bal_norun]),  # combinatorial, not data
      rho_alt_rt                = unname(ct$estimate),
      rho_alt_rt_p              = ct$p.value
    )
    seqs[[tr]] <- sq
    
    # ======================================================
    # PART B -- RT and pattern completion
    # ======================================================
    b_all <- d[!run & side %in% c("H", "T")]
    b_all[, extend := as.integer(side != last)]
    
    for (subset in c("balanced no-run (12)", "all no-run (32)")) {
      b <- if (subset == "balanced no-run (12)") b_all[bal_norun == TRUE] else copy(b_all)
      if (nrow(b) < 50L) next
      # within participant, then within sequence
      b[, lrt_w := lrt - mean(lrt), by = pid]
      b[, lrt_w := lrt_w - mean(lrt_w), by = seq]
      
      # (1) LPM with participant + sequence FE, clustered SE
      lpm <- lm(extend ~ lrt_w + factor(pid) + factor(seq), data = b)
      V   <- sandwich::vcovCL(lpm, cluster = ~ pid)
      ct1 <- lmtest::coeftest(lpm, vcov. = V)["lrt_w", ]
      ci1 <- lmtest::coefci(lpm, vcov. = V)["lrt_w", ]
      
      # (2) logistic GLMM
      g <- tryCatch(
        suppressWarnings(suppressMessages(   # console stays clean; singularity is recorded below
          lme4::glmer(extend ~ lrt_w + (1 | pid) + (1 | seq), data = b, family = binomial,
                      control = lme4::glmerControl(optimizer = "bobyqa")))),
        error = function(e) NULL)
      ct2 <- if (!is.null(g)) summary(g)$coefficients["lrt_w", ] else rep(NA_real_, 4)
      
      pat[[length(pat) + 1]] <- data.table(
        treatment = tr, subset = subset,
        n_trials = nrow(b), n_participants = uniqueN(b$pid), n_sequences = uniqueN(b$seq),
        share_extend = mean(b$extend),
        lpm_beta_per_logrt = ct1[1], lpm_se = ct1[2], lpm_p = ct1[4],
        lpm_ci_low = ci1[1], lpm_ci_high = ci1[2],
        lpm_pp_per_50pct_slower = 100 * ct1[1] * log(1.5),
        glmm_logodds_per_logrt = ct2[1], glmm_se = ct2[2], glmm_p = ct2[4],
        glmm_singular = if (!is.null(g)) lme4::isSingular(g) else NA
      )
      
      b[, rt_tercile := cut(lrt_w, quantile(lrt_w, c(0, 1/3, 2/3, 1)),
                            labels = c("fast", "middle", "slow"), include.lowest = TRUE)]
      terc[[length(terc) + 1]] <- b[, .(n_trials = .N, share_extend = mean(extend)),
                                    by = rt_tercile][order(rt_tercile)][
                                      , `:=`(treatment = tr, subset = subset)]
    }
  }
  
  save <- function(x, name, label) {
    f <- file.path(path_out, name)
    if (!should_skip(f, cfg, "output", label)) { fwrite(rbindlist(x), f); msg("Saved: ", f) }
  }
  save(summ, "synthesis_reaction_times_summary.csv",         "RT by sequence summary")
  save(seqs, "synthesis_reaction_times_sequences.csv",       "RT by sequence, per sequence")
  save(pat,  "synthesis_reaction_times_pattern_test.csv",    "RT and pattern completion")
  save(terc, "synthesis_reaction_times_pattern_terciles.csv","RT terciles and pattern completion")
  
  invisible(list(summary = rbindlist(summ), sequences = rbindlist(seqs),
                 pattern_test = rbindlist(pat), pattern_terciles = rbindlist(terc)))
}