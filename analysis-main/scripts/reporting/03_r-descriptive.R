# ============================================================
# 03_r-descriptive.R
#
# PURPOSE
#   Descriptive analysis of the CRRA parameter r_i (posterior mean of the
#   primary individual-level Holt-Laury model) in two parts.
#   Part 1 (r_descriptive_summary) summarises the distribution of r_i across
#   participants: for each treatment, for all, consistent and inconsistent
#   participants, and within each group for all participants, women and men,
#   it reports the number of participants, the minimum, quartiles, maximum,
#   mean and standard deviation.
#   Part 2 (r_consistency_sex) relates the Holt-Laury consistency and r_i to
#   sex, separately for each treatment. It reports, for women and men, the
#   number and share of participants with an inconsistent choice pattern with
#   Fisher's exact test, and compares the mean r_i of consistent women and men
#   (difference of means with a Welch 95% confidence interval and a Wilcoxon
#   rank-sum test). It also describes the sex composition of the full sample
#   of each treatment and of its consistent participants, which is the
#   confirmatory sample of the FN treatment, and tests whether the sex
#   difference in inconsistency depends on the treatment.
#   The figures of the distribution are produced by 04_r-distribution.R, the
#   comparison of estimates and treatments by 05_r-robustness.R.
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
#   r-consistency-by-sex.csv    one row per treatment: group sizes, numbers and
#                               shares of inconsistent participants by sex,
#                               Fisher's exact p-value and conditional odds
#                               ratio, mean r_i of consistent women and men,
#                               their difference with a 95% CI, and the
#                               Wilcoxon rank-sum p-value
#   r-sample-composition.csv    one row per treatment and sample (all,
#                               consistent): numbers of women, men and
#                               participants of other or unreported sex, and
#                               the ratio of women to men
#   r-sex-treatment-interaction.csv
#                               one row per test: (i) the interaction of sex
#                               and treatment in a logistic regression of
#                               inconsistency, reported as the ratio of the
#                               odds ratios (women against men) of the second
#                               and the first treatment with a Wald 95% CI and
#                               the p-value of a likelihood-ratio test; (ii)
#                               the difference between the treatments in the
#                               inconsistency of women and of men, and (iii)
#                               in the sex composition (Fisher's exact tests
#                               with conditional odds ratios)
#
# NOTES
#   - The sex group "All" in the summary includes participants who do not
#     report their sex (four in the FN treatment), so the rows for women and
#     men need not add up to it
#   - Standard deviations are across participants (not posterior SDs) and are
#     reported only for groups with at least two participants
#   - Estimates of inconsistent participants are bounded above by the
#     indifference point of the last row of the list (r = 1.70), so their
#     distribution accumulates near this ceiling
#   - Tests of the sex differences use participants who report their sex as
#     female or male; the remaining participants appear only in the
#     composition table
#   - Posterior means of r_i are treated as data in the comparison of consistent
#     women and men; their posterior SD is small (see 05_r-robustness.R)
#   - Both tests are two-sided
#   - The interaction test is exploratory: it is not preregistered, and the
#     cells of the second treatment are small
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

r_consistency_sex <- function(cfg) {
  
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
  dt[, sex_label := fcase(sex == "F", "Women", sex == "M", "Men", default = "Other")]
  
  # ---- Inconsistency and risk aversion by sex ----
  res_sex <- rbindlist(lapply(tr_vec, function(tr) {
    d <- dt[treat == tr & sex_label %in% c("Women", "Men")]
    
    # rows: sex; columns: inconsistent (1) and consistent (0)
    tab <- table(factor(d$sex_label,    levels = c("Women", "Men")),
                 factor(d$inconsistent, levels = c(1L, 0L)))
    ft  <- fisher.test(tab)
    
    r_w <- d[sex_label == "Women" & inconsistent == 0L, r_mean]
    r_m <- d[sex_label == "Men"   & inconsistent == 0L, r_mean]
    tt  <- t.test(r_w, r_m)                    # Welch
    wt  <- wilcox.test(r_w, r_m, exact = FALSE)
    
    data.table(
      treatment               = tr,
      n_women                 = sum(tab["Women", ]),
      n_men                   = sum(tab["Men", ]),
      n_inconsistent_women    = tab["Women", "1"],
      n_inconsistent_men      = tab["Men", "1"],
      share_inconsistent_women = tab["Women", "1"] / sum(tab["Women", ]),
      share_inconsistent_men   = tab["Men", "1"]   / sum(tab["Men", ]),
      fisher_p                = ft$p.value,
      odds_ratio              = unname(ft$estimate),
      odds_ratio_lo           = ft$conf.int[1],
      odds_ratio_hi           = ft$conf.int[2],
      n_consistent_women      = length(r_w),
      n_consistent_men        = length(r_m),
      mean_r_consistent_women = mean(r_w),
      mean_r_consistent_men   = mean(r_m),
      diff_r_women_minus_men  = mean(r_w) - mean(r_m),
      diff_r_lo               = tt$conf.int[1],
      diff_r_hi               = tt$conf.int[2],
      wilcoxon_p              = wt$p.value
    )
  }))
  
  f_sex <- file.path(path_csv, "r-consistency-by-sex.csv")
  fwrite(res_sex, f_sex)
  msg("Saved: ", f_sex)
  
  # ---- Sex composition of the full and the consistent sample ----
  res_comp <- rbindlist(lapply(tr_vec, function(tr) {
    rbindlist(lapply(c("all", "consistent"), function(smp) {
      d <- dt[treat == tr]
      if (smp == "consistent") d <- d[inconsistent == 0L]
      data.table(
        treatment = tr, sample = smp, n = nrow(d),
        n_women = sum(d$sex_label == "Women"),
        n_men   = sum(d$sex_label == "Men"),
        n_other = sum(d$sex_label == "Other")
      )
    }))
  }))
  res_comp[, women_per_man := n_women / n_men]
  
  f_comp <- file.path(path_csv, "r-sample-composition.csv")
  fwrite(res_comp, f_comp)
  msg("Saved: ", f_comp)
  
  # ---- Does the sex difference in inconsistency depend on the treatment? ----
  stopifnot(length(tr_vec) == 2L)
  d_int <- dt[sex_label %in% c("Women", "Men")]
  d_int[, female := as.integer(sex_label == "Women")]
  d_int[, tr     := factor(treat, levels = tr_vec)]
  
  fit0 <- glm(inconsistent ~ female + tr, family = binomial, data = d_int)
  fit1 <- glm(inconsistent ~ female * tr, family = binomial, data = d_int)
  lrt  <- anova(fit0, fit1, test = "LRT")
  term <- paste0("female:tr", tr_vec[2])
  cf   <- summary(fit1)$coefficients[term, ]
  
  one_row <- function(test, estimate, lo, hi, p, method, note) {
    data.table(test = test, estimate = estimate, lo = lo, hi = hi, p_value = p,
               method = method, note = note)
  }
  
  # conditional odds ratio of a 2x2 table with Fisher's exact test
  fisher_row <- function(test, tab, note) {
    ft <- fisher.test(tab)
    one_row(test, unname(ft$estimate), ft$conf.int[1], ft$conf.int[2], ft$p.value,
            "Fisher's exact", note)
  }
  
  tr_tab <- function(d, col, levels_col)
    table(factor(d$treat, levels = tr_vec), factor(d[[col]], levels = levels_col))
  
  res_int <- rbindlist(list(
    one_row("sex_by_treatment_interaction_inconsistency",
            exp(cf["Estimate"]),
            exp(cf["Estimate"] - qnorm(0.975) * cf["Std. Error"]),
            exp(cf["Estimate"] + qnorm(0.975) * cf["Std. Error"]),
            lrt$`Pr(>Chi)`[2],
            "logistic regression, likelihood-ratio test",
            paste0("ratio of the odds ratios (women against men) of ", tr_vec[2], " and ", tr_vec[1])),
    fisher_row("treatment_difference_inconsistency_women",
               tr_tab(d_int[female == 1L], "inconsistent", c(1L, 0L)),
               paste0("odds ratio of inconsistency, ", tr_vec[1], " against ", tr_vec[2])),
    fisher_row("treatment_difference_inconsistency_men",
               tr_tab(d_int[female == 0L], "inconsistent", c(1L, 0L)),
               paste0("odds ratio of inconsistency, ", tr_vec[1], " against ", tr_vec[2])),
    fisher_row("treatment_difference_sex_composition",
               tr_tab(d_int, "female", c(1L, 0L)),
               paste0("odds ratio of being a woman, ", tr_vec[1], " against ", tr_vec[2]))
  ))
  
  f_int <- file.path(path_csv, "r-sex-treatment-interaction.csv")
  fwrite(res_int, f_int)
  msg("Saved: ", f_int)
  
  invisible(list(by_sex = res_sex, composition = res_comp, interaction = res_int))
}

# ---- Both parts ----
r_descriptive <- function(cfg) {
  r_descriptive_summary(cfg)
  r_consistency_sex(cfg)
  invisible(NULL)
}
