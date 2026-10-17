# ============================================================
# 96_synthesis-debrief-validation.R
#
# PURPOSE
#   Validation of the model-inferred side-choice heterogeneity classes
#   (with, against or uncertain, for the imbalance and the run cue) against
#   the debriefing questionnaire, in three descriptive passes:
#     (1) used_strategy (yes / no / unsure) by class, all participants, with
#         a Fisher test (no or unsure vs. yes, by uncertain vs. directional)
#         and an omnibus chi-square test, for both cues.
#     (2) The hand-coded free-text strategy category (anti-trend, pro-trend,
#         no-trend, unclear, misc) by class, restricted to participants who
#         answered "yes" and left a codable description. Three directional
#         Fisher tests: anti-trend by against_imb, pro-trend by with_imb,
#         no-trend by uncertain.
#     (3) Likert items by class, all participants: Kruskal-Wallis and
#         Mann-Whitney tests for belief_independence and
#         reliance_on_sequence, and pairwise class contrasts (mean
#         difference with a percentile-bootstrap 95% interval).
#   The passes are descriptive. The coding scheme of pass (2) was fixed after
#   the classification existed, so the p-values are not confirmatory.
#
# SAMPLE
#   Full sample of each treatment in cfg$run$treatment (the sample of the
#   joint model "synthesis-joint").
#
# INPUT
#   data/output/csv/synthesis-<tr>-joint-{imb,run}-slope-participants.csv
#                                  (95_synthesis-joint-heterogeneity.R)
#   data/source/debriefing.csv     questionnaire, one row per participant
#   raw-data/manual/manual_debrief_strategies.csv
#                                  hand-coded strategies, semicolon-delimited,
#                                  columns participant.code, treatment,
#                                  strategy_text, category (coded blind to
#                                  the model-derived classes)
#
# OUTPUT (data/output/csv/)
#   synthesis-debrief-used-strategy-crosstab-{imb,run}.csv
#   synthesis-debrief-used-strategy-tests.csv
#   synthesis-debrief-used-strategy-coverage.csv
#   synthesis-debrief-merged.csv
#   synthesis-debrief-crosstab-{imb,run}-<tr>.csv, ...-all.csv
#   synthesis-debrief-fisher-tests.csv
#   synthesis-debrief-likert-by-class.csv
#   synthesis-debrief-likert-tests.csv
#   synthesis-debrief-likert-contrasts.csv
#   synthesis-debrief-likert-by-category.csv
#
# NOTES
#   - Pass (2) is restricted to participants who answered "yes" and left a
#     codable description, so that it answers whether what participants
#     described tracks the class among those who reported a strategy.
#     Participants who answered "no" or "unsure" but left codable text enter
#     pass (1) only; the overlap is in the coverage file
#   - The pass (2) tests are run on the pooled sample only, since the cells
#     are too thin per treatment; per-treatment crosstabs are still written
#   - In the Likert tests mw_p compares the collapsed directional group with
#     the uncertain group; single-class contrasts are in the contrasts file
#   - The run-cue class is cross-tabulated but not tested directionally,
#     because with_run is nearly empty
#   - The alternative of the Fisher tests, the coverage warning threshold and
#     the number of bootstrap draws are set in settings.R
#     (cfg$synthesis$debrief); the bootstrap is seeded by cfg$run$seed
#   - Output files are overwritten at every run
# ============================================================

synthesis_debrief_validation <- function(cfg) {
  
  tr_vec <- unique(as.character(cfg$run$treatment))
  
  cat_levels <- c("anti-trend", "pro-trend", "no-trend", "unclear", "misc")
  imb_levels <- c("with_imb", "uncertain", "against_imb")
  run_levels <- c("with_run", "uncertain", "against_run")
  
  alternative  <- cfg$synthesis$debrief$alternative
  min_coverage <- cfg$synthesis$debrief$min_coverage
  
  # ============================================================
  # STEP 0 -- model-derived per-participant classes, every treatment,
  # read once and reused by all three passes below.
  # ============================================================
  full_universe_list <- list()
  for (tr in tr_vec) {
    
    f_imb <- file.path(path_csv, paste0("synthesis-", tr, "-joint-imb-slope-participants.csv"))
    f_run <- file.path(path_csv, paste0("synthesis-", tr, "-joint-run-slope-participants.csv"))
    if (!file.exists(f_imb) || !file.exists(f_run))
      stop("Missing ", f_imb, " and/or ", f_run, " (run 95_synthesis-joint-heterogeneity.R first).", call. = FALSE)

    imb_all <- fread(f_imb)[, .(pid, imb_median = median, imb_class = class)]
    run_all <- fread(f_run)[, .(pid, run_median = median, run_class = class)]
    stopifnot(setequal(imb_all$pid, run_all$pid))  # same participant universe, same joint model
    
    # Participants are matched by name with merge()
    full_tr <- merge(imb_all, run_all, by = "pid")
    full_tr[, treatment := tr]
    full_universe_list[[tr]] <- full_tr
  }
  if (!length(full_universe_list)) stop("No model-derived classes found.", call. = FALSE)
  full_universe <- rbindlist(full_universe_list)  # pid, imb_*, run_*, treatment
  
  # ============================================================
  # PASS (1) -- used_strategy (yes/no/unsure) x class, all participants
  # ============================================================
  f_full <- file.path(path_src, "debriefing.csv")
  if (!file.exists(f_full)) stop("Debriefing questionnaire not found: ", f_full)
  
  debrief_q <- fread(f_full, encoding = "UTF-8")
  if (anyDuplicated(debrief_q$pid)) stop("Duplicate pid(s) in ", f_full)
  debrief_q[, used_strategy := tolower(trimws(used_strategy))]
  bad_used <- setdiff(unique(debrief_q$used_strategy), c("yes", "no", "unsure"))
  if (length(bad_used)) {
    stop("Unexpected used_strategy value(s) in ", f_full, ": ", paste(bad_used, collapse = ", "))
  }
  
  # LEFT join off full_universe (the participants with a fitted
  # personal slope) so a questionnaire row for someone outside that set
  # -- or a model participant missing a questionnaire row -- shows up as
  # an explicit NA/warning rather than silently vanishing.
  full_q <- merge(full_universe, debrief_q, by = "pid", all.x = TRUE)
  no_q <- full_q[is.na(used_strategy), pid]
  if (length(no_q)) {
    warning("Model participant(s) with no row in ", f_full, ": ", paste(no_q, collapse = ", "))
  }
  orphan_q <- setdiff(debrief_q$pid, full_universe$pid)
  if (length(orphan_q)) {
    warning("Questionnaire pid(s) with no model-derived class: ", paste(orphan_q, collapse = ", "))
  }
  full_q <- full_q[!is.na(used_strategy)]
  
  ct_used <- function(class_col, class_levels) {
    tt <- table(factor(full_q$used_strategy, levels = c("yes", "no", "unsure")),
                factor(full_q[[class_col]], levels = class_levels))
    out <- as.data.table(tt)
    setnames(out, c("used_strategy", "class", "n"))
    out <- dcast(out, used_strategy ~ class, value.var = "n")
    setcolorder(out, c("used_strategy", class_levels))
    out[, N := rowSums(.SD), .SDcols = class_levels]
    out
  }
  used_imb_tbl <- ct_used("imb_class", imb_levels)
  used_run_tbl <- ct_used("run_class", run_levels)
  f_used_imb <- file.path(path_csv, "synthesis-debrief-used-strategy-crosstab-imb.csv")
  f_used_run <- file.path(path_csv, "synthesis-debrief-used-strategy-crosstab-run.csv")
  fwrite(used_imb_tbl, f_used_imb)
  fwrite(used_run_tbl, f_used_run)
  msg("Saved: ", f_used_imb)
  msg("Saved: ", f_used_run)
  
  # "no/unsure participants are mostly uncertain" -- tested, not assumed:
  # 2x2 exact test (no/unsure vs. yes) x (uncertain vs. not), plus the
  # omnibus 3x3 chi-square across all three used_strategy levels, for
  # both classes. 
  used_test <- function(class_col, class_levels) {
    no_or_unsure <- factor(full_q$used_strategy != "yes", levels = c(FALSE, TRUE))
    is_uncertain <- factor(full_q[[class_col]] == "uncertain", levels = c(FALSE, TRUE))
    ft <- fisher.test(table(no_or_unsure, is_uncertain))
    ct3 <- table(factor(full_q$used_strategy, levels = c("yes", "no", "unsure")),
                 factor(full_q[[class_col]], levels = class_levels))
    chi <- suppressWarnings(chisq.test(ct3))
    data.table(class_dim = class_col,
               odds_ratio_2x2 = unname(ft$estimate), p_2x2 = ft$p.value,
               chisq_3x3 = unname(chi$statistic), df_3x3 = unname(chi$parameter), p_3x3 = chi$p.value)
  }
  used_test_tbl <- rbindlist(list(
    used_test("imb_class", imb_levels),
    used_test("run_class", run_levels)
  ))
  f_used_tests <- file.path(path_csv, "synthesis-debrief-used-strategy-tests.csv")
  fwrite(used_test_tbl, f_used_tests)
  msg("Saved: ", f_used_tests)
  
  yes_pid <- full_q[used_strategy == "yes", pid]
  
  # ============================================================
  # PASS (2) -- free-text category x class, restricted to
  # used_strategy == "yes" AND codable (see RESTRICTION in the header)
  # ============================================================
  f_debrief <- file.path(path_manual, "manual_debrief_strategies.csv")
  if (!file.exists(f_debrief)) stop("Debrief file not found: ", f_debrief)
  
  debrief <- fread(f_debrief, sep = ";", encoding = "UTF-8")
  setnames(debrief, old = "participant.code", new = "pid")
  
  # The raw category column has inconsistent spacing/hyphenation from
  # hand-coding ("no trend" vs "no-trend", trailing whitespace on
  # "unclear "). Normalized to a canonical hyphenated lowercase form
  # rather than trusted verbatim -- fread would otherwise silently
  # split what should be one category into two distinct factor levels.
  debrief[, category := tolower(trimws(category))]
  debrief[, category := gsub("\\s+", "-", category)]
  
  bad_cat <- setdiff(unique(debrief$category), cat_levels)
  if (length(bad_cat)) {
    stop(
      "Unexpected debrief category label(s) after normalization: ",
      paste(bad_cat, collapse = ", "),
      ". Expected one of: ", paste(cat_levels, collapse = ", "),
      ". Fix the raw label(s) in ", f_debrief, " (participant(s): ",
      paste(debrief[category %in% bad_cat, pid], collapse = ", "), ")."
    )
  }
  debrief[, category := factor(category, levels = cat_levels)]
  
  dup_pid <- debrief[, .N, by = .(pid, treatment)][N > 1, pid]
  if (length(dup_pid)) {
    stop("Duplicate (pid, treatment) rows in debrief file: ", paste(dup_pid, collapse = ", "))
  }
  
  # yes/codable overlap, logged as a file (all codable pids, including
  # the two who answered "no"/"unsure" and are excluded from the
  # crosstabs below) rather than assumed to be a perfect match.
  coded_pid     <- debrief$pid
  yes_not_coded <- setdiff(yes_pid, coded_pid)
  coded_not_yes <- setdiff(coded_pid, yes_pid)
  coverage_tbl <- data.table(
    n_used_strategy_yes = length(yes_pid),
    n_codable_category  = length(coded_pid),
    n_yes_and_coded     = length(intersect(yes_pid, coded_pid)),
    n_yes_not_coded     = length(yes_not_coded),
    n_coded_not_yes     = length(coded_not_yes),
    coded_not_yes_pids  = paste(coded_not_yes, collapse = "; ")
  )
  f_used_cov <- file.path(path_csv, "synthesis-debrief-used-strategy-coverage.csv")
  fwrite(coverage_tbl, f_used_cov)
  msg("Saved: ", f_used_cov)
  
  build_crosstab <- function(dt, class_col, class_levels) {
    tt <- table(
      factor(dt$category, levels = cat_levels),
      factor(dt[[class_col]], levels = class_levels)
    )
    out <- as.data.table(tt)
    setnames(out, c("category", "class", "n"))
    out <- dcast(out, category ~ class, value.var = "n")
    setcolorder(out, c("category", class_levels))
    out[, N := rowSums(.SD), .SDcols = class_levels]
    out
  }
  
  merged_list <- list()
  crosstabs_imb <- list()
  crosstabs_run <- list()
  
  for (tr in tr_vec) {
    
    full_tr <- full_universe_list[[tr]]
    if (is.null(full_tr)) next  # already msg()'d as skipped in STEP 0
    
    deb_tr <- debrief[treatment == tr & pid %in% yes_pid]  # <-- pass-(2) restriction
    orphan <- setdiff(deb_tr$pid, full_tr$pid)
    if (length(orphan)) {
      stop("Debrief pid(s) not found among [", tr, "] model participants: ",
           paste(orphan, collapse = ", "))
    }
    deb_cols <- deb_tr[, .(pid, strategy_text, category)]  # drop treatment: already known (== tr), and
    # keeping it would collide with full_tr's own
    # treatment column in the merge below
    
    m <- merge(full_tr, deb_cols, by = "pid", all.x = TRUE)  # LEFT join: keep every participant,
    # NA-fill strategy_text/category for the
    # ones with no codable "yes" free text
    n_total <- nrow(full_tr)
    n_coded <- sum(!is.na(m$category))
    coverage <- n_coded / n_total
    
    if (coverage < min_coverage) {
      warning("Debrief coverage [", tr, "] (used_strategy == 'yes' & codable) is ",
              round(100 * coverage, 1), "%, below the ", round(100 * min_coverage),
              "% sanity threshold (cfg$synthesis$debrief$min_coverage) -- check the ",
              "join keys and the raw debrief file before trusting the crosstabs below.")
    }
    
    m_coded <- m[!is.na(category)]
    merged_list[[tr]] <- m_coded[, .(pid, treatment = tr, category, strategy_text,
                                     imb_class, imb_median, run_class, run_median)]
    
    crosstabs_imb[[tr]] <- build_crosstab(m_coded, "imb_class", imb_levels)
    crosstabs_run[[tr]] <- build_crosstab(m_coded, "run_class", run_levels)
    
    f_ct_imb_tr <- file.path(path_csv, paste0("synthesis-debrief-crosstab-imb-", tr, ".csv"))
    f_ct_run_tr <- file.path(path_csv, paste0("synthesis-debrief-crosstab-run-", tr, ".csv"))
    fwrite(crosstabs_imb[[tr]], f_ct_imb_tr)
    fwrite(crosstabs_run[[tr]], f_ct_run_tr)
    msg("Saved: ", f_ct_imb_tr)
    msg("Saved: ", f_ct_run_tr)
  }
  
  merged_all <- rbindlist(merged_list)
  f_merged <- file.path(path_csv, "synthesis-debrief-merged.csv")
  fwrite(merged_all, f_merged)
  msg("Saved: ", f_merged)
  
  crosstab_imb_all <- build_crosstab(merged_all, "imb_class", imb_levels)
  crosstab_run_all <- build_crosstab(merged_all, "run_class", run_levels)
  f_ct_imb_all <- file.path(path_csv, "synthesis-debrief-crosstab-imb-all.csv")
  f_ct_run_all <- file.path(path_csv, "synthesis-debrief-crosstab-run-all.csv")
  fwrite(crosstab_imb_all, f_ct_imb_all)
  fwrite(crosstab_run_all, f_ct_run_all)
  msg("Saved: ", f_ct_imb_all)
  msg("Saved: ", f_ct_run_all)
  
  # Directional Fisher tests (pooled -- see NOTES in the header). Each row:
  # is the named category enriched for the named class, against
  # everyone else in the pooled yes-and-codable sample. 2x2 exact test
  # (not chi-square, because several cells are small),
  # `alternative` defaults to two-sided (matches how these numbers were
  # the default); set
  # cfg$synthesis$debrief$alternative = "greater" for the
  # one-sided version of the same a-priori-directional hypotheses.
  # The anti-trend and pro-trend categories are small, so point
  # estimates and intervals should be read alongside n_in_category.
  fisher_row <- function(cat_val, class_col, class_val) {
    # factor(..., levels = c(FALSE, TRUE)) forces table() to emit a full
    # 2x2 even when one combination is empty (e.g. no-trend x with_imb =
    # 0/12) -- cleaner and less error-prone than padding a table() result
    # after the fact, which silently mis-labels rows/columns if rbind/cbind
    # don't line up with what fisher.test() expects.
    in_cat   <- factor(merged_all$category == cat_val, levels = c(FALSE, TRUE))
    in_class <- factor(merged_all[[class_col]] == class_val, levels = c(FALSE, TRUE))
    tt <- table(in_cat, in_class)
    ft <- fisher.test(tt, alternative = alternative)
    data.table(
      category = cat_val, class_dim = class_col, target_class = class_val,
      n_in_category = sum(in_cat == TRUE),
      n_in_category_and_class = tt["TRUE", "TRUE"],
      n_outside_category_and_class = tt["FALSE", "TRUE"],
      odds_ratio = unname(ft$estimate), p_value = ft$p.value,
      ci_low = ft$conf.int[1], ci_high = ft$conf.int[2], alternative = alternative
    )
  }
  
  fisher_tbl <- rbindlist(list(
    fisher_row("anti-trend", "imb_class", "against_imb"),
    fisher_row("pro-trend",  "imb_class", "with_imb"),
    fisher_row("no-trend",   "imb_class", "uncertain")
  ))
  f_fisher <- file.path(path_csv, "synthesis-debrief-fisher-tests.csv")
  fwrite(fisher_tbl, f_fisher)
  msg("Saved: ", f_fisher)
  
  # ============================================================
  # PASS (3) -- Likert items x class, all participants (uses full_q from pass (1)
  # directly -- independent of the used_strategy/category restriction
  # in pass (2))
  # ============================================================
  likert_items <- c("belief_independence", "reliance_on_sequence", "action_seeking",
                    "self_risk_tolerance", "perceived_win_probability", "perceived_risk")
  
  # long format via melt -- avoids get()/dynamic-column NSE patterns entirely
  long_imb <- melt(full_q[, c("pid", "imb_class", likert_items), with = FALSE],
                   id.vars = c("pid", "imb_class"), variable.name = "item", value.name = "value")
  setnames(long_imb, "imb_class", "class")
  long_imb[, class_dim := "imb_class"]
  
  long_run <- melt(full_q[, c("pid", "run_class", likert_items), with = FALSE],
                   id.vars = c("pid", "run_class"), variable.name = "item", value.name = "value")
  setnames(long_run, "run_class", "class")
  long_run[, class_dim := "run_class"]
  
  long_all <- rbindlist(list(long_imb, long_run))
  likert_by_class <- long_all[, .(mean = mean(value), median = median(as.numeric(value)), n = .N),
                              by = .(class_dim, item, class)]
  setorder(likert_by_class, class_dim, item)
  f_likert_class <- file.path(path_csv, "synthesis-debrief-likert-by-class.csv")
  fwrite(likert_by_class, f_likert_class)
  msg("Saved: ", f_likert_class)
  
  # Kruskal-Wallis (3 groups) + Mann-Whitney (directional vs. uncertain,
  # collapsed) for the two theory-relevant items. The collapsed test is the more direct check of "does
  # having ANY credible personal slope track the item", but the omnibus
  # 3-group test is what actually motivated collapsing in the first
  # place, and both are kept because for belief_independence x
  # imb_class they tell different parts of the story -- with_imb is the
  # low outlier on both items, while against_imb sits close to
  # "uncertain" on belief_independence despite being highest of the
  # three on reliance_on_sequence.
  # NB: mw_p is the COLLAPSED directional-vs-uncertain test; single-class
  # contrasts with CIs are in synthesis-debrief-likert-contrasts.csv.
  likert_tests <- rbindlist(lapply(c("belief_independence", "reliance_on_sequence"), function(item) {
    rbindlist(lapply(list(list(cc = "imb_class", lv = imb_levels),
                          list(cc = "run_class", lv = run_levels)), function(spec) {
                            x  <- full_q[[item]]
                            g  <- full_q[[spec$cc]]
                            kw <- kruskal.test(x, g)
                            directional <- g %in% setdiff(spec$lv, "uncertain")
                            wt <- wilcox.test(x[directional], x[!directional])
                            data.table(item = item, class_dim = spec$cc,
                                       kw_stat = unname(kw$statistic), kw_df = unname(kw$parameter), kw_p = kw$p.value,
                                       directional_mean = mean(x[directional]), uncertain_mean = mean(x[!directional]),
                                       mw_p = wt$p.value)
                          }))
  }))
  f_likert_tests <- file.path(path_csv, "synthesis-debrief-likert-tests.csv")
  fwrite(likert_tests, f_likert_tests)
  msg("Saved: ", f_likert_tests)
  
  # Pairwise class contrasts: mean difference with a percentile-bootstrap
  # 95% CI (groups resampled separately). Source for the intervals quoted in the text.
  set.seed(as.integer(cfg$run$seed))
  n_boot <- as.integer(cfg$synthesis$debrief$n_boot)
  
  contrast_specs <- list(
    c("imb_class", "with_imb",    "uncertain"),
    c("imb_class", "against_imb", "uncertain"),
    c("imb_class", "with_imb",    "against_imb"),
    c("run_class", "against_run", "uncertain")
  )
  likert_contrasts <- rbindlist(lapply(c("belief_independence", "reliance_on_sequence"), function(item) {
    rbindlist(lapply(contrast_specs, function(s) {
      x1 <- as.numeric(full_q[[item]][full_q[[s[1]]] == s[2]])
      x2 <- as.numeric(full_q[[item]][full_q[[s[1]]] == s[3]])
      d  <- replicate(n_boot, mean(sample(x1, replace = TRUE)) - mean(sample(x2, replace = TRUE)))
      data.table(item = item, class_dim = s[1], group1 = s[2], group2 = s[3],
                 n1 = length(x1), n2 = length(x2), mean1 = mean(x1), mean2 = mean(x2),
                 mean_diff = mean(x1) - mean(x2),
                 ci_low = unname(quantile(d, 0.025)), ci_high = unname(quantile(d, 0.975)),
                 mw_p = suppressWarnings(wilcox.test(x1, x2)$p.value))
    }))
  }))
  f_likert_contrasts <- file.path(path_csv, "synthesis-debrief-likert-contrasts.csv")
  fwrite(likert_contrasts, f_likert_contrasts)
  msg("Saved: ", f_likert_contrasts)
  
  # Likert items by the pass-(2) 5-way self-report category (# yes-and-codable only -- matches merged_all's restriction)
  cat_likert <- merge(merged_all[, .(pid, category)],
                      full_q[, c("pid", likert_items), with = FALSE], by = "pid")
  long_cat <- melt(cat_likert, id.vars = c("pid", "category"), variable.name = "item", value.name = "value")
  likert_by_cat <- long_cat[, .(mean = mean(value), median = median(as.numeric(value)), n = .N),
                            by = .(item, category)]
  setorder(likert_by_cat, item)
  f_likert_cat <- file.path(path_csv, "synthesis-debrief-likert-by-category.csv")
  fwrite(likert_by_cat, f_likert_cat)
  msg("Saved: ", f_likert_cat)
  
  invisible(list(used_strategy_crosstab_imb = used_imb_tbl, used_strategy_crosstab_run = used_run_tbl,
                 used_strategy_tests = used_test_tbl, merged = merged_all,
                 crosstab_imb = crosstabs_imb, crosstab_run = crosstabs_run,
                 crosstab_imb_all = crosstab_imb_all, crosstab_run_all = crosstab_run_all,
                 fisher = fisher_tbl, full_questionnaire = full_q,
                 likert_by_class = likert_by_class, likert_tests = likert_tests,
                 likert_contrasts = likert_contrasts,
                 likert_by_category = likert_by_cat))
}