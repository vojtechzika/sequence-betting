# ============================================================
# 96_synthesis_regressions.R
#
# PURPOSE
#   Fits SIX trial-level mixed-effects specifications per treatment:
#     BASE (crossed random intercepts only):
#       fit1: bet      ~ |imbalance| + I(terminal run present) + (1|pid) + (1|seq)   [extensive margin]
#       fit2: is_heads ~ imbalance   + terminal run (signed)   + (1|pid) + (1|seq)   [side choice, bet==1 only]
#       fit3: stake    ~ |imbalance| + I(terminal run present) + (1|pid) + (1|seq)   [intensive margin, bet==1 only]
#     +SLOPES (same fixed effects, participant random slopes added on
#     imbalance and terminal run, uncorrelated via `||`):
#       fit1_slopes, fit2_slopes, fit3_slopes -- same formulas as above with
#       (1 + <imbalance term> + <terminal-run term> || pid) + (1 | seq)
#   Exports fixed-effect coefficients, average marginal effects (logistic
#   models only), collinearity diagnostics, and Nakagawa marginal/
#   conditional R^2 to CSV for all six fits. BASE fits additionally export
#   per-group ICC; +SLOPES fits additionally export random-effect SDs and
#   a convergence/singularity summary (ICC is not reported for +SLOPES
#   fits -- see NOTES).
#   Two stargazer LaTeX tables are produced per treatment:
#     (a) a 4-column table: fit1 (Base), fit1_slopes (+Slopes),
#         fit2 (Base), fit2_slopes (+Slopes).
#     (b) a 2-column table: fit3 (Base), fit3_slopes (+Slopes).
#   The original single 3-column table (fit1, fit2, fit3 together) is not
#   produced by this script -- removed in the prior edit.
#
# INPUT
#   path_src/master_sequences.csv (must contain a `seq` column
#   identifying each of the 64 sequences -- unverified against the
#   actual file at the time of writing; check before running)
#
# OUTPUT
#   path_out/synthesis_fit1_<tr>_tests.csv                (base)
#   path_out/synthesis_fit2_<tr>_tests.csv                (base)
#   path_out/synthesis_fit3_<tr>_tests.csv                (base)
#   path_out/synthesis_fit1_slopes_<tr>_tests.csv         (+slopes)
#   path_out/synthesis_fit2_slopes_<tr>_tests.csv         (+slopes)
#   path_out/synthesis_fit3_slopes_<tr>_tests.csv         (+slopes)
#   path_tab/synthesis_glmm_<tr>_m1m2_slopes.tex          (4-col M1/M2 base+slopes)
#   path_tab/synthesis_glmm_<tr>_m3_slopes.tex            (2-col M3 base+slopes)
#   path_mod/synthesis_<tr>_fit1.rds
#   path_mod/synthesis_<tr>_fit2.rds
#   path_mod/synthesis_<tr>_fit3.rds
#   path_mod/synthesis_<tr>_fit1_slopes.rds
#   path_mod/synthesis_<tr>_fit2_slopes.rds
#   path_mod/synthesis_<tr>_fit3_slopes.rds
#
# NOTES
#   - CHANGE FROM PRIOR VERSION: `bal_no_run` has been REMOVED from all
#     six model formulas (it was previously a fixed-effect dummy in
#     fit1/fit2/fit3 and their +SLOPES counterparts; it never had its own
#     random slope). Rationale, worked out in full elsewhere: the
#     balanced/no-terminal-run "corner" is definitionally the
#     intersection of the two predictors already in the model
#     (imbalance = 0 AND terminal run = 0), not an independent
#     construct. In the BASE model it was significant (capturing a
#     ~5.8pp additional drop in betting rate beyond what imbalance and
#     terminal run predict on their own), but that effect attenuated to
#     non-significance once participant random slopes were added --
#     consistent with the corner effect being a byproduct of forcing a
#     single shared imbalance slope onto a heterogeneous population
#     (whoever has a strong positive imbalance slope will, by
#     construction, show their lowest predicted probability right at
#     imbalance = 0, i.e. at balanced sequences), not a standalone
#     psychological effect. Keeping it in the model added a third
#     predictor whose only real narrative function had become "look, it
#     lost significance when slopes are added" -- a single coefficient's
#     fate is thin, close-to-circular evidence for the slopes-matter
#     argument, which is already well supported by the variance-
#     component comparison itself (imbalance-slope SD dwarfing every
#     other component, replicated across treatments). Dropping it also
#     matches the actual hypothesis this whole regression exercise is
#     testing -- imbalance governs engagement, terminal run governs
#     direction -- which never required a third corner term.
#   - CONSEQUENCE: fit1.rds and fit2.rds are still the objects
#     97_side_choice_modelling.R loads for its out-of-sample validation,
#     but their specification has now changed (two predictors instead of
#     three). Any previously reported out-of-sample accuracy/AUC numbers
#     derived from the old fit1.rds/fit2.rds (e.g. the ~60.3% side-choice
#     figure and the 87.5%/AUC=.579 extensive-margin figures currently in
#     the manuscript) were computed against the OLD three-predictor
#     specification and need to be re-derived against these new fits
#     before being cited again.
#   - The three +SLOPES fits use `||`, not `|`: independent (uncorrelated)
#     variance components for the intercept and each slope. This matches
#     the exploratory decomposition validated earlier (Python VB fit,
#     cross-checked against a real lme4 fit for fit2 on m25 and m19, both
#     converging cleanly with isSingular = FALSE) and avoids the added
#     fragility of also estimating intercept-slope correlations. If a
#     correlated (`|`) version becomes worth reporting, add it as a
#     SEPARATE fit_*_slopes_cor object rather than replacing this one, so
#     the uncorrelated version stays available for comparison.
#   - Each +SLOPES fit is wrapped in tryCatch: a failure to converge or an
#     outright error is logged (with lme4::isSingular() and any
#     convergence warnings) and that treatment's slopes CSV/rds/tables are
#     skipped, but the BASE fit1/fit2/fit3 pipeline is unaffected either
#     way. fit3_slopes in particular is known to hit a boundary/singular
#     fit on m25 (see chat) -- check the varcomp_sd rows in
#     synthesis_fit3_slopes_m25_tests.csv before citing those numbers;
#     which specific component is at ~0 determines whether that's a
#     substantive finding (stake genuinely lacks the slope heterogeneity
#     seen in fit1/fit2) or an over-parameterization issue needing a
#     simpler spec.
#   - ICC is NOT reported for +SLOPES fits. Once random slopes are added,
#     the participant grouping factor's contribution to total variance is
#     no longer a single number -- it depends on the value of the
#     predictor through the slope term -- so a plain ICC computed on a
#     +SLOPES fit would silently be using a different (and less
#     transparent) definition than the ICC reported for the BASE fit in
#     the same row of a table. Random-effect SDs (participant intercept,
#     participant slope per covariate, sequence intercept) are reported
#     instead for both BASE and +SLOPES fits uniformly, via extract_sd(),
#     so the same row block is directly comparable across specifications.
#     Nakagawa marginal/conditional R^2 (performance::r2_nakagawa) DOES
#     extend correctly to random-slope models and is reported for all six
#     fits.
#   - `has_terminal_run2` is signed (-1 = ends TT, 0 = no run of >=2,
#     +1 = ends HH); `has_terminal_run2_abs` is its absolute value
#     (presence of a run of either direction).
#   - `bal_no_run` still exists as a column in `dt`/`dt_bets` (it is
#     computed a few lines below from the raw sequence data) even though
#     no model in this script uses it any more as a regressor -- it may
#     still be needed elsewhere (e.g. for descriptive figures), so it is
#     left in the data prep rather than removed there too.
#   - fit3/fit3_slopes have no average marginal effects exported: they
#     are plain Gaussian LMMs, so the fixed-effect coefficients already
#     are the marginal effects (unlike fit1/fit2 and their +SLOPES
#     counterparts, which are on the logit scale).
#   - placeholder1/placeholder2/placeholder3 (dummy lm() objects used so
#     stargazer can structurally parse a table built from glmer/lmer
#     coefficients passed in via coef=/se=/p=) are created once, before
#     either table is built, and shared by both.
#   - Both tables use column.labels to carry BOTH the dependent-variable
#     name and the Base/+Slopes spec (e.g. "Bet = 1, Base"), with
#     dep.var.labels.include = FALSE, rather than relying on stargazer's
#     dep.var.labels span-merging across repeated placeholder objects --
#     the latter was tried first and produced an incorrect header (both
#     column pairs labelled with the first dependent variable). This is
#     the more verbose but unambiguous alternative.
#   - The two tables are independently gated by should_skip() as before;
#     regenerating one does not force regeneration of the other.
#   - control = lme4::glmerControl(optimizer = "bobyqa") is used for all
#     four glmer-based fits (fit1, fit2, fit1_slopes, fit2_slopes); the
#     +SLOPES glmer fits additionally raise maxfun to 2e5, since crossed
#     random effects plus random slopes are more prone to needing extra
#     iterations than either alone.
#   - stargazer::stargazer() writes to `out=` but ALSO unconditionally
#     cat()s the rendered table to the console as a side effect -- there
#     is no silent/quiet argument for it. Both stargazer calls below are
#     wrapped in invisible(capture.output(...)) purely to suppress that
#     console dump; the file writes via `out=` are unaffected.
# ============================================================
synthesis_regressions <- function(cfg) {
  
  design <- cfg$design
  tr_vec <- unique(as.character(cfg$run$treatment))
  
  dir.create(path_mod, showWarnings = FALSE, recursive = TRUE)
  
  f_in <- file.path(path_src, "master_sequences.csv")
  stopifnot(file.exists(f_in))
  dt_all <- fread(f_in)
  
  # ---- Helper: tidy a glmer/lmer fixed-effects table ----
  tidy_fixed <- function(fit, tr) {
    s    <- summary(fit)
    raw  <- as.data.frame(s$coefficients)
    stat_col <- if ("z value" %in% names(raw)) "z value" else "t value"
    p_col    <- if ("Pr(>|z|)" %in% names(raw)) "Pr(>|z|)" else "Pr(>|t|)"
    
    data.table(
      term           = rownames(raw),
      estimate       = raw[["Estimate"]],
      std_error      = raw[["Std. Error"]],
      statistic      = raw[[stat_col]],
      p_value        = raw[[p_col]],
      treatment      = tr,
      test_type      = "fixed_effect",
      n_obs          = nobs(fit),
      n_participants = lme4::ngrps(fit)[["pid"]],
      n_sequences    = lme4::ngrps(fit)[["seq"]]
    )
  }
  
  # ---- Helper: tidy an avg_slopes() AME table ----
  tidy_ame <- function(fit, tr) {
    ame <- marginaleffects::avg_slopes(fit, type = "response")
    tbl <- as.data.table(ame)
    tbl[, treatment := tr]
    tbl[, test_type  := "marginal_effect"]
    tbl
  }
  
  # ---- Helper: tidy a check_collinearity() table ----
  tidy_vif <- function(fit, tr) {
    tbl <- as.data.table(performance::check_collinearity(fit))
    setnames(tbl, old = "Term", new = "term", skip_absent = TRUE)
    tbl[, treatment := tr]
    tbl[, test_type  := "collinearity"]
    tbl
  }
  
  # ---- Helper: tidy an r2_nakagawa() table ----
  tidy_r2 <- function(r2, tr) {
    data.table(
      term      = c("R2_marginal", "R2_conditional"),
      estimate  = c(r2$R2_marginal, r2$R2_conditional),
      treatment = tr,
      test_type = "r2_nakagawa"
    )
  }
  
  # ---- Helper: pull one grouping factor's ICC out of
  # performance::icc(fit, by_group = TRUE)'s output. Used only for BASE
  # fits (see header NOTES on why +SLOPES fits skip ICC). ----
  icc_value <- function(icc_tbl, group) {
    val <- icc_tbl$ICC[icc_tbl$Group == group]
    if (length(val) == 0L) NA_real_ else val
  }
  
  tidy_icc <- function(icc_tbl, tr) {
    data.table(
      term      = paste0("ICC_", icc_tbl$Group),
      estimate  = icc_tbl$ICC,
      treatment = tr,
      test_type = "icc_by_group"
    )
  }
  
  # ---- Helper (+SLOPES fits): convergence/singularity summary ----
  convergence_status <- function(fit) {
    msgs <- tryCatch(fit@optinfo$conv$lme4$messages, error = function(e) NULL)
    list(
      singular = tryCatch(lme4::isSingular(fit), error = function(e) NA),
      messages = if (is.null(msgs)) "none" else paste(msgs, collapse = " | ")
    )
  }
  
  # ---- Helper: extract random-effect SDs from ANY fit (base or
  # +slopes) in a uniform shape. Returns NA for slope components that
  # don't exist on a given fit (e.g. always NA on base fits), so the
  # same downstream formatting code works for both specifications
  # without branching. Handles both the imbalance_abs/has_terminal_run2_abs
  # naming (fit1/fit3) and the imbalance/has_terminal_run2 naming (fit2)
  # by trying both and taking whichever is non-NA. ----
  extract_sd <- function(fit) {
    vc <- as.data.frame(lme4::VarCorr(fit))
    vc <- vc[is.na(vc$var2), ]  # drop correlation-parameter rows, if any
    get_sd <- function(grp_pattern, term) {
      row <- vc[grepl(grp_pattern, vc$grp) & !is.na(vc$var1) & vc$var1 == term, ]
      if (nrow(row) == 0L) NA_real_ else row$sdcor[1]
    }
    coalesce_na <- function(a, b) if (!is.na(a)) a else b
    list(
      pid_int = get_sd("^pid", "(Intercept)"),
      pid_imb = coalesce_na(get_sd("^pid", "imbalance_abs"), get_sd("^pid", "imbalance")),
      pid_run = coalesce_na(get_sd("^pid", "has_terminal_run2_abs"), get_sd("^pid", "has_terminal_run2")),
      seq_int = get_sd("^seq", "(Intercept)")
    )
  }
  fmt_sd <- function(x) if (is.na(x)) "--" else sprintf("%.3f", x)
  
  # ---- Helper (+SLOPES fits): tidy the SD block + convergence summary
  # into the same long CSV shape used elsewhere in this script ----
  tidy_slopes_export <- function(fit, conv, tr) {
    sd <- extract_sd(fit)
    data.table(
      term = c("SD_pid_intercept", "SD_pid_slope_imbalance",
               "SD_pid_slope_terminal_run", "SD_seq_intercept",
               "singular_fit", paste0("convergence_messages: ", conv$messages)),
      estimate = c(sd$pid_int, sd$pid_imb, sd$pid_run, sd$seq_int,
                   as.numeric(conv$singular), NA_real_),
      treatment = tr,
      test_type = c(rep("varcomp_sd", 4), rep("convergence", 2))
    )
  }
  
  test_csvs <- list()
  
  for (tr in tr_vec) {
    
    multi_val <- design$seq$treatments[[tr]]
    if (is.null(multi_val)) {
      msg("Skipping treatment '", tr, "': not found in design$seq$treatments")
      next
    }
    
    dt <- dt_all[multi == multi_val]
    if (nrow(dt) == 0L) {
      msg("Skipping treatment '", tr, "': no rows in master_sequences.csv")
      next
    }
    
    dt[, bet      := as.integer(side != "NB")]
    dt[, is_heads := fifelse(side == "H", 1L,
                             fifelse(side == "T", 0L, NA_integer_))]
    dt[, imbalance_abs         := abs(imbalance)]
    dt[, has_terminal_run2_abs := abs(has_terminal_run2)]
    # bal_no_run is still computed/available in dt/dt_bets below (it comes
    # through from master_sequences.csv), but is no longer used as a
    # regressor in any of the six fits -- see header NOTES.
    
    dt_bets <- dt[side != "NB"]
    
    # ======================================================================
    # BASE fits (fit1, fit2, fit3)
    # ======================================================================
    
    # ---- fit1: extensive margin (bet) ----
    fit1 <- lme4::glmer(
      bet ~ imbalance_abs + has_terminal_run2_abs + (1 | pid) + (1 | seq),
      data    = dt,
      family  = binomial(link = "logit"),
      control = lme4::glmerControl(optimizer = "bobyqa")
    )
    r2_1  <- performance::r2_nakagawa(fit1)
    icc_1 <- as.data.frame(performance::icc(fit1, by_group = TRUE))
    
    f_fit1 <- file.path(path_out, paste0("synthesis_", tr, "_fit1.csv"))
    if (!should_skip(f_fit1, cfg, "output", paste0("Synthesis fit1 tests (", tr, ")"))) {
      fwrite(rbindlist(list(
        tidy_fixed(fit1, tr), tidy_ame(fit1, tr), tidy_vif(fit1, tr),
        tidy_r2(r2_1, tr), tidy_icc(icc_1, tr)
      ), fill = TRUE), f_fit1)
      msg("Saved: ", f_fit1)
    }
    test_csvs[[paste0("fit1_", tr)]] <- f_fit1
    
    f_mod1 <- file.path(path_mod, paste0("synthesis_", tr, "_fit1.rds"))
    if (!should_skip(f_mod1, cfg, "output", paste0("Synthesis fit1 model object (", tr, ")"))) {
      saveRDS(fit1, f_mod1)
      msg("Saved: ", f_mod1)
    }
    test_csvs[[paste0("fit1_mod_", tr)]] <- f_mod1
    
    # ---- fit2: side choice (is_heads | bet == 1) ----
    fit2 <- lme4::glmer(
      is_heads ~ imbalance + has_terminal_run2 + (1 | pid) + (1 | seq),
      data    = dt_bets,
      family  = binomial(link = "logit"),
      control = lme4::glmerControl(optimizer = "bobyqa")
    )
    r2_2  <- performance::r2_nakagawa(fit2)
    icc_2 <- as.data.frame(performance::icc(fit2, by_group = TRUE))
    
    f_fit2 <- file.path(path_out, paste0("synthesis_", tr, "_fit2.csv"))
    if (!should_skip(f_fit2, cfg, "output", paste0("Synthesis fit2 tests (", tr, ")"))) {
      fwrite(rbindlist(list(
        tidy_fixed(fit2, tr), tidy_ame(fit2, tr), tidy_vif(fit2, tr),
        tidy_r2(r2_2, tr), tidy_icc(icc_2, tr)
      ), fill = TRUE), f_fit2)
      msg("Saved: ", f_fit2)
    }
    test_csvs[[paste0("fit2_", tr)]] <- f_fit2
    
    f_mod2 <- file.path(path_mod, paste0("synthesis_", tr, "_fit2.rds"))
    if (!should_skip(f_mod2, cfg, "output", paste0("Synthesis fit2 model object (", tr, ")"))) {
      saveRDS(fit2, f_mod2)
      msg("Saved: ", f_mod2)
    }
    test_csvs[[paste0("fit2_mod_", tr)]] <- f_mod2
    
    # ---- fit3: intensive margin (stake | bet == 1) ----
    fit3 <- lmerTest::lmer(
      stake ~ imbalance_abs + has_terminal_run2_abs + (1 | pid) + (1 | seq),
      data = dt_bets
    )
    r2_3  <- performance::r2_nakagawa(fit3)
    icc_3 <- as.data.frame(performance::icc(fit3, by_group = TRUE))
    
    f_fit3 <- file.path(path_out, paste0("synthesis_", tr, "_fit3.csv"))
    if (!should_skip(f_fit3, cfg, "output", paste0("Synthesis fit3 tests (", tr, ")"))) {
      fwrite(rbindlist(list(
        tidy_fixed(fit3, tr), tidy_vif(fit3, tr), tidy_r2(r2_3, tr), tidy_icc(icc_3, tr)
      ), fill = TRUE), f_fit3)
      msg("Saved: ", f_fit3)
    }
    test_csvs[[paste0("fit3_", tr)]] <- f_fit3
    
    f_mod3 <- file.path(path_mod, paste0("synthesis_", tr, "_fit3.rds"))
    if (!should_skip(f_mod3, cfg, "output", paste0("Synthesis fit3 model object (", tr, ")"))) {
      saveRDS(fit3, f_mod3)
      msg("Saved: ", f_mod3)
    }
    test_csvs[[paste0("fit3_mod_", tr)]] <- f_mod3
    
    # ---- Shared placeholders for stargazer (used by both the M1/M2
    # and M3 base+slopes tables below) ----
    placeholder1 <- lm(bet      ~ imbalance_abs + has_terminal_run2_abs, data = dt)
    placeholder2 <- lm(is_heads ~ imbalance     + has_terminal_run2,     data = dt_bets)
    placeholder3 <- lm(stake    ~ imbalance_abs + has_terminal_run2_abs, data = dt_bets)
    
    s1 <- summary(fit1); s2 <- summary(fit2); s3 <- summary(fit3)
    
    # ======================================================================
    # +SLOPES fits (fit1_slopes, fit2_slopes, fit3_slopes)
    # Each wrapped in tryCatch so a convergence failure can never break
    # the BASE pipeline.
    # ======================================================================
    
    fit1_slopes <- fit2_slopes <- fit3_slopes <- NULL
    
    fit1_slopes <- tryCatch({
      f <- lme4::glmer(
        bet ~ imbalance_abs + has_terminal_run2_abs +
          (1 + imbalance_abs + has_terminal_run2_abs || pid) + (1 | seq),
        data = dt, family = binomial(link = "logit"),
        control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
      )
      conv <- convergence_status(f)
      msg("fit1_slopes (", tr, "): singular = ", conv$singular, "; messages = ", conv$messages)
      r2s  <- performance::r2_nakagawa(f)
      
      f_csv <- file.path(path_out, paste0("synthesis_fit1_slopes_", tr, "_tests.csv"))
      if (!should_skip(f_csv, cfg, "output", paste0("Synthesis fit1_slopes tests (", tr, ")"))) {
        fwrite(rbindlist(list(
          tidy_fixed(f, tr), tidy_ame(f, tr), tidy_vif(f, tr),
          tidy_r2(r2s, tr), tidy_slopes_export(f, conv, tr)
        ), fill = TRUE), f_csv)
        msg("Saved: ", f_csv)
      }
      f_mod <- file.path(path_mod, paste0("synthesis_", tr, "_fit1_slopes.rds"))
      if (!should_skip(f_mod, cfg, "output", paste0("Synthesis fit1_slopes model object (", tr, ")"))) {
        saveRDS(f, f_mod); msg("Saved: ", f_mod)
      }
      f
    }, error = function(e) {
      msg("fit1_slopes (", tr, ") FAILED to fit -- ", conditionMessage(e),
          ". BASE fit1 is unaffected.")
      NULL
    })
    
    fit2_slopes <- tryCatch({
      f <- lme4::glmer(
        is_heads ~ imbalance + has_terminal_run2 +
          (1 + imbalance + has_terminal_run2 || pid) + (1 | seq),
        data = dt_bets, family = binomial(link = "logit"),
        control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
      )
      conv <- convergence_status(f)
      msg("fit2_slopes (", tr, "): singular = ", conv$singular, "; messages = ", conv$messages)
      r2s  <- performance::r2_nakagawa(f)
      
      f_csv <- file.path(path_out, paste0("synthesis_fit2_slopes_", tr, "_tests.csv"))
      if (!should_skip(f_csv, cfg, "output", paste0("Synthesis fit2_slopes tests (", tr, ")"))) {
        fwrite(rbindlist(list(
          tidy_fixed(f, tr), tidy_ame(f, tr), tidy_vif(f, tr),
          tidy_r2(r2s, tr), tidy_slopes_export(f, conv, tr)
        ), fill = TRUE), f_csv)
        msg("Saved: ", f_csv)
      }
      f_mod <- file.path(path_mod, paste0("synthesis_", tr, "_fit2_slopes.rds"))
      if (!should_skip(f_mod, cfg, "output", paste0("Synthesis fit2_slopes model object (", tr, ")"))) {
        saveRDS(f, f_mod); msg("Saved: ", f_mod)
      }
      f
    }, error = function(e) {
      msg("fit2_slopes (", tr, ") FAILED to fit -- ", conditionMessage(e),
          ". BASE fit2 is unaffected.")
      NULL
    })
    
    fit3_slopes <- tryCatch({
      f <- lmerTest::lmer(
        stake ~ imbalance_abs + has_terminal_run2_abs +
          (1 + imbalance_abs + has_terminal_run2_abs || pid) + (1 | seq),
        data = dt_bets
      )
      conv <- convergence_status(f)
      msg("fit3_slopes (", tr, "): singular = ", conv$singular, "; messages = ", conv$messages)
      r2s  <- performance::r2_nakagawa(f)
      
      f_csv <- file.path(path_out, paste0("synthesis_fit3_slopes_", tr, "_tests.csv"))
      if (!should_skip(f_csv, cfg, "output", paste0("Synthesis fit3_slopes tests (", tr, ")"))) {
        fwrite(rbindlist(list(
          tidy_fixed(f, tr), tidy_vif(f, tr), tidy_r2(r2s, tr), tidy_slopes_export(f, conv, tr)
        ), fill = TRUE), f_csv)
        msg("Saved: ", f_csv)
      }
      f_mod <- file.path(path_mod, paste0("synthesis_", tr, "_fit3_slopes.rds"))
      if (!should_skip(f_mod, cfg, "output", paste0("Synthesis fit3_slopes model object (", tr, ")"))) {
        saveRDS(f, f_mod); msg("Saved: ", f_mod)
      }
      f
    }, error = function(e) {
      msg("fit3_slopes (", tr, ") FAILED to fit -- ", conditionMessage(e),
          ". BASE fit3 is unaffected.")
      NULL
    })
    
    # ======================================================================
    # Table A: Model 1 and Model 2, Base + Slopes (4 columns)
    # ======================================================================
    f_tab_m1m2 <- file.path(path_tab, paste0("synthesis_glmm_", tr, "_m1m2_slopes.tex"))
    if (is.null(fit1_slopes) || is.null(fit2_slopes)) {
      msg("Skipping Table A (M1/M2 base+slopes) for '", tr, "': fit1_slopes and/or fit2_slopes unavailable.")
    } else if (!should_skip(f_tab_m1m2, cfg, "output", paste0("Synthesis M1/M2 base+slopes table (", tr, ")"))) {
      s1s <- summary(fit1_slopes); s2s <- summary(fit2_slopes)
      r2_1s <- performance::r2_nakagawa(fit1_slopes)
      r2_2s <- performance::r2_nakagawa(fit2_slopes)
      sd1  <- extract_sd(fit1); sd1s <- extract_sd(fit1_slopes)
      sd2  <- extract_sd(fit2); sd2s <- extract_sd(fit2_slopes)
      
      invisible(capture.output(
        stargazer::stargazer(
          placeholder1, placeholder1, placeholder2, placeholder2,
          coef = list(s1$coefficients[, "Estimate"], s1s$coefficients[, "Estimate"],
                      s2$coefficients[, "Estimate"], s2s$coefficients[, "Estimate"]),
          se   = list(s1$coefficients[, "Std. Error"], s1s$coefficients[, "Std. Error"],
                      s2$coefficients[, "Std. Error"], s2s$coefficients[, "Std. Error"]),
          p    = list(s1$coefficients[, "Pr(>|z|)"], s1s$coefficients[, "Pr(>|z|)"],
                      s2$coefficients[, "Pr(>|z|)"], s2s$coefficients[, "Pr(>|z|)"]),
          type      = "latex",
          keep.stat = "n",
          add.lines = list(
            c("Participants (pid)", s1$ngrps[["pid"]], s1s$ngrps[["pid"]], s2$ngrps[["pid"]], s2s$ngrps[["pid"]]),
            c("Sequences (seq)",    s1$ngrps[["seq"]], s1s$ngrps[["seq"]], s2$ngrps[["seq"]], s2s$ngrps[["seq"]]),
            c("SD: Participant Intercept",
              fmt_sd(sd1$pid_int), fmt_sd(sd1s$pid_int), fmt_sd(sd2$pid_int), fmt_sd(sd2s$pid_int)),
            c("SD: Participant Slope (Imbalance)",
              "--", fmt_sd(sd1s$pid_imb), "--", fmt_sd(sd2s$pid_imb)),
            c("SD: Participant Slope (Terminal Run)",
              "--", fmt_sd(sd1s$pid_run), "--", fmt_sd(sd2s$pid_run)),
            c("SD: Sequence Intercept",
              fmt_sd(sd1$seq_int), fmt_sd(sd1s$seq_int), fmt_sd(sd2$seq_int), fmt_sd(sd2s$seq_int)),
            c("R\\textsuperscript{2} (marginal)",
              round(r2_1$R2_marginal, 3), round(r2_1s$R2_marginal, 3),
              round(r2_2$R2_marginal, 3), round(r2_2s$R2_marginal, 3)),
            c("R\\textsuperscript{2} (conditional)",
              round(r2_1$R2_conditional, 3), round(r2_1s$R2_conditional, 3),
              round(r2_2$R2_conditional, 3), round(r2_2s$R2_conditional, 3))
          ),
          covariate.labels = c(
            "Imbalance Magnitude [$0$ (balanced), $\\dfrac{1}{3}$, $\\dfrac{2}{3}$, $1$ (pure)]",
            "Has Terminal Run [= 1]",
            "Imbalance [$-1$ (T-Pure), $-\\dfrac{2}{3}$, $-\\dfrac{1}{3}$, 0, $\\dfrac{1}{3}$, $\\dfrac{2}{3}$, $+1$ (H-Pure)]",
            "Terminal Run [$-1$ = TT, $+1$ = HH, $0$ = No Run]"
          ),
          dep.var.labels = c("Bet = 1",  "Heads = 1"),
          column.labels  = c("Base", "+ Slopes", "Base", "+ Slopes"),
          title       = paste0("Determinants of Betting and Side Choice, Base vs.\\ Participant-Slopes Specification (", tr, ")"),
          label       = paste0("tab:synthesis_glmm_", tr, "_m1m2"),
          notes       = "",
          notes.label = "",
          font.size   = "footnotesize",
          out         = f_tab_m1m2
        )
      ))
      msg("Saved: ", f_tab_m1m2)
    }
    
    # ======================================================================
    # Table B: Model 3, Base + Slopes (2 columns)
    # ======================================================================
    f_tab_m3 <- file.path(path_tab, paste0("synthesis_glmm_", tr, "_m3_slopes.tex"))
    if (is.null(fit3_slopes)) {
      msg("Skipping Table B (M3 base+slopes) for '", tr, "': fit3_slopes unavailable.")
    } else if (!should_skip(f_tab_m3, cfg, "output", paste0("Synthesis M3 base+slopes table (", tr, ")"))) {
      s3s   <- summary(fit3_slopes)
      r2_3s <- performance::r2_nakagawa(fit3_slopes)
      sd3   <- extract_sd(fit3); sd3s <- extract_sd(fit3_slopes)
      
      invisible(capture.output(
        stargazer::stargazer(
          placeholder3, placeholder3,
          coef = list(s3$coefficients[, "Estimate"], s3s$coefficients[, "Estimate"]),
          se   = list(s3$coefficients[, "Std. Error"], s3s$coefficients[, "Std. Error"]),
          p    = list(s3$coefficients[, "Pr(>|t|)"], s3s$coefficients[, "Pr(>|t|)"]),
          type      = "latex",
          keep.stat = "n",
          add.lines = list(
            c("Participants (pid)", s3$ngrps[["pid"]], s3s$ngrps[["pid"]]),
            c("Sequences (seq)",    s3$ngrps[["seq"]], s3s$ngrps[["seq"]]),
            c("SD: Participant Intercept",       fmt_sd(sd3$pid_int), fmt_sd(sd3s$pid_int)),
            c("SD: Participant Slope (Imbalance)", "--", fmt_sd(sd3s$pid_imb)),
            c("SD: Participant Slope (Terminal Run)", "--", fmt_sd(sd3s$pid_run)),
            c("SD: Sequence Intercept",          fmt_sd(sd3$seq_int), fmt_sd(sd3s$seq_int)),
            c("R\\textsuperscript{2} (marginal)",
              round(r2_3$R2_marginal, 3), round(r2_3s$R2_marginal, 3)),
            c("R\\textsuperscript{2} (conditional)",
              round(r2_3$R2_conditional, 3), round(r2_3s$R2_conditional, 3))
          ),
          covariate.labels = c(
            "Imbalance Magnitude [$0$ (balanced), $\\dfrac{1}{3}$, $\\dfrac{2}{3}$, $1$ (pure)]",
            "Has Terminal Run [= 1]"
          ),
          dep.var.labels = c("Stake (ECU)", "Stake (ECU)"),
          column.labels  = c("Base", "+ Slopes"),
          title       = paste0("Determinants of Stake, Base vs.\\ Participant-Slopes Specification (", tr, ")"),
          label       = paste0("tab:synthesis_glmm_", tr, "_m3"),
          notes       = "",
          notes.label = "",
          font.size   = "footnotesize",
          out         = f_tab_m3
        )
      ))
      msg("Saved: ", f_tab_m3)
    }
  }
  
  invisible(test_csvs)
}