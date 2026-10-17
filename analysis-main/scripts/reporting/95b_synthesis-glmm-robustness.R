# ============================================================
# 95b_synthesis-glmm-robustness.R
#
# PURPOSE
#   Frequentist robustness check of the joint Bayesian model of the
#   Synthesis section (91_synthesis-stan-joint.R): trial-level mixed-effects
#   models of the same three decisions, with and without participant random
#   slopes, estimated by maximum likelihood (lme4), per treatment.
#     fit1  bet ~ |imbalance| + terminal run present + (1|pid) + (1|seq)
#           extensive margin, all trials, logit
#     fit2  heads ~ imbalance + signed terminal run + (1|pid) + (1|seq)
#           side choice, bets only, logit
#     fit3  stake ~ |imbalance| + terminal run present + (1|pid) + (1|seq)
#           intensive margin, bets only, linear
#     fitk_slopes  the same fixed effects with participant random slopes on
#           the two cues, uncorrelated (the || syntax), plus (1|seq)
#   Reported for every fit: fixed effects, collinearity (VIF), Nakagawa
#   marginal and conditional R2, and the random-effect SDs. The base fits
#   also report the ICC by group (not defined for the slope models, where
#   the participant share of variance depends on the predictor). The logit
#   fits also report average marginal effects. Each slope model is wrapped
#   in tryCatch, so a convergence failure does not affect the base fits; its
#   singularity flag and convergence messages are exported. Two LaTeX tables
#   per treatment compare the base and slope specifications.
#
# SAMPLE
#   synthesis-full of each treatment in cfg$run$treatment. The cues are
#   computed from the sequence string: signed imbalance (#Heads - #Tails)/6
#   and signed terminal run (+1 for a run of at least two Heads, -1 for
#   Tails, 0 for none).
#
# INPUT
#   model "mpl": sample/analysis-trials.csv (through analysis_trials())
#
# OUTPUT
#   data/output/csv/synthesis-<tr>-glmm-<fit>.csv        fit in fit1, fit2, fit3,
#                                                        fit1-slopes, fit2-slopes, fit3-slopes
#   data/output/tables/synthesis-glmm-<tr>-m1m2.tex      bet and side choice, base vs. slopes
#   data/output/tables/synthesis-glmm-<tr>-m3.tex        stake, base vs. slopes
#   data/derived/synthesis-<tr>-glmm-<fit>.rds           the fitted model objects
#
# NOTES
#   - Requires the packages lme4, lmerTest, performance, marginaleffects and
#     stargazer
#   - Optimiser: bobyqa for all logit fits, with 2e5 function evaluations for
#     the slope models
#   - stargazer writes the table to file; its console echo is suppressed
#   - Output files are overwritten at every run
# ============================================================

synthesis_glmm_robustness <- function(cfg) {

  design <- cfg$design
  tr_vec <- unique(as.character(cfg$run$treatment))

  heads <- as.character(design$seq$side_labels$heads)
  tails <- as.character(design$seq$side_labels$tails)

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
  # fits (the ICC is not defined for the slope models, see the header). ----
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

    dt <- analysis_trials(tr, "synthesis-full")
    dt[, `:=`(pid = as.character(pid), seq = as.character(seq))]
    dt[, bet      := as.integer(side %in% c(heads, tails) & is.finite(stake) & stake > 0)]
    dt[, is_heads := fifelse(bet == 1L, as.integer(side == heads), NA_integer_)]

    # Cues from the sequence string: signed imbalance (#Heads - #Tails) / 6 and
    # signed terminal run (+1: ends in >= 2 Heads, -1: >= 2 Tails, 0: none)
    cue <- function(x) {
      ch <- strsplit(x, "", fixed = TRUE)[[1L]]
      n_h <- sum(ch == heads); last <- ch[6L]; k <- 1L
      while (k < 6L && ch[6L - k] == last) k <- k + 1L
      list(imbalance = (2 * n_h - 6) / 6,
           terminal_run = if (k >= 2L) (if (last == heads) 1 else -1) else 0)
    }
    cues <- unique(dt[, .(seq)])
    cues[, c("imbalance", "terminal_run") := rbindlist(lapply(seq, cue))]
    dt <- merge(dt, cues, by = "seq", sort = FALSE)
    dt[, imbalance_abs         := abs(imbalance)]
    dt[, has_terminal_run2     := terminal_run]       # signed: 1 = HH, -1 = TT, 0 = none
    dt[, has_terminal_run2_abs := abs(terminal_run)]  # presence of a run of either direction

    dt_bets <- dt[bet == 1L]

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

    f_fit1 <- file.path(path_csv, paste0("synthesis-", tr, "-glmm-fit1.csv"))
    {
      fwrite(rbindlist(list(
        tidy_fixed(fit1, tr), tidy_ame(fit1, tr), tidy_vif(fit1, tr),
        tidy_r2(r2_1, tr), tidy_icc(icc_1, tr)
      ), fill = TRUE), f_fit1)
      msg("Saved: ", f_fit1)
    }
    test_csvs[[paste0("fit1_", tr)]] <- f_fit1

    f_mod1 <- file.path(path_der, paste0("synthesis-", tr, "-glmm-fit1.rds"))
    {
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

    f_fit2 <- file.path(path_csv, paste0("synthesis-", tr, "-glmm-fit2.csv"))
    {
      fwrite(rbindlist(list(
        tidy_fixed(fit2, tr), tidy_ame(fit2, tr), tidy_vif(fit2, tr),
        tidy_r2(r2_2, tr), tidy_icc(icc_2, tr)
      ), fill = TRUE), f_fit2)
      msg("Saved: ", f_fit2)
    }
    test_csvs[[paste0("fit2_", tr)]] <- f_fit2

    f_mod2 <- file.path(path_der, paste0("synthesis-", tr, "-glmm-fit2.rds"))
    {
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

    f_fit3 <- file.path(path_csv, paste0("synthesis-", tr, "-glmm-fit3.csv"))
    {
      fwrite(rbindlist(list(
        tidy_fixed(fit3, tr), tidy_vif(fit3, tr), tidy_r2(r2_3, tr), tidy_icc(icc_3, tr)
      ), fill = TRUE), f_fit3)
      msg("Saved: ", f_fit3)
    }
    test_csvs[[paste0("fit3_", tr)]] <- f_fit3

    f_mod3 <- file.path(path_der, paste0("synthesis-", tr, "-glmm-fit3.rds"))
    {
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

      f_csv <- file.path(path_csv, paste0("synthesis-", tr, "-glmm-fit1-slopes.csv"))
      {
        fwrite(rbindlist(list(
          tidy_fixed(f, tr), tidy_ame(f, tr), tidy_vif(f, tr),
          tidy_r2(r2s, tr), tidy_slopes_export(f, conv, tr)
        ), fill = TRUE), f_csv)
        msg("Saved: ", f_csv)
      }
      f_mod <- file.path(path_der, paste0("synthesis-", tr, "-glmm-fit1-slopes.rds"))
      {
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

      f_csv <- file.path(path_csv, paste0("synthesis-", tr, "-glmm-fit2-slopes.csv"))
      {
        fwrite(rbindlist(list(
          tidy_fixed(f, tr), tidy_ame(f, tr), tidy_vif(f, tr),
          tidy_r2(r2s, tr), tidy_slopes_export(f, conv, tr)
        ), fill = TRUE), f_csv)
        msg("Saved: ", f_csv)
      }
      f_mod <- file.path(path_der, paste0("synthesis-", tr, "-glmm-fit2-slopes.rds"))
      {
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

      f_csv <- file.path(path_csv, paste0("synthesis-", tr, "-glmm-fit3-slopes.csv"))
      {
        fwrite(rbindlist(list(
          tidy_fixed(f, tr), tidy_vif(f, tr), tidy_r2(r2s, tr), tidy_slopes_export(f, conv, tr)
        ), fill = TRUE), f_csv)
        msg("Saved: ", f_csv)
      }
      f_mod <- file.path(path_der, paste0("synthesis-", tr, "-glmm-fit3-slopes.rds"))
      {
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
    f_tab_m1m2 <- file.path(path_tab, paste0("synthesis-glmm-", tr, "-m1m2.tex"))
    if (is.null(fit1_slopes) || is.null(fit2_slopes)) {
      msg("Skipping Table A (M1/M2 base+slopes) for '", tr, "': fit1_slopes and/or fit2_slopes unavailable.")
    } else {
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
    f_tab_m3 <- file.path(path_tab, paste0("synthesis-glmm-", tr, "-m3.tex"))
    if (is.null(fit3_slopes)) {
      msg("Skipping Table B (M3 base+slopes) for '", tr, "': fit3_slopes unavailable.")
    } else {
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
