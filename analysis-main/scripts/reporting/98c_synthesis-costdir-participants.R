# ============================================================
# 98c_synthesis-costdir-participants.R
#
# PURPOSE
#   Synthesis, participant level: are the participants who lose welfare the
#   same as those whose side choice responds to the sequence? Participants
#   are a sample, so each posterior draw is combined with Bayesian-bootstrap
#   weights over participants (Rubin 1981); intervals reflect estimation and
#   sampling uncertainty. Intervals from the posterior alone are also
#   reported.
#   Cue responsiveness comes from the personal side-choice slopes of the
#   joint model:
#     run responsiveness  = -(beta_side[2] + u[5, i]), higher = bets more
#                           strongly against the terminal run
#     imbalance response  = |beta_side[3] + u[6, i]|, higher = side choice more
#                           driven by the imbalance, in either direction
#   B1  Spearman correlation of each responsiveness with welfare loss, betting
#       rate and absolute stake deviation; also for participants with at
#       least cfg$synthesis$costdir$min_bets bets.
#   B2  Weighted OLS of log welfare loss on both responsiveness measures,
#       then adding the number of bets, side preference |mu_h_i - hbar| and
#       the betting rate. Written as a LaTeX table.
#   B3  Outcomes by slope class. Classes follow the 95% interval of each
#       personal slope over the K draws (credibly > 0 = with, credibly < 0 =
#       against, else uncertain). Per draw and weight vector, the weighted
#       mean and median of each outcome in two groups and their difference
#       are computed (welfare loss, betting rate and absolute stake deviation
#       in percentage points; and the raw descriptives within-person SD of
#       positive stakes, share of focal stakes 50 or 100, mean positive
#       stake, median decision time in seconds, and the CRRA coefficient
#       r_mean of the MPL model). Descriptive, no model.
#
# SAMPLE
#   Treatments with a confirmatory sample (cfg$design$a_flags$betting_normative):
#   both the confirmatory and the full fits of RQ1-RQ3; other treatments: the
#   full fits. The joint model and RQ4 are always from the full sample. The
#   LaTeX table uses the confirmatory (or, if none, full) fits.
#
# INPUT
#   models selected under rq1, rq2, rq3, rq4, synthesis, mpl
#     (mu_b_i, mu_a_i, mu_c_i, mu_h_i and hbar; beta_side and u; scored-<tr>.csv)
#   data/derived/master-sequences.csv   stakes and decision times per trial
#
# OUTPUT
#   data/output/csv/synthesis-<tr>-costdir-pid-spearman.csv
#   data/output/csv/synthesis-<tr>-costdir-pid-ols.csv
#   data/output/csv/synthesis-<tr>-costdir-pid-groups.csv
#   data/output/tables/synthesis-<tr>-costdir-pid.tex
#
# NOTES
#   - Participant effects from separately estimated models are shrunk
#     independently, which biases cross-model correlations toward zero, so
#     the estimates are conservative in magnitude
#   - The K draws (cfg$synthesis$costdir) are subsampled independently from
#     each fit and paired by row index
#   - Seeded by cfg$run$seed
#   - Output files are overwritten at every run
# ============================================================

synthesis_costdir_participants <- function(cfg) {
  
  tr_vec <- unique(as.character(cfg$run$treatment))
  set.seed(as.integer(cfg$run$seed))
  tl <- synthesis_draw_tools(cfg)
  K  <- tl$K; min_bets <- tl$min_bets
  wspear <- tl$wspear; wmean <- tl$wmean; wmed <- tl$wmed
  
  master <- fread(file.path(path_der, "master-sequences.csv"), encoding = "UTF-8")
  master[, `:=`(pid = as.character(pid), treat = as.character(treat),
                stake = fifelse(is.na(as.numeric(stake)), 0, as.numeric(stake)))]
  
  out <- list()
  
  for (tr in tr_vec) {
    
    conf_ok <- isTRUE(cfg$design$a_flags$betting_normative[[tr]])
    tg_main <- if (conf_ok) "confirmatory" else "full"
    tgs     <- if (conf_ok) c("confirmatory", "full") else "full"
    
    R4  <- tl$load_rq4(tr)
    SYN <- tl$load_syn(tr)
    L   <- list()
    for (tg in tgs)
      L[[tg]] <- list(rq1 = tl$load_rq("rq1", tr, tg, "mu_b", "mu_b_i"),
                      rq2 = tl$load_rq("rq2", tr, tg, "mu_a", "mu_a_i"),
                      rq3 = tl$load_rq("rq3", tr, tg, "mu_c", "mu_c_i"))
    if (is.null(R4) || is.null(SYN) ||
        any(vapply(L[[tg_main]], is.null, logical(1L)))) {
      warning("Synthesis cost vs direction (participants): fits missing for ", tr, "; skipped.",
              call. = FALSE)
      next
    }
    
    nb_all <- master[treat == tr, .(n_bets = sum(stake > 0)), by = pid]
    
    b_sp <- list(); b_ols <- list(); b_grp <- list()
    for (tg in tgs) {
      S_ <- L[[tg]]
      if (any(vapply(S_, is.null, logical(1L)))) next
      P <- Reduce(intersect, list(colnames(SYN$run), colnames(R4$pid),
                                  colnames(S_$rq1$pid), colnames(S_$rq2$pid), colnames(S_$rq3$pid)))
      P <- sort(P); n <- length(P)
      run_r <- -SYN$run[, P]
      imb_r <- abs(SYN$imb[, P])
      lossP <- S_$rq3$pid[, P]
      betP  <- S_$rq1$pid[, P]
      stkP  <- abs(S_$rq2$pid[, P])
      sideP <- abs(R4$pid[, P] - R4$hbar)
      nbets <- nb_all[match(P, pid), n_bets]
      keepB <- nbets >= min_bets
      
      W <- t(replicate(K, tl$dirichlet(n)))                  # K x n bootstrap weights
      U <- matrix(1 / n, K, n)                               # equal weights (posterior only)
      
      # ---- B1: Spearman correlations ----
      for (unc in c("posterior+bootstrap", "posterior")) {
        Wm <- if (unc == "posterior") U else W
        for (x in c("run_resp", "imb_resp")) {
          Xm <- if (x == "run_resp") run_r else imb_r
          for (y in c("loss", "bet", "abs_stake")) {
            Ym <- switch(y, loss = lossP, bet = betP, abs_stake = stkP)
            r  <- sapply(1:K, function(k) wspear(Xm[k, ], Ym[k, ], Wm[k, ]))
            b_sp[[length(b_sp) + 1L]] <- tl$sumrow(r, treatment = tr, sample = tg, subset = "all",
                                                   uncertainty = unc, x = x, y = y, n = n)
          }
          r <- sapply(1:K, function(k) {
            w <- Wm[k, keepB]; wspear(Xm[k, keepB], lossP[k, keepB], w / sum(w)) })
          b_sp[[length(b_sp) + 1L]] <- tl$sumrow(r, treatment = tr, sample = tg,
                                                 subset = paste0("n_bets >= ", min_bets),
                                                 uncertainty = unc, x = x, y = "loss", n = sum(keepB))
        }
      }
      
      # ---- B2: weighted OLS of log welfare loss ----
      specs <- list(base          = character(0),
                    n_bets        = "n_bets",
                    n_bets_side   = c("n_bets", "side_pref"),
                    n_bets_side_b = c("n_bets", "side_pref", "bet"))
      for (sp in names(specs)) {
        cf <- sapply(1:K, function(k) {
          X <- cbind(1, run_resp = run_r[k, ], imb_resp = imb_r[k, ],
                     n_bets = nbets, side_pref = sideP[k, ], bet = betP[k, ])
          X <- X[, c(1, 2, 3, match(specs[[sp]], colnames(X))), drop = FALSE]
          colnames(X)[1] <- "(Intercept)"
          lm.wfit(X, log(pmax(lossP[k, ], 1e-6)), W[k, ])$coefficients
        })
        b_ols[[length(b_ols) + 1L]] <- rbindlist(lapply(rownames(cf), function(t)
          tl$sumrow(cf[t, ], treatment = tr, sample = tg, spec = sp, term = t, n = n)))
      }
      
      # ---- B3: outcomes by slope class ----
      # For the imbalance slope, "against" = counters the majority
      slope_class <- function(M) {
        lo <- apply(M, 2, quantile, probs = .025)
        hi <- apply(M, 2, quantile, probs = .975)
        ifelse(lo > 0, "with", ifelse(hi < 0, "against", "uncertain"))
      }
      run_cls <- slope_class(SYN$run[, P])
      imb_cls <- slope_class(SYN$imb[, P])
      contrasts <- list(
        list(lab = "against imbalance vs. uncertain on imbalance",
             A = imb_cls == "against", B = imb_cls == "uncertain"),
        list(lab = "against imbalance vs. with imbalance",
             A = imb_cls == "against", B = imb_cls == "with"),
        list(lab = "against imbalance vs. all others",
             A = imb_cls == "against", B = imb_cls != "against"),
        list(lab = "against imbalance and against run vs. all others",
             A = imb_cls == "against" & run_cls == "against",
             B = !(imb_cls == "against" & run_cls == "against")),
        list(lab = "against run vs. uncertain on run",
             A = run_cls == "against", B = run_cls == "uncertain"),
        list(lab = "uncertain on both vs. all others",
             A = imb_cls == "uncertain" & run_cls == "uncertain",
             B = !(imb_cls == "uncertain" & run_cls == "uncertain")))
      # model-based outcomes (percentage points)
      b_out <- list(loss = 100 * lossP, bet = 100 * betP, abs_stake = 100 * stkP)
      # raw participant descriptives (identical across draws); NA when fewer
      # than two bets for the stake SD
      raw <- master[treat == tr & pid %in% P, {
        st <- stake[stake > 0]
        .(stake_sd   = if (length(st) >= 2) sd(st) else NA_real_,
          focal      = if (length(st) >= 1) 100 * mean(st %in% c(50, 100)) else NA_real_,
          stake_mean = if (length(st) >= 1) mean(st) else NA_real_,
          rt_s       = median(screen_ms, na.rm = TRUE) / 1000)
      }, by = pid][match(P, pid)]
      for (v in c("stake_sd", "focal", "stake_mean", "rt_s"))
        b_out[[v]] <- matrix(raw[[v]], K, n, byrow = TRUE)
      # CRRA coefficient of the MPL model; NA without an MPL estimate
      f_mpl <- tryCatch(model_path("mpl", paste0("scored-", tr, ".csv")),
                        error = function(e) NA_character_)
      if (!is.na(f_mpl)) {
        mpl_r <- fread(f_mpl, encoding = "UTF-8")[, .(pid = as.character(pid), r_mean)]
        b_out[["r_mean"]] <- matrix(mpl_r$r_mean[match(P, mpl_r$pid)], K, n, byrow = TRUE)
      }
      for (unc in c("posterior+bootstrap", "posterior")) {
        Wm <- if (unc == "posterior") U else W
        for (oc in names(b_out)) {
          Ym <- b_out[[oc]]
          for (ct in contrasts) {
            if (sum(ct$A) < 3L || sum(ct$B) < 3L) next
            res <- sapply(1:K, function(k) {
              ia <- ct$A & !is.na(Ym[k, ]); ib <- ct$B & !is.na(Ym[k, ])
              xa <- Ym[k, ia]; wa <- Wm[k, ia]
              xb <- Ym[k, ib]; wb <- Wm[k, ib]
              c(mean_a = wmean(xa, wa), mean_b = wmean(xb, wb),
                mean_diff = wmean(xa, wa) - wmean(xb, wb),
                median_a = wmed(xa, wa), median_b = wmed(xb, wb),
                median_diff = wmed(xa, wa) - wmed(xb, wb))
            })
            b_grp[[length(b_grp) + 1L]] <- rbindlist(lapply(rownames(res), function(q)
              tl$sumrow(res[q, ], treatment = tr, sample = tg, uncertainty = unc, outcome = oc,
                        contrast = ct$lab, quantity = q, n_a = sum(ct$A & !is.na(Ym[1, ])),
                        n_b = sum(ct$B & !is.na(Ym[1, ])), n = n)))
          }
        }
      }
    }
    b_sp <- rbindlist(b_sp); b_ols <- rbindlist(b_ols); b_grp <- rbindlist(b_grp)
    b_sp[, n_draws := K]
    b_ols[, `:=`(n_draws = K, uncertainty = "posterior+bootstrap")]
    
    f1 <- file.path(path_csv, paste0("synthesis-", tr, "-costdir-pid-spearman.csv"))
    f2 <- file.path(path_csv, paste0("synthesis-", tr, "-costdir-pid-ols.csv"))
    fwrite(b_sp,  f1);  msg("Saved: ", f1)
    fwrite(b_ols, f2);  msg("Saved: ", f2)
    if (nrow(b_grp)) {
      b_grp[, n_draws := K]
      f3 <- file.path(path_csv, paste0("synthesis-", tr, "-costdir-pid-groups.csv"))
      fwrite(b_grp, f3);  msg("Saved: ", f3)
    }
    
    # ---- LaTeX table (main sample) ----
    O <- b_ols[sample == tg_main]
    if (nrow(O)) {
      sp_names <- c("base", "n_bets", "n_bets_side", "n_bets_side_b")
      tlab <- c(run_resp = "Run responsiveness", imb_resp = "Imbalance responsiveness",
                n_bets = "Number of bets", side_pref = "Side preference $|\\hat\\mu^h_i-\\bar h|$",
                bet = "Betting rate $\\hat\\mu^b_i$")
      fmt3 <- function(x) sprintf("%.3f", ifelse(abs(x) < 5e-4, 0, x))
      body <- unlist(lapply(names(tlab), function(t) {
        cs <- lapply(sp_names, function(s) {
          r <- O[spec == s & term == t]
          if (nrow(r)) tl$cell(r, fmt3) else c("", "")
        })
        c(paste(tlab[[t]], "&", paste(sapply(cs, `[`, 1), collapse = " & "), "\\\\"),
          paste(" &", paste(sapply(cs, `[`, 2), collapse = " & "), "\\\\[3pt]"))
      }))
      Sp <- b_sp[sample == tg_main & subset == "all" & y == "loss" & uncertainty == "posterior+bootstrap"]
      rho_line <- function(xx, lab) {
        r <- Sp[Sp$x == xx]
        paste0(lab, " & \\multicolumn{4}{c}{$", tl$fmt(r$median), "$ \\; $[", tl$fmt(r$q025), ", ",
               tl$fmt(r$q975), "]$} \\\\")
      }
      tl$write_tex(
        file.path(path_tab, paste0("synthesis-", tr, "-costdir-pid.tex")),
        caption = "Participant-level cue responsiveness and welfare loss",
        label   = paste0("Tab:synthesis-costdir-pid", if (tr != tr_vec[1L]) paste0("-", tr)),
        colspec = "l*{4}{c}",
        header  = c(" & \\multicolumn{4}{c}{Dependent variable: $\\log \\hat\\mu^c_i$} \\\\",
                    "\\cmidrule(lr){2-5}", " & (1) & (2) & (3) & (4) \\\\"),
        body    = body,
        footer  = c(rho_line("run_resp", "Spearman $\\rho$(run resp., loss)"),
                    rho_line("imb_resp", "Spearman $\\rho$(imb. resp., loss)"),
                    paste("Participants &", paste(rep(O$n[1], 4), collapse = " & "), "\\\\")),
        note    = paste0("Run responsiveness is minus the participant's personal side-choice slope on the ",
                         "signed terminal-run cue; imbalance responsiveness is the absolute personal slope on ",
                         "signed imbalance (both from the joint synthesis model). Welfare loss $\\hat\\mu^c_i$ ",
                         "(share of the endowment), betting rate, and side preference come from the RQ3, RQ1, and ",
                         "RQ4 models. Each statistic is computed within each of ", K, " posterior draws, combined ",
                         "with Bayesian-bootstrap weights over participants; cells report the posterior median and ",
                         "95\\% credible interval, and an asterisk denotes an interval excluding zero. Because the ",
                         "inputs come from separately estimated models, correlations are attenuated toward zero. ",
                         if (conf_ok) "Confirmatory" else "Full", " sample."))
    }
    
    out[[tr]] <- list(pid_spearman = b_sp, pid_ols = b_ols, pid_groups = b_grp)
  }
  
  invisible(out)
}
