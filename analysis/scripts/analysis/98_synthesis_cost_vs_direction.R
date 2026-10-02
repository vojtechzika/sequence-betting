# ============================================================
# 98_synthesis_cost_vs_direction.R
#
# PURPOSE
#   Answers the question that opens the Synthesis section: are the
#   sequences (and participants) that are economically costly the same as
#   those whose side choice is steered by the sequence?
#   Complements z96_synthesis_regressions.R (trial-level GLMMs), which has
#   no welfare-loss equation and no direct cost-vs-direction test.
#
#   Everything is computed WITHIN POSTERIOR DRAWS of the fitted models and
#   summarised as posterior median + 95% credible interval, so posterior
#   uncertainty in every input quantity is propagated. No quantity is
#   treated as known.
#
#   PART A -- SEQUENCE LEVEL (the 64 sequences are the full population
#   of interest, so posterior uncertainty is the only uncertainty)
#     A1  Direct link: side bias |mu_h_s - hbar| vs welfare loss mu_c_s.
#         Per draw: Spearman rho (all 64; excluding the two pure
#         sequences) and OLS slope of loss on bias, alone / + imbalance
#         class / + terminal run.
#     A2  Feature regressions: side bias, betting rate, stake deviation,
#         welfare loss each regressed on a terminal-run indicator and
#         imbalance-CLASS dummies (non-linear on purpose: loss and stake
#         effects are concentrated at 0-6). -> LaTeX table (OA).
#
#   PART B -- PARTICIPANT LEVEL (participants are a sample, so each draw
#   is combined with Bayesian-bootstrap weights over participants
#   (Rubin 1981); intervals then reflect estimation AND sampling
#   uncertainty. Posterior-only intervals are also reported.)
#     Cue responsiveness from the synthesis model's personal slopes:
#       run_resp = -(beta_side[2] + u[5,i])  higher = bets more strongly
#                                            against the terminal run
#       imb_resp = |beta_side[3] + u[6,i]|   higher = side choice more
#                                            imbalance-driven (either way)
#     B1  Spearman rho of each with welfare loss, betting rate, and
#         |stake deviation|; also on participants with >= min_bets bets.
#     B2  OLS log(welfare loss) ~ run_resp + imb_resp, then + n_bets,
#         + side preference |mu_h_i - hbar|, + betting rate.
#         -> LaTeX table (OA).
#     B3  Welfare loss by slope CLASS: do participants who counter the
#         imbalance ("balance-seekers") or bet against the run lose more
#         than the others, in particular than participants without a
#         discernible response? Classes follow the synthesis model's
#         personal slopes (95% interval credibly < 0 = "against",
#         credibly > 0 = "with", else "uncertain"), computed here from the
#         K subsampled draws (so they can differ marginally from the
#         classes in synthesis_stan_*_slope_participants_<tr>.csv, which
#         use all draws). Per draw and per Bayesian-bootstrap weight
#         vector, the weighted mean and weighted median welfare loss (pp
#         of the endowment) of the two groups and their difference are
#         computed; cells report posterior median and 95% credible
#         interval. Descriptive, no model. Besides welfare loss, betting
#         rate and |stake deviation|, the same contrasts are run on raw
#         descriptives from master_sequences.csv: within-person SD of
#         positive stakes (stake_sd), share of focal stakes 50/100 (focal),
#         mean positive stake (stake_mean), median decision time in
#         seconds (rt_s), and individual CRRA risk aversion r_i (r_mean,
#         from mpl_scored_<tr>.csv); participants with < 2 bets have NA stake_sd.
#
#   PART C -- WHO vs WHAT (runs before B; needs only the RQ1-RQ4 fits)
#     For RQ1, RQ3 and RQ4 (RQ2 excluded: standardised within participant,
#     so its participant variance is zero by construction), the posterior
#     of the participant- and
#     sequence-level random-intercept SDs, their ratio, and the participant
#     share sigma_u^2 / (sigma_u^2 + sigma_s^2). Scales differ across models
#     (logit, standardised, log), so only the within-model ratio and share
#     are comparable across outcomes. Tests the closing claim that the
#     participant drives cost while the sequence drives direction.
#     -> LaTeX table.
#
# DRAW PAIRING
#   Each model is estimated separately, so posteriors are independent
#   across models. n_draw draws are sampled independently from each fit
#   and paired by row index; within a model, draws stay joint across
#   sequences / participants.
#
# ATTENUATION NOTE
#   Participant effects from separately estimated models are shrunk
#   independently, which biases cross-model correlations toward zero.
#   Part B estimates are therefore conservative in magnitude.
#
# INPUT (path_mod unless noted)
#   rq1_fit_sequences_<tr>_<tag><sfx>.rds  (+ _pid_levels_, _seq_levels_)
#   rq2_fit_sequences_<tr>_<tag><sfx>.rds  (+ levels, rq2_prepared_*)
#   rq3_fit_sequences_<tr>_<tag><sfx>.rds  (+ levels)
#   rq4_fit_sequences_<tr>_full.rds        (+ levels)
#   synthesis_fit_<tr>.rds
#   path_out/rq{1,2,3}_diagnostics.csv     (selected model per tr/tag)
#   path_src/master_sequences.csv          (number of bets per participant)
#   <tag> = confirmatory where betting is normative (FN), else full.
#
# OUTPUT
#   path_out/synthesis_costdir_seq_link_<tr>.csv
#   path_out/synthesis_costdir_seq_features_<tr>.csv
#   path_out/synthesis_costdir_pid_spearman_<tr>.csv
#   path_out/synthesis_costdir_pid_ols_<tr>.csv
#   path_out/synthesis_costdir_pid_groups_<tr>.csv   -- NEW (B3): welfare loss by slope class,
#                                                       group means/medians and differences (pp)
#   path_tab/synthesis_features_<tr>.tex
#   path_tab/synthesis_costdir_pid_<tr>.tex
#   path_out/synthesis_costdir_varcomp_<tr>.csv
#   path_tab/synthesis_costdir_varcomp_<tr>.tex
#
# CONFIG (all optional): cfg$synthesis$n_draw (default 1000),
#   cfg$synthesis$min_bets (default 30)
# ============================================================
synthesis_cost_vs_direction <- function(cfg) {
  
  `%||%` <- function(a, b) if (is.null(a)) b else a
  
  tr_vec   <- unique(as.character(cfg$run$treatment))
  scfg     <- cfg$synthesis %||% list()
  K        <- as.integer(scfg$n_draw   %||% 1000L)
  min_bets <- as.integer(scfg$min_bets %||% 30L)
  set.seed(as.integer(cfg$run$seed %||% 12345L))
  
  # ---------------- generic helpers ----------------
  ext <- function(fit, pars) {            # plain lists pass through (testing)
    if (inherits(fit, "stanfit")) rstan::extract(fit, pars = pars) else fit[pars]
  }
  pick <- function(n) sample.int(n, K, replace = n < K)
  # rstan returns scalar parameters as 1-d arrays -> coerce to plain vectors
  rows <- function(x, idx) if (length(dim(x)) <= 1L) as.vector(x)[idx] else
    if (length(dim(x)) == 2L) x[idx, , drop = FALSE] else x[idx, , , drop = FALSE]
  qs <- function(x) c(median = median(x, na.rm = TRUE),
                      q025 = unname(quantile(x, .025, na.rm = TRUE)),
                      q975 = unname(quantile(x, .975, na.rm = TRUE)),
                      p_gt0 = mean(x > 0, na.rm = TRUE))
  sumrow <- function(.v, ...) data.table(..., t(qs(.v)))
  save_csv <- function(dt, f, label) {
    if (!should_skip(f, cfg, "output", label)) { fwrite(dt, f); msg("Saved: ", f) }
  }
  
  n_H <- function(s) nchar(gsub("T", "", s))
  seq_features <- function(s) data.table(
    sequence  = s,
    imb_class = factor(abs(n_H(s) - 3), levels = 0:3, labels = c("3-3", "2-4", "1-5", "0-6")),
    run       = as.integer(substr(s, 5, 5) == substr(s, 6, 6)))
  
  wcor <- function(x, y, w) {             # weighted Pearson
    mx <- sum(w * x); my <- sum(w * y)
    sum(w * (x - mx) * (y - my)) / sqrt(sum(w * (x - mx)^2) * sum(w * (y - my)^2))
  }
  wspear <- function(x, y, w) wcor(rank(x), rank(y), w)
  dirichlet <- function(n) { g <- rgamma(n, 1); g / sum(g) }
  wmean <- function(x, w) sum(w * x) / sum(w)           # weighted mean
  wmed  <- function(x, w) {                             # weighted median
    o <- order(x); cw <- cumsum(w[o]) / sum(w); unname(x[o][which(cw >= 0.5)[1]])
  }
  
  # ---------------- model loaders (return K subsampled draws) ----------------
  sel_suffix <- function(rq, tr, tg) {
    f <- file.path(path_out, paste0(rq, "_diagnostics.csv"))
    if (!file.exists(f)) return("")
    d <- fread(f); r <- d[d$treatment == tr & d$tag == tg]
    if (nrow(r) == 0L) return("")
    switch(as.character(r$selected_model[1]),
           alternative = "_alt", gamma_only = "_gamma", "")
  }
  rq_files <- function(rq, tr, tg, sfx) {
    stem <- paste0("_", tr, "_", tg, sfx, ".rds")
    list(fit  = file.path(path_mod, paste0(rq, "_fit_sequences", stem)),
         pid  = file.path(path_mod, paste0(rq, "_pid_levels",    stem)),
         seq  = file.path(path_mod, paste0(rq, "_seq_levels",    stem)),
         prep = file.path(path_mod, paste0(rq, "_prepared",      stem)))
  }
  named <- function(m, nm) { colnames(m) <- nm; m }
  
  load_rq13 <- function(rq, tr, tg) {      # rq1: betting, rq3: welfare loss
    f <- rq_files(rq, tr, tg, sel_suffix(rq, tr, tg))
    if (!all(file.exists(c(f$fit, f$pid, f$seq)))) return(NULL)
    pids <- as.character(readRDS(f$pid)); seqs <- as.character(readRDS(f$seq))
    v <- if (rq == "rq1") c("mu_b", "mu_b_i") else c("mu_c", "mu_c_i")
    p <- ext(readRDS(f$fit), v); idx <- pick(nrow(p[[1]]))
    list(seq = named(rows(p[[1]], idx), seqs), pid = named(rows(p[[2]], idx), pids))
  }
  load_rq2 <- function(tr, tg) {           # stake deviation, as in 23_rq2_tables.R
    f <- rq_files("rq2", tr, tg, sel_suffix("rq2", tr, tg))
    if (!all(file.exists(unlist(f)))) return(NULL)
    pids <- as.character(readRDS(f$pid)); seqs <- as.character(readRDS(f$seq))
    e    <- as.numeric(cfg$design$seq$endowment)
    prep <- as.data.table(readRDS(f$prep))[, pid := as.character(pid)]
    pm   <- prep[, .(delta_bar = delta_bar[1], sd_star = sd_star[1]), by = pid]
    pm   <- pm[match(pids, pm$pid)]
    p    <- ext(readRDS(f$fit), c("alpha", "u", "b")); idx <- pick(length(p$alpha))
    eta  <- sweep(rows(p$u, idx), 1, rows(p$alpha, idx), "+")                 # K x N
    b    <- rows(p$b, idx)                                               # K x S
    map  <- function(x) sweep(sweep(x, 2, pm$sd_star, "*"), 2, pm$delta_bar, "+") / e
    mu_s <- sapply(seq_along(seqs), function(s) rowMeans(map(eta + b[, s])))
    list(seq = named(mu_s, seqs), pid = named(map(eta), pids))
  }
  load_rq4 <- function(tr) {               # heads share, full sample only
    f <- rq_files("rq4", tr, "full", "")
    if (!all(file.exists(c(f$fit, f$pid, f$seq)))) return(NULL)
    pids <- as.character(readRDS(f$pid)); seqs <- as.character(readRDS(f$seq))
    p <- ext(readRDS(f$fit), c("mu_h", "hbar", "mu_h_i")); idx <- pick(length(p$hbar))
    list(seq = named(rows(p$mu_h, idx), seqs), hbar = rows(p$hbar, idx),
         pid = named(rows(p$mu_h_i, idx), pids))
  }
  load_syn <- function(tr) {               # personal side-choice slopes
    f <- file.path(path_mod, paste0("synthesis_fit_", tr, ".rds"))
    if (!file.exists(f)) return(NULL)
    obj <- readRDS(f); p <- ext(obj$fit, c("beta_side", "u"))
    idx <- pick(nrow(p$beta_side)); u <- rows(p$u, idx); bs <- rows(p$beta_side, idx)
    list(run = named(bs[, 2] + u[, 5, ], obj$pid_levels),
         imb = named(bs[, 3] + u[, 6, ], obj$pid_levels))
  }
  
  # ---------------- LaTeX writer ----------------
  fmt <- function(x) { x <- ifelse(abs(x) < 0.005, 0, x); sprintf("%.2f", x) }
  cell <- function(r, digits_fmt = fmt, stars = TRUE) {
    star <- if (stars && (r$q025 > 0 || r$q975 < 0)) "^{*}" else ""
    c(sprintf("$%s%s$", digits_fmt(r$median), star),
      sprintf("$[%s, %s]$", digits_fmt(r$q025), digits_fmt(r$q975)))
  }
  write_tex <- function(f, caption, label, colspec, header, body, footer, note, what) {
    if (should_skip(f, cfg, "output", what)) return(invisible(NULL))
    writeLines(c("\\begin{table}[htbp]", "\\centering", paste0("\\caption{", caption, "}"),
                 paste0("\\label{", label, "}"), "\\begin{threeparttable}",
                 paste0("\\begin{tabular}{", colspec, "}"), "\\toprule", header, "\\midrule",
                 body, "\\midrule", footer, "\\bottomrule", "\\end{tabular}",
                 "\\begin{tablenotes}[flushleft]", "\\footnotesize",
                 paste0("\\item \\textit{Notes:} ", note), "\\end{tablenotes}",
                 "\\end{threeparttable}", "\\end{table}"), f)
    msg("Saved: ", f)
  }
  
  master <- fread(file.path(path_src, "master_sequences.csv"), encoding = "UTF-8")
  master[, `:=`(pid = as.character(pid), treat = as.character(treat),
                stake = fifelse(is.na(as.numeric(stake)), 0, as.numeric(stake)))]
  
  out <- list()
  
  for (tr in tr_vec) {
    
    conf_ok <- isTRUE(cfg$design$a_flags$betting_normative[[tr]])
    tg_main <- if (conf_ok) "confirmatory" else "full"
    tgs     <- if (conf_ok) c("confirmatory", "full") else "full"
    
    L <- list(rq4 = load_rq4(tr), syn = load_syn(tr))
    for (tg in tgs) L[[tg]] <- list(rq1 = load_rq13("rq1", tr, tg), rq2 = load_rq2(tr, tg),
                                    rq3 = load_rq13("rq3", tr, tg))
    M <- L[[tg_main]]
    if (is.null(L$rq4) || any(vapply(M, is.null, TRUE))) next
    
    # ======================================================
    # PART A -- SEQUENCE LEVEL
    # ======================================================
    seqs <- sort(colnames(L$rq4$seq))
    Fe   <- seq_features(seqs)
    bias  <- 100 * abs(L$rq4$seq[, seqs] - L$rq4$hbar)       # K x 64, pp
    bet   <- 100 * M$rq1$seq[, seqs]
    stake <- 100 * M$rq2$seq[, seqs]
    loss  <- 100 * M$rq3$seq[, seqs]
    nopure <- Fe$imb_class != "0-6"
    
    X_imb    <- model.matrix(~ imb_class, Fe)
    X_imbrun <- model.matrix(~ imb_class + run, Fe)
    slope <- function(Xbase, k) {
      X <- cbind(Xbase, bias = bias[k, ]); qr.coef(qr(X), loss[k, ])[["bias"]]
    }
    a1_draws <- data.table(
      rho_all    = sapply(1:K, function(k) cor(bias[k, ], loss[k, ], method = "spearman")),
      rho_nopure = sapply(1:K, function(k) cor(bias[k, nopure], loss[k, nopure], method = "spearman")),
      b_alone    = sapply(1:K, function(k) slope(matrix(1, 64, 1, dimnames = list(NULL, "(Intercept)")), k)),
      b_imb      = sapply(1:K, function(k) slope(X_imb, k)),
      b_imbrun   = sapply(1:K, function(k) slope(X_imbrun, k)))
    a1 <- rbindlist(list(
      sumrow(a1_draws$rho_all,    treatment = tr, statistic = "spearman_rho", spec = "all 64"),
      sumrow(a1_draws$rho_nopure, treatment = tr, statistic = "spearman_rho", spec = "excl. pure"),
      sumrow(a1_draws$b_alone,    treatment = tr, statistic = "ols_slope_loss_on_bias", spec = "alone"),
      sumrow(a1_draws$b_imb,      treatment = tr, statistic = "ols_slope_loss_on_bias", spec = "+ imbalance class"),
      sumrow(a1_draws$b_imbrun,   treatment = tr, statistic = "ols_slope_loss_on_bias", spec = "+ imbalance class + run"),
      sumrow(10 * a1_draws$b_imb, treatment = tr, statistic = "loss_pp_per_10pp_bias", spec = "+ imbalance class")))
    a1[, `:=`(sample_bet_stake_loss = tg_main, n_draws = K)]
    save_csv(a1, file.path(path_out, paste0("synthesis_costdir_seq_link_", tr, ".csv")),
             paste0("Synthesis cost-direction link (", tr, ")"))
    
    X_feat <- model.matrix(~ run + imb_class, Fe)
    qrF    <- qr(X_feat)
    Y      <- list(bias = bias, bet = bet, stake = stake, loss = loss)
    a2 <- rbindlist(lapply(names(Y), function(y) {
      B <- qr.coef(qrF, t(Y[[y]]))                              # p x K
      rbindlist(lapply(rownames(B), function(t) sumrow(B[t, ], treatment = tr, outcome = y, term = t)))
    }))
    a2[, sample_bet_stake_loss := tg_main]
    save_csv(a2, file.path(path_out, paste0("synthesis_costdir_seq_features_", tr, ".csv")),
             paste0("Synthesis sequence features (", tr, ")"))
    
    terms  <- c("run", "imb_class2-4", "imb_class1-5", "imb_class0-6", "(Intercept)")
    labels <- c("Terminal run", "Imbalance 2--4", "Imbalance 1--5", "Imbalance 0--6",
                "Constant (balanced, no run)")
    body <- unlist(lapply(seq_along(terms), function(i) {
      cs <- lapply(names(Y), function(y) cell(a2[outcome == y & term == terms[i]],
                                              function(x) sprintf("%.1f", ifelse(abs(x) < 0.05, 0, x))))
      c(paste(labels[i], "&", paste(sapply(cs, `[`, 1), collapse = " & "), "\\\\"),
        paste(" &", paste(sapply(cs, `[`, 2), collapse = " & "),
              if (i < length(terms)) "\\\\[3pt]" else "\\\\"))
    }))
    write_tex(
      file.path(path_tab, paste0("synthesis_features_", tr, ".tex")),
      caption = "Sequence features and outcomes across the four margins",
      label   = paste0("Tab:synthesis-features", if (tr != "m25") paste0("-", tr)),
      colspec = "l*{4}{c}",
      header  = c(" & Direction & \\multicolumn{3}{c}{Engagement and welfare} \\\\",
                  "\\cmidrule(lr){2-2}\\cmidrule(lr){3-5}",
                  " & Side bias & Betting rate & Stake deviation & Welfare loss \\\\",
                  " & $|\\hat\\mu^h_s-\\bar h|$ & $\\hat\\mu^b_s$ & $\\hat\\mu^a_s$ & $\\hat\\mu^c_s$ \\\\"),
      body    = body,
      footer  = c(paste("Sequences & 64 & 64 & 64 & 64 \\\\"),
                  paste("Posterior draws &", paste(rep(K, 4), collapse = " & "), "\\\\")),
      note    = paste0("Each column regresses a sequence-level quantity on a terminal-run indicator ",
                       "(last two symbols identical) and imbalance-class indicators (reference: balanced 3--3), ",
                       "separately within each posterior draw. Cells report the posterior median of the coefficient ",
                       "(percentage points) with its 95\\% credible interval in brackets; an asterisk denotes an ",
                       "interval excluding zero. Side bias uses the full sample; betting rate, stake deviation, and ",
                       "welfare loss (shares of the endowment) use the ",
                       if (conf_ok) "confirmatory" else "full", " sample. ",
                       "The 0--6 class contains only the two pure sequences."),
      what = paste0("Synthesis features table (", tr, ")"))
    
    # ======================================================
    # PART C -- WHO vs WHAT: participant vs sequence variance components
    # ======================================================
    vc_rows <- list()
    add_vc <- function(outcome, component, rq, tg, u_par, s_par) {
      sfx <- if (rq == "rq4") "" else sel_suffix(rq, tr, tg)
      f   <- rq_files(rq, tr, tg, sfx)$fit
      if (!file.exists(f)) return(invisible(NULL))
      p <- ext(readRDS(f), c(u_par, s_par))
      if (any(vapply(p, is.null, TRUE))) return(invisible(NULL))
      su <- as.vector(p[[u_par]]); ss <- as.vector(p[[s_par]])
      idx <- pick(length(su)); su <- su[idx]; ss <- ss[idx]
      d <- list(sd_participant = su, sd_sequence = ss, ratio = su / ss,
                share_participant = su^2 / (su^2 + ss^2))
      vc_rows[[length(vc_rows) + 1]] <<- rbindlist(lapply(names(d), function(q)
        sumrow(d[[q]], treatment = tr, outcome = outcome, component = component,
               model = paste0(rq, sfx), sample = tg, quantity = q)))
    }
    add_vc("Betting", "logit", "rq1", tg_main, "sigma_u", "sigma_s")
    # RQ2 (stake deviation) is deliberately excluded: its outcome is
    # standardised within participant around the participant's own mean
    # (delta_bar, sd_star; see 23_rq2_tables.R), so between-participant
    # variance is removed before estimation and sigma_u ~ 0 by construction.
    switch(sel_suffix("rq3", tr, tg_main),
           "_alt"   = add_vc("Welfare loss", "Gaussian", "rq3", tg_main, "sigma_u", "sigma_s"),
           "_gamma" = add_vc("Welfare loss", "positive-loss mean", "rq3", tg_main, "sup", "sbp"),
           { add_vc("Welfare loss", "zero-loss hurdle",  "rq3", tg_main, "su0", "sb0")
             add_vc("Welfare loss", "positive-loss mean", "rq3", tg_main, "sup", "sbp") })
    add_vc("Side choice", "logit", "rq4", "full", "sigma_u", "sigma_s")
    vc <- rbindlist(vc_rows)
    if (nrow(vc)) {
      vc[, n_draws := K]
      save_csv(vc, file.path(path_out, paste0("synthesis_costdir_varcomp_", tr, ".csv")),
               paste0("Synthesis variance components (", tr, ")"))
      
      rowkeys <- unique(vc[, .(outcome, component, sample)])
      body <- unlist(lapply(seq_len(nrow(rowkeys)), function(i) {
        k  <- rowkeys[i]
        cs <- lapply(c("sd_participant", "sd_sequence", "ratio", "share_participant"), function(q)
          cell(vc[outcome == k$outcome & component == k$component & quantity == q], stars = FALSE))
        lab <- paste0(k$outcome, ", ", k$component)
        c(paste(lab, "&", paste(sapply(cs, `[`, 1), collapse = " & "), "\\\\"),
          paste(" &", paste(sapply(cs, `[`, 2), collapse = " & "),
                if (i < nrow(rowkeys)) "\\\\[3pt]" else "\\\\"))
      }))
      write_tex(
        file.path(path_tab, paste0("synthesis_costdir_varcomp_", tr, ".tex")),
        caption = "Participant versus sequence heterogeneity, by outcome",
        label   = paste0("Tab:synthesis-varcomp", if (tr != "m25") paste0("-", tr)),
        colspec = "l*{4}{c}",
        header  = c(" & $\\sigma_{\\text{participant}}$ & $\\sigma_{\\text{sequence}}$ & Ratio & Participant share \\\\"),
        body    = body,
        footer  = paste("Posterior draws & \\multicolumn{4}{c}{", K, "} \\\\"),
        note    = paste0("Standard deviations of the participant- and sequence-level random intercepts in the ",
                         "preregistered model for each outcome (RQ1 betting, RQ3 welfare loss, RQ4 side choice), on each model's linear-predictor scale. Because scales differ across ",
                         "models, only the within-model ratio $\\sigma_{\\text{participant}}/\\sigma_{\\text{sequence}}$ ",
                         "and the participant share $\\sigma^2_{\\text{participant}}/(\\sigma^2_{\\text{participant}}+",
                         "\\sigma^2_{\\text{sequence}})$ are comparable across rows. Cells report posterior medians ",
                         "and 95\\% credible intervals. Side choice uses the full sample; the other outcomes the ",
                         if (conf_ok) "confirmatory" else "full", " sample."),
        what = paste0("Synthesis variance components table (", tr, ")"))
    }
    
    # ======================================================
    # PART B -- PARTICIPANT LEVEL
    # ======================================================
    if (is.null(L$syn)) next
    nb_all <- master[treat == tr, .(n_bets = sum(stake > 0)), by = pid]
    
    b_sp <- list(); b_ols <- list(); b_grp <- list()
    for (tg in tgs) {
      S_ <- L[[tg]]
      if (any(vapply(S_, is.null, TRUE))) next
      P <- Reduce(intersect, list(colnames(L$syn$run), colnames(L$rq4$pid),
                                  colnames(S_$rq1$pid), colnames(S_$rq2$pid), colnames(S_$rq3$pid)))
      P <- sort(P); n <- length(P)
      run_r <- -L$syn$run[, P]
      imb_r <- abs(L$syn$imb[, P])
      lossP <- S_$rq3$pid[, P]
      betP  <- S_$rq1$pid[, P]
      stkP  <- abs(S_$rq2$pid[, P])
      sideP <- abs(L$rq4$pid[, P] - L$rq4$hbar)
      nbets <- nb_all[match(P, pid), n_bets]
      keepB <- nbets >= min_bets
      
      W <- t(replicate(K, dirichlet(n)))                        # K x n
      U <- matrix(1 / n, K, n)
      
      for (unc in c("posterior+bootstrap", "posterior")) {
        Wm <- if (unc == "posterior") U else W
        for (x in c("run_resp", "imb_resp")) {
          Xm <- if (x == "run_resp") run_r else imb_r
          for (y in c("loss", "bet", "abs_stake")) {
            Ym <- switch(y, loss = lossP, bet = betP, abs_stake = stkP)
            r  <- sapply(1:K, function(k) wspear(Xm[k, ], Ym[k, ], Wm[k, ]))
            b_sp[[length(b_sp) + 1]] <- sumrow(r, treatment = tr, sample = tg, subset = "all",
                                               uncertainty = unc, x = x, y = y, n = n)
          }
          r <- sapply(1:K, function(k) {
            w <- Wm[k, keepB]; wspear(Xm[k, keepB], lossP[k, keepB], w / sum(w)) })
          b_sp[[length(b_sp) + 1]] <- sumrow(r, treatment = tr, sample = tg,
                                             subset = paste0("n_bets >= ", min_bets), uncertainty = unc,
                                             x = x, y = "loss", n = sum(keepB))
        }
      }
      
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
        b_ols[[length(b_ols) + 1]] <- rbindlist(lapply(rownames(cf), function(t)
          sumrow(cf[t, ], treatment = tr, sample = tg, spec = sp, term = t, n = n)))
      }
      
      # ---- B3: welfare loss, betting rate, and stake deviation by slope class ----
      # Class of each participant from the 95% interval of their personal
      # side-choice slope across the K subsampled draws: credibly > 0
      # ("with"), credibly < 0 ("against"), else "uncertain". For the
      # imbalance slope, "against" = counters the majority (balance-seeking).
      slope_class <- function(M) {               # M: K x P signed slope draws
        lo <- apply(M, 2, quantile, probs = .025)
        hi <- apply(M, 2, quantile, probs = .975)
        ifelse(lo > 0, "with", ifelse(hi < 0, "against", "uncertain"))
      }
      run_cls <- slope_class(L$syn$run[, P])
      imb_cls <- slope_class(L$syn$imb[, P])
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
      # outcomes: welfare loss, betting rate, absolute stake deviation (pp; as in
      # the Spearman block), plus raw descriptives added below
      b_out <- list(loss = 100 * lossP, bet = 100 * betP, abs_stake = 100 * stkP)
      # raw participant-level descriptives (no posterior; identical across draws):
      # within-person SD of positive stakes, share of focal stakes (50 or 100),
      # mean positive stake, and median decision time (s). NA when < 2 bets.
      raw <- master[treat == tr & pid %in% P, {
        st <- stake[stake > 0]
        .(stake_sd   = if (length(st) >= 2) sd(st) else NA_real_,
          focal      = if (length(st) >= 1) 100 * mean(st %in% c(50, 100)) else NA_real_,
          stake_mean = if (length(st) >= 1) mean(st) else NA_real_,
          rt_s       = median(screen_ms, na.rm = TRUE) / 1000)
      }, by = pid][match(P, pid)]
      for (v in c("stake_sd", "focal", "stake_mean", "rt_s"))
        b_out[[v]] <- matrix(raw[[v]], K, n, byrow = TRUE)
      # individual CRRA risk aversion r_i (posterior mean from the MPL model,
      # 02_stan_r_from_mpl.R); NA for participants without an MPL estimate
      f_mpl <- file.path(path_out, paste0("mpl_scored_", tr, ".csv"))
      if (file.exists(f_mpl)) {
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
            b_grp[[length(b_grp) + 1]] <- rbindlist(lapply(rownames(res), function(q)
              sumrow(res[q, ], treatment = tr, sample = tg, uncertainty = unc, outcome = oc,
                     contrast = ct$lab, quantity = q, n_a = sum(ct$A & !is.na(Ym[1, ])),
                     n_b = sum(ct$B & !is.na(Ym[1, ])), n = n)))
          }
        }
      }
    }
    b_sp <- rbindlist(b_sp); b_ols <- rbindlist(b_ols); b_grp <- rbindlist(b_grp)
    b_sp[, n_draws := K]; b_ols[, `:=`(n_draws = K, uncertainty = "posterior+bootstrap")]
    save_csv(b_sp,  file.path(path_out, paste0("synthesis_costdir_pid_spearman_", tr, ".csv")),
             paste0("Synthesis participant Spearman (", tr, ")"))
    save_csv(b_ols, file.path(path_out, paste0("synthesis_costdir_pid_ols_", tr, ".csv")),
             paste0("Synthesis participant OLS (", tr, ")"))
    if (nrow(b_grp)) {
      b_grp[, n_draws := K]
      save_csv(b_grp, file.path(path_out, paste0("synthesis_costdir_pid_groups_", tr, ".csv")),
               paste0("Synthesis participant loss by slope class (", tr, ")"))
    }
    
    # ---- participant LaTeX table (main sample) ----
    O <- b_ols[sample == tg_main]
    if (nrow(O)) {
      sp_names <- c("base", "n_bets", "n_bets_side", "n_bets_side_b")
      tlab <- c(run_resp = "Run responsiveness", imb_resp = "Imbalance responsiveness",
                n_bets = "Number of bets", side_pref = "Side preference $|\\hat\\mu^h_i-\\bar h|$",
                bet = "Betting rate $\\hat\\mu^b_i$")
      body <- unlist(lapply(names(tlab), function(t) {
        cs <- lapply(sp_names, function(s) {
          r <- O[spec == s & term == t]
          if (nrow(r)) cell(r, function(x) sprintf("%.3f", ifelse(abs(x) < 5e-4, 0, x))) else c("", "")
        })
        c(paste(tlab[[t]], "&", paste(sapply(cs, `[`, 1), collapse = " & "), "\\\\"),
          paste(" &", paste(sapply(cs, `[`, 2), collapse = " & "), "\\\\[3pt]"))
      }))
      Sp <- b_sp[sample == tg_main & subset == "all" & y == "loss" & uncertainty == "posterior+bootstrap"]
      rho_line <- function(xx, lab) {
        r <- Sp[Sp$x == xx]
        paste0(lab, " & \\multicolumn{4}{c}{$", fmt(r$median), "$ \\; $[", fmt(r$q025), ", ",
               fmt(r$q975), "]$} \\\\")
      }
      write_tex(
        file.path(path_tab, paste0("synthesis_costdir_pid_", tr, ".tex")),
        caption = "Participant-level cue responsiveness and welfare loss",
        label   = paste0("Tab:synthesis-costdir-pid", if (tr != "m25") paste0("-", tr)),
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
                         if (conf_ok) "Confirmatory" else "Full", " sample."),
        what = paste0("Synthesis participant table (", tr, ")"))
    }
    
    out[[tr]] <- list(seq_link = a1, seq_features = a2, varcomp = vc, pid_spearman = b_sp, pid_ols = b_ols,
                      pid_groups = b_grp)
  }
  
  invisible(out)
}