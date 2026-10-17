# ============================================================
# 51_ex1-chi-draws.R
#
# PURPOSE
#   Constructs the posterior draws of the directional indices of EX1
#   (Gamblerish-Hot-handish Index, GHI) from the draws of the RQ4 model.
#   No model is fitted here: the script transforms posterior draws, and
#   the draws are stored in a model folder so that the later EX1 steps
#   (regression on the GHI, tables, figures) use the same set of draws.
#
#   For posterior draw k of the RQ4 model, theta_s = mu_h[s] is the
#   population Heads probability of sequence s and theta_0 = hbar the
#   population baseline. The pure sequences H = HHHHHH and T = TTTTTT are
#   the anchors. With tolerance eps and stabiliser eta:
#       I_a(k,s) = 1{ |theta_s - theta_a| < eps },   a in {H, T, 0}
#       w_a(k,s) = (I_a(k,s) + eta) / (I_H + I_T + I_0 + 3 eta)
#       z_is(k)  = sign_is * (w_H(k,s) - w_T(k,s))
#       chi_s(k) = mean of z_is(k) over the betting trials on sequence s
#       chi_i(k) = mean of z_is(k) over the betting trials of participant i
#   where sign_is = +1 if participant i chose Heads on sequence s and -1 if
#   Tails. The preregistration collapses the similarities across draws before
#   this step; the draw-level version is used here so that the classification
#   probabilities are well defined (preregistration addendum C.6).
#
#   Because w_H - w_T depends on the draw and the sequence only, chi_s and
#   chi_i are computed exactly by matrix products, not trial by trial:
#       chi_s(k) = (w_H - w_T)(k,s) * mean_sign_s
#       chi_i(k) = sum_s (w_H - w_T)(k,s) * A[s,i],  A[s,i] = (sum of signs of i on s) / n_i
#
# SAMPLE
#   rq4-full for each treatment in cfg$run$treatment (betting trials, side
#   H or T), from analysis_trials() (09_analysis-sample.R). The participants
#   are those of the RQ4 fit.
#
# INPUT
#   model "rq4": fit-<tr>-full.rds, pid-levels-<tr>-full.rds,
#                seq-levels-<tr>-full.rds   (mu_h, hbar)
#   model "mpl": sample/analysis-trials.csv, sample/sample-definitions.csv
#   cfg$design$ex1 (eps, eta), cfg$design$seq$anchor_labels, side_labels
#
# OUTPUT (model folder path_mod/ex1-chi-draws/)
#   chi-draws-<tr>-eps<xx>.rds   one file per treatment and tolerance:
#       list(treatment, eps, eta, seq_levels, pid_levels, n_s, n_i,
#            chi_s [draws x sequences], chi_i [draws x participants])
#   diagnostics.csv              sample sizes against the definition, number of
#                                non-finite draws, anchors present
#   inputs.csv, run-info.txt
#
# NOTES
#   - The model folder is built as ex1-chi-draws.tmp and renamed when
#     complete; an existing folder is never overwritten
#   - eps in cfg$design$ex1$eps (first value = main, others sensitivity);
#     the file name carries eps with two decimals (0.05 -> eps005)
#   - A sequence without a betting trial would have an undefined chi_s; such
#     a case is counted in the diagnostics
# ============================================================

ex1_chi_draws <- function(cfg) {

  run <- model_begin("ex1-chi-draws")
  if (is.null(run)) return(invisible(NULL))

  design <- cfg$design
  tr_vec <- unique(as.character(cfg$run$treatment))

  eps_vec <- as.numeric(design$ex1$eps)
  eta     <- as.numeric(design$ex1$eta)
  stopifnot(length(eps_vec) >= 1L, all(is.finite(eps_vec)), all(eps_vec > 0),
            length(eta) == 1L, is.finite(eta), eta > 0)

  lab_heads  <- as.character(design$seq$side_labels$heads)
  lab_tails  <- as.character(design$seq$side_labels$tails)
  pure_heads <- as.character(design$seq$anchor_labels$pure_heads)
  pure_tails <- as.character(design$seq$anchor_labels$pure_tails)

  eps_nm <- function(eps) gsub("\\.", "", sprintf("%.2f", eps))

  diag <- list()

  for (tr in tr_vec) {

    fit_file <- function(prefix)
      model_path("rq4", paste0(prefix, "-", tr, "-full.rds"))
    fit        <- readRDS(fit_file("fit"))
    seq_levels <- as.character(readRDS(fit_file("seq-levels")))
    pid_levels <- as.character(readRDS(fit_file("pid-levels")))

    post <- rstan::extract(fit, pars = c("mu_h", "hbar"))
    theta_s <- post$mu_h
    theta_0 <- as.numeric(post$hbar)
    K <- nrow(theta_s); S <- length(seq_levels); N <- length(pid_levels)
    stopifnot(ncol(theta_s) == S, length(theta_0) == K)

    idx_H <- match(pure_heads, seq_levels)
    idx_T <- match(pure_tails, seq_levels)
    if (is.na(idx_H) || is.na(idx_T))
      stop("EX1: anchor sequence missing in the RQ4 fit (", tr, ").", call. = FALSE)
    theta_H <- theta_s[, idx_H]
    theta_T <- theta_s[, idx_T]

    # ---- Betting trials with their side ----
    d <- analysis_trials(tr, "rq4-full")
    stopifnot(all(d$side %in% c(lab_heads, lab_tails)), setequal(unique(d$pid), pid_levels))
    d[, sgn   := ifelse(side == lab_heads, 1, -1)]
    d[, pid_i := match(pid, pid_levels)]
    d[, sid_s := match(seq, seq_levels)]
    stopifnot(!anyNA(d$pid_i), !anyNA(d$sid_s))

    n_s <- tabulate(d$sid_s, S)
    n_i <- tabulate(d$pid_i, N)
    sg_s <- numeric(S)
    s_tab <- d[, .(sg = sum(sgn)), by = sid_s]
    sg_s[s_tab$sid_s] <- s_tab$sg
    m_s <- ifelse(n_s > 0L, sg_s / pmax(n_s, 1L), NA_real_)   # mean sign per sequence

    # A[s, i] = (sum of signs of participant i on sequence s) / n_i
    agg <- d[, .(sg = sum(sgn)), by = .(pid_i, sid_s)]
    A <- matrix(0, nrow = S, ncol = N)
    A[cbind(agg$sid_s, agg$pid_i)] <- agg$sg / n_i[agg$pid_i]

    for (eps in eps_vec) {

      I_H <- abs(theta_s - theta_H) < eps          # recycled by column: draw k uses theta_H[k]
      I_T <- abs(theta_s - theta_T) < eps
      I_0 <- abs(theta_s - theta_0) < eps
      denom  <- I_H + I_T + I_0 + 3 * eta
      w_H    <- (I_H + eta) / denom
      w_T    <- (I_T + eta) / denom
      w_0    <- (I_0 + eta) / denom
      diff_w <- w_H - w_T

      chi_s <- sweep(diff_w, 2, m_s, "*")
      chi_i <- diff_w %*% A

      saveRDS(list(treatment = tr, eps = eps, eta = eta, seq_levels = seq_levels,
                   pid_levels = pid_levels, n_s = n_s, n_i = n_i,
                   chi_s = chi_s, chi_i = chi_i),
              file.path(run, paste0("chi-draws-", tr, "-eps", eps_nm(eps), ".rds")))

      tag <- paste(tr, "full", paste0("eps", eps_nm(eps)), sep = "_")
      diag[[tag]] <- data.table(
        tag = tag, treatment = tr, sample = "full",
        check = c("chi_s_nonfinite", "chi_i_nonfinite", "weights_sum_to_one"),
        value = c(sum(!is.finite(chi_s)), sum(!is.finite(chi_i)),
                  max(abs(w_H + w_T + w_0 - 1))),
        threshold = c(0, 0, 1e-8)
      )
      diag[[tag]][, pass := value <= threshold]
    }

    # sample size against the definition of the sample
    def <- sample_definition(tr, "rq4-full")
    tag <- paste(tr, "full", "sample", sep = "_")
    diag[[tag]] <- data.table(
      tag = tag, treatment = tr, sample = "full",
      check = c("rq4-full_participants", "rq4-full_trials"),
      value = c(N, nrow(d)), threshold = c(def$n_participants, def$n_trials)
    )
    diag[[tag]][, pass := value == threshold]

    msg("EX1 chi draws: ", tr, " | draws=", K, " | S=", S, " | N=", N)
  }

  write_diagnostics(run, rbindlist(diag))
  model_end(run, cfg, settings = list(eps = eps_vec, eta = eta), key = "ex1-chi")
}
