# ============================================================
# _rq3-welfare-loss.R
#
# PURPOSE
#   Computes the trial-level welfare loss of RQ3 for one treatment and
#   sample, shared by the three RQ3 likelihoods and by the floor
#   robustness check. The loss of a bet a on trial t of participant i is
#       y_t = E_k[ max( CE(a*_ik) - CE(a_t), 0 ) ] / e
#   where CE(.) is the certainty equivalent of a stake under CRRA utility
#   with risk parameter r_ik, a*_ik is the expected-utility optimal stake
#   for the same draw, and k indexes posterior draws. The loss is expressed
#   as a proportion of the endowment e. A trial with y_t = 0 is an
#   expected-utility optimal choice in every sampled draw.
#
# FUNCTION
#   rq3_welfare_loss(cfg, tr, sample, xmin = NULL)
#     tr      treatment ("m25", "m19")
#     sample  "full" or "confirmatory"
#     xmin    wealth floor at the zero-wealth state of an all-in stake;
#             defaults to cfg$design$seq$xmin
#   Returns list(d, y, pid_levels, seq_levels, krep) with d the trial table
#   (pid, seq, block_c, stake, pid_i, sid_s) in the order of y.
#
# INPUT
#   model "mpl": r-draws-<tr>.rds        list(pid, r_draws), draws x pid
#                a-star/a-star-draws-<tr>-full.rds   list(pid, a_star_draws)
#   analysis_trials(tr, "rq3-<sample>")  (09_analysis-sample.R)
#   ce_stake_vec_r()                     (scripts/config/02_eu-functions.R)
#   cfg$design$seq (endowment, xmin, coin_prob, treatments), cfg$model$ppc$rq3_k
#
# NOTES
#   - A random subset of cfg$model$ppc$rq3_k posterior draws (seeded with
#     cfg$run$seed) is used; the same draw index is used for r and a*, so a*
#     is the optimum for the r of that draw
#   - r and a* come from the full-sample fits: the risk parameters of the
#     Holt-Laury model are not pooled across participants, so the draws of a
#     participant do not depend on the sample
#   - The optimal certainty equivalent depends on the participant alone and
#     the actual one on the participant and the stake; each is evaluated
#     once per distinct combination
# ============================================================

rq3_welfare_loss <- function(cfg, tr, sample, xmin = NULL) {

  design <- cfg$design
  e      <- as.numeric(design$seq$endowment)
  p_win  <- as.numeric(design$seq$coin_prob)
  if (is.null(xmin)) xmin <- as.numeric(design$seq$xmin)
  m_map  <- design$seq$treatments
  stopifnot(tr %in% names(m_map))
  m <- as.numeric(m_map[[tr]])
  stopifnot(length(e) == 1L, is.finite(e), e > 0,
            length(xmin) == 1L, is.finite(xmin), xmin > 0,
            length(p_win) == 1L, p_win >= 0, p_win <= 1,
            is.finite(m), m > 1)

  seed <- as.integer(cfg$run$seed)
  krep <- as.integer(cfg$model$ppc$rq3_k)
  stopifnot(length(krep) == 1L, krep >= 10L)

  # ---- Posterior draws of r and a* ----
  r_obj <- readRDS(model_path("mpl", paste0("r-draws-", tr, ".rds")))
  a_obj <- readRDS(model_path("mpl", file.path("a-star", paste0("a-star-draws-", tr, "-full.rds"))))
  r_all <- r_obj$r_draws
  a_all <- a_obj$a_star_draws
  K_all <- nrow(r_all)
  stopifnot(K_all >= 10L, nrow(a_all) == K_all)

  krep <- min(krep, K_all)
  set.seed(seed)
  k_idx <- sort(sample.int(K_all, krep, replace = FALSE))

  # ---- Trials of the sample ----
  d <- analysis_trials(tr, paste0("rq3-", sample))[, .(pid, seq, stake, block_c)]
  d[, pid := as.character(pid)]
  stopifnot(!anyNA(d$stake), all(d$stake >= 0 & d$stake <= e))

  pid_levels <- sort(unique(d$pid))
  seq_levels <- sort(unique(d$seq))
  idx_r <- match(pid_levels, as.character(r_obj$pid))
  idx_a <- match(pid_levels, as.character(a_obj$pid))
  if (anyNA(idx_r)) stop("RQ3: participants missing from the r draws (", tr, ").", call. = FALSE)
  if (anyNA(idx_a)) stop("RQ3: participants missing from the a* draws (", tr, ").", call. = FALSE)

  r_draws <- r_all[k_idx, idx_r, drop = FALSE]
  a_draws <- a_all[k_idx, idx_a, drop = FALSE]

  d[, pid_i := match(pid, pid_levels)]
  d[, sid_s := match(seq, seq_levels)]

  # ---- Optimal certainty equivalent: once per participant ----
  ce_opt <- lapply(seq_along(pid_levels), function(i)
    ce_stake_vec_r(a = a_draws[, i], r = r_draws[, i], m = m, e = e,
                   p_win = p_win, xmin = xmin))

  # ---- Welfare loss: once per distinct (participant, stake) ----
  cells <- unique(d[, .(pid_i, stake)])
  cells[, loss := vapply(seq_len(.N), function(j) {
    i <- pid_i[j]
    ce_act <- ce_stake_vec_r(a = rep(stake[j], krep), r = r_draws[, i], m = m, e = e,
                             p_win = p_win, xmin = xmin)
    mean(pmax(ce_opt[[i]] - ce_act, 0)) / e
  }, numeric(1))]

  y <- cells$loss[match(paste(d$pid_i, d$stake), paste(cells$pid_i, cells$stake))]
  if (!all(is.finite(y))) stop("RQ3: non-finite welfare loss (", tr, "/", sample, ").", call. = FALSE)

  list(d = d, y = y, pid_levels = pid_levels, seq_levels = seq_levels, krep = krep)
}
