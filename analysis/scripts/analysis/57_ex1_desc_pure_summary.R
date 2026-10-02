# ============================================================
# 57_ex1_desc_pure_summary.R
#
# EX1.III fallback: descriptive associations between pure-sequence
# classification (hot-hand consistent / gambler consistent /
# inconsistent / missing) and the three EX1.II covariates:
#   Zopt  standardized LOT-R (optimism)
#   Zrt   standardized log mean RT on betting trials
#   Zr    standardized posterior mean risk parameter r
#
# For each treatment x covariate x type cell, reports:
#   n, mean, SD, 95% bootstrap CI (2000 resamples)
#
# INPUTS
#   path_out/ex1_pure_seq_summary.csv   (from ex1_pure_seq_summary)
#   path_out/lotr_scored.csv
#   path_mod/mpl_r_draws_<tr>.rds
#   path_src/master_sequences.csv
#
# OUTPUTS
#   path_out/ex1_pure_seq_covariates.csv
#   path_fig/ex1_pure_seq_covariates.pdf
#   path_fig/ex1_pure_seq_covariates.png
# ============================================================

ex1_pure_summary <- function(cfg) {
  
  lab_heads  <- as.character(cfg$design$seq$side_labels$heads)
  lab_tails  <- as.character(cfg$design$seq$side_labels$tails)
  pure_H     <- as.character(cfg$design$seq$anchor_labels$pure_heads)
  pure_T     <- as.character(cfg$design$seq$anchor_labels$pure_tails)
  tr_vec     <- unique(as.character(cfg$run$treatment))
  
  dt <- fread(file.path(path_src, "master_sequences.csv"), encoding = "UTF-8")
  dt[, pid   := as.character(pid)]
  dt[, treat := as.character(treat)]
  dt[, seq   := as.character(seq)]
  dt[, stake := as.numeric(stake)]
  dt[, side  := as.character(side)]
  dt[is.na(stake), stake := 0]
  
  d <- dt[
    seq %in% c(pure_H, pure_T) &
      treat %in% tr_vec &
      is.finite(stake) & stake > 0 &
      side %in% c(lab_heads, lab_tails)
  ]
  
  d[, h := as.integer(side == lab_heads)]
  
  d_agg <- d[, .(
    n_trials  = .N,
    n_heads   = sum(h),
    p_heads   = mean(h),
    side_mode = ifelse(mean(h) > 0.5, lab_heads,
                       ifelse(mean(h) < 0.5, lab_tails, NA_character_))
  ), by = .(treat, pid, seq)]
  
  d_H <- d_agg[seq == pure_H, .(treat, pid,
                                h_pureH = n_heads >= 1L,
                                n_pureH = n_trials)]
  
  d_T <- d_agg[seq == pure_T, .(treat, pid,
                                h_pureT = n_heads >= 1L,
                                n_pureT = n_trials)]
  
  pid_tbl <- merge(d_H, d_T, by = c("treat", "pid"), all = TRUE)
  
  pid_tbl[, chi_pure := fcase(
    h_pureH == TRUE  & h_pureT == FALSE,  1L,
    h_pureH == FALSE & h_pureT == TRUE,  -1L,
    h_pureH == TRUE  & h_pureT == TRUE,   0L,
    h_pureH == FALSE & h_pureT == FALSE,  0L,
    default = NA_integer_
  )]
  
  pid_tbl[, type := fcase(
    chi_pure ==  1L, "hothand_consistent",
    chi_pure == -1L, "gambler_consistent",
    chi_pure ==  0L, "inconsistent",
    is.na(chi_pure), "missing_pure_seq"
  )]
  
  setorder(pid_tbl, treat, pid)
  
  counts <- pid_tbl[, .(n = .N), by = .(treat, type)]
  counts[, pct := round(100 * n / sum(n), 1), by = treat]
  counts[, type := factor(type, levels = c("hothand_consistent",
                                           "gambler_consistent",
                                           "inconsistent",
                                           "missing_pure_seq"))]
  setorder(counts, treat, type)
  
  msg("\n=== Pure sequence summary ===")
  for (tr in tr_vec) {
    msg("\nTreatment: ", tr)
    print(counts[treat == tr])
  }
  
  f_pid <- file.path(path_out, "ex1_pure_seq_summary.csv")
  f_cnt <- file.path(path_out, "ex1_pure_seq_counts.csv")
  fwrite(pid_tbl, f_pid)
  fwrite(counts,  f_cnt)
  msg("\nSaved: ", f_pid)
  msg("Saved: ", f_cnt)
  
  # ==========================================================
  # COVARIATE ASSOCIATIONS
  # ==========================================================
  
  set.seed(as.integer(cfg$run$seed))
  n_boot <- 2000L
  
  f_lotr <- file.path(path_out, "lotr_scored.csv")
  if (!file.exists(f_lotr))
    stop("ex1_pure_seq_summary: missing ", f_lotr)
  
  lotr <- fread(f_lotr, encoding = "UTF-8")
  lotr[, pid        := as.character(pid)]
  lotr[, lotr_score := as.numeric(lotr_score)]
  lotr  <- lotr[, .(pid, lotr_score)]
  
  master <- fread(file.path(path_src, "master_sequences.csv"),
                  encoding = "UTF-8")
  master[, pid       := as.character(pid)]
  master[, treat     := as.character(treat)]
  master[, stake     := as.numeric(stake)]
  master[is.na(stake), stake := 0]
  master[, screen_ms := as.numeric(screen_ms)]
  
  rt_pid <- master[
    treat %in% tr_vec & stake > 0 &
      is.finite(screen_ms) & screen_ms > 0,
    .(rt_log_mean = mean(log(screen_ms))),
    by = .(treat, pid)
  ]
  
  r_list <- lapply(tr_vec, function(tr) {
    f_r <- file.path(path_mod, paste0("mpl_r_draws_", tr, ".rds"))
    if (!file.exists(f_r))
      stop("ex1_pure_seq_summary: missing ", f_r)
    obj <- readRDS(f_r)
    data.table(
      treat  = tr,
      pid    = as.character(obj$pid),
      r_mean = apply(obj$r_draws, 2, mean)
    )
  })
  r_dt <- rbindlist(r_list)
  
  cov <- Reduce(
    function(a, b) merge(a, b, by = intersect(names(a), names(b)), all.x = TRUE),
    list(pid_tbl[, .(treat, pid, type)], lotr, r_dt, rt_pid)
  )
  
  z_score <- function(x) {
    m <- mean(x, na.rm = TRUE)
    s <- sd(x,   na.rm = TRUE)
    if (!is.finite(s) || s <= 0) return(rep(NA_real_, length(x)))
    (x - m) / s
  }
  
  cov[, Zopt := z_score(lotr_score),  by = treat]
  cov[, Zrt  := z_score(rt_log_mean), by = treat]
  cov[, Zr   := z_score(r_mean),      by = treat]
  
  ci_t <- function(x, alpha = 0.05) {
    x <- x[is.finite(x)]
    n <- length(x)
    if (n < 2L) return(c(lo = NA_real_, hi = NA_real_))
    se <- sd(x) / sqrt(n)
    t  <- qt(1 - alpha / 2, df = n - 1)
    setNames(mean(x) + c(-1, 1) * t * se, c("lo", "hi"))
  }
  
  covariates <- c(Zopt = "Optimism (LOT-R)",
                  Zrt  = "Log RT",
                  Zr   = "Risk aversion (r)")
  
  type_levels <- c("hothand_consistent", "gambler_consistent",
                   "inconsistent",       "missing_pure_seq")
  
  rows <- lapply(tr_vec, function(tr) {
    d_tr <- cov[treat == tr]
    lapply(names(covariates), function(cv) {
      lapply(type_levels, function(tp) {
        x  <- d_tr[type == tp, get(cv)]
        x  <- x[is.finite(x)]
        n  <- length(x)
        ci <- ci_t(x)
        data.table(
          treatment   = tr,
          covariate   = covariates[cv],
          covariate_z = cv,
          type        = tp,
          n           = n,
          mean        = if (n >= 1L) mean(x) else NA_real_,
          sd          = if (n >= 2L) sd(x)   else NA_real_,
          ci_lo       = ci["lo"],
          ci_hi       = ci["hi"]
        )
      }) |> rbindlist()
    }) |> rbindlist()
  }) |> rbindlist()
  
  rows[, type := factor(type, levels = type_levels)]
  setorder(rows, treatment, covariate, type)
  
  f_cov <- file.path(path_out, "ex1_pure_seq_covariates.csv")
  fwrite(rows, f_cov)
  msg("Saved: ", f_cov)
  
  msg("\n=== EX1 pure-sequence covariate associations ===")
  for (tr in tr_vec) {
    msg("\nTreatment: ", tr)
    print(rows[treatment == tr,
               .(covariate, type, n,
                 mean  = round(mean,  3),
                 sd    = round(sd,    3),
                 ci_lo = round(ci_lo, 3),
                 ci_hi = round(ci_hi, 3))],
          nrows = Inf)
  }
  
  # ==========================================================
  # PLOT
  # ==========================================================
  
  type_labels <- c(
    hothand_consistent = "Hot-hand consistent",
    gambler_consistent = "Gambler consistent",
    inconsistent       = "Inconsistent",
    missing_pure_seq   = "Missing"
  )
  
  type_colors <- c(
    hothand_consistent = "#185FA5",
    gambler_consistent = "#D85A30",
    inconsistent       = "#888780",
    missing_pure_seq   = "#B4B2A9"
  )
  
  treat_labels <- setNames(
    paste0(tr_vec, " (m = ",
           sapply(tr_vec, function(tr) cfg$design$seq$treatments[[tr]]), ")"),
    tr_vec
  )
  
  plot_dt <- copy(rows)
  plot_dt[, type_label := factor(type_labels[as.character(type)],
                                 levels = type_labels)]
  plot_dt[, treat_label := treat_labels[treatment]]
  
  # n label for each cell
  n_labels <- plot_dt[, .(
    n_label = paste0("n=", n),
    covariate, type_label, treat_label,
    mean, ci_lo, ci_hi
  )]
  
  p <- ggplot(plot_dt,
              aes(x      = type_label,
                  y      = mean,
                  ymin   = ci_lo,
                  ymax   = ci_hi,
                  colour = as.character(type))) +
    geom_hline(yintercept = 0, linetype = "dashed",
               colour = "grey60", linewidth = 0.4) +
    geom_errorbar(width = 0.25, linewidth = 0.6) +
    geom_point(size = 2.5) +
    geom_text(data = n_labels,
              aes(x = type_label, y = ci_lo - 0.12, label = n_label),
              size = 2.8, colour = "grey50", inherit.aes = FALSE) +
    scale_colour_manual(values = type_colors, guide = "none") +
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 12)) +
    facet_grid(covariate ~ treat_label, scales = "free_y") +
    labs(
      x       = NULL,
      y       = "Standardized mean (95% bootstrap CI)",
      caption = paste0("Covariates standardized within treatment. ",
                       "CIs based on ", n_boot, " bootstrap resamples.")
    ) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      strip.background   = element_rect(fill = "grey95", colour = "grey80"),
      strip.text         = element_text(size = 9),
      axis.text.x        = element_text(size = 8),
      plot.caption       = element_text(size = 7, colour = "grey50", hjust = 0)
    )
  

  f_png <- file.path(path_fig, "ex1_pure_seq_covariates.png")
  ggsave(f_png, p, width = 3.5 * length(tr_vec), height = 7, dpi = 300)
  msg("Saved: ", f_png)
  
  invisible(list(participants = pid_tbl, counts = counts,
                 covariates   = rows,    plot    = p))
}