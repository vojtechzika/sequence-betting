# ============================================================
# 58_ex1_desc_pure_vs_covariates.R
#
# Descriptive check: does dispositional optimism (LOT-R)
# correlate with response time across participants?
#
# Computes Pearson and Spearman correlations between
# standardized LOT-R score and mean log RT on betting trials,
# separately per treatment. Reports point estimates and
# 95% bootstrap CIs.
#
# INPUTS
#   path_out/lotr_scored.csv
#   path_src/master_sequences.csv
#
# OUTPUTS
#   path_out/ex1_optimism_rt.csv
# ============================================================

ex1_pure_vs_covariates <- function(cfg) {
  
  set.seed(as.integer(cfg$run$seed))
  tr_vec <- unique(as.character(cfg$run$treatment))
  n_boot <- 2000L
  
  # ---- LOT-R ----
  f_lotr <- file.path(path_out, "lotr_scored.csv")
  if (!file.exists(f_lotr))
    stop("ex1_optimism_rt: missing ", f_lotr)
  
  lotr <- fread(f_lotr, encoding = "UTF-8")
  lotr[, pid        := as.character(pid)]
  lotr[, lotr_score := as.numeric(lotr_score)]
  lotr  <- lotr[, .(pid, lotr_score)]
  
  # ---- RT ----
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
  
  # ---- Merge ----
  d <- merge(rt_pid, lotr, by = "pid", all.x = TRUE)
  d <- d[is.finite(lotr_score) & is.finite(rt_log_mean)]
  
  # ---- Bootstrap CI helper ----
  boot_cor <- function(x, y, method, B = n_boot, alpha = 0.05) {
    n <- length(x)
    if (n < 5L) return(c(lo = NA_real_, hi = NA_real_))
    vals <- replicate(B, {
      idx <- sample.int(n, n, replace = TRUE)
      cor(x[idx], y[idx], method = method)
    })
    quantile(vals, c(alpha / 2, 1 - alpha / 2), names = FALSE) |>
      setNames(c("lo", "hi"))
  }
  
  # ---- Compute ----
  rows <- lapply(tr_vec, function(tr) {
    d_tr <- d[treat == tr]
    n    <- nrow(d_tr)
    
    lapply(c("pearson", "spearman"), function(method) {
      r  <- cor(d_tr$lotr_score, d_tr$rt_log_mean, method = method)
      ci <- boot_cor(d_tr$lotr_score, d_tr$rt_log_mean, method = method)
      data.table(
        treatment = tr,
        method    = method,
        n         = n,
        r         = r,
        ci_lo     = ci["lo"],
        ci_hi     = ci["hi"]
      )
    }) |> rbindlist()
  }) |> rbindlist()
  
  # ---- Save ----
  f_out <- file.path(path_out, "ex1_optimism_rt.csv")
  fwrite(rows, f_out)
  msg("Saved: ", f_out)
  
  # ---- Console ----
  msg("\n=== Optimism x log RT correlations ===")
  print(rows[, .(treatment, method, n,
                 r     = round(r,     3),
                 ci_lo = round(ci_lo, 3),
                 ci_hi = round(ci_hi, 3))],
        nrows = Inf)
  
  # ---- Plot ----
  treat_labels <- setNames(
    paste0(tr_vec, " (m = ",
           sapply(tr_vec, function(tr) cfg$design$seq$treatments[[tr]]), ")"),
    tr_vec
  )
  
  d[, treat_label := treat_labels[treat]]
  
  p <- ggplot(d, aes(x = lotr_score, y = rt_log_mean)) +
    geom_point(size = 1.5, alpha = 0.5,
               colour = "#185FA5") +
    geom_smooth(method = "lm", se = TRUE,
                colour = "#185FA5", fill = "#B5D4F4",
                linewidth = 0.7) +
    geom_text(
      data = rows[, .(treat_label = treat_labels[treatment],
                      method, r, ci_lo, ci_hi)],
      aes(label = paste0(tools::toTitleCase(method), "  r = ",
                         round(r, 2),
                         " [", round(ci_lo, 2), ", ",
                         round(ci_hi, 2), "]")),
      x = -Inf, y = Inf, hjust = -0.05, vjust = 1.4 + 1.4 * (method == "spearman"),
      size = 3, colour = "grey40", inherit.aes = FALSE
    ) +
    facet_wrap(~ treat_label, scales = "free") +
    labs(
      x       = "LOT-R score (optimism)",
      y       = "Mean log RT on betting trials",
      caption = paste0("CIs based on ", n_boot, " bootstrap resamples.")
    ) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95", colour = "grey80"),
      plot.caption     = element_text(size = 7, colour = "grey50", hjust = 0)
    )
  

  f_png <- file.path(path_fig, "ex1_optimism_rt.png")
  ggsave(f_pdf, p, width = 3.5 * length(tr_vec), height = 4)
  ggsave(f_png, p, width = 3.5 * length(tr_vec), height = 4, dpi = 300)
  msg("Saved: ", f_pdf)
  msg("Saved: ", f_png)
  
  invisible(list(data = d, correlations = rows, plot = p))
}