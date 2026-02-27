# ------------------------------------------------------------
# Add-on: group-level comparison of r between treatments
#   - reports mean(r_mean) by treatment
#   - reports delta = mean_m25 - mean_m19
#   - bootstrap CI over participants (diagnostic only)
# ------------------------------------------------------------

compare_group_r_by_treat <- function(cfg, tr_a = "m25", tr_b = "m19", B = 5000L) {
  
  stopifnot(is.list(cfg), !is.null(cfg$run), !is.null(cfg$run$dataset))
  ds <- as.character(cfg$run$dataset)
  
  out_dir <- path_out_ds(ds)
  
  f_a <- file.path(path_clean_ds(ds), paste0("mpl_scored_", tr_a, ".csv"))
  f_b <- file.path(path_clean_ds(ds), paste0("mpl_scored_", tr_b, ".csv"))
  stopifnot(file.exists(f_a), file.exists(f_b))
  
  da <- data.table::fread(f_a)
  db <- data.table::fread(f_b)
  
  stopifnot(all(c("pid", "r_mean") %in% names(da)))
  stopifnot(all(c("pid", "r_mean") %in% names(db)))
  
  da <- da[!is.na(r_mean)]
  db <- db[!is.na(r_mean)]
  
  ma <- mean(da$r_mean)
  mb <- mean(db$r_mean)
  delta <- ma - mb
  
  # bootstrap over participants (diagnostic)
  set.seed(as.integer(cfg$run$seed) + 4101L)
  Ba <- nrow(da); Bb <- nrow(db)
  stopifnot(Ba >= 2L, Bb >= 2L)
  
  deltas <- numeric(B)
  for (b in seq_len(B)) {
    sa <- da[sample.int(Ba, Ba, replace = TRUE), mean(r_mean)]
    sb <- db[sample.int(Bb, Bb, replace = TRUE), mean(r_mean)]
    deltas[b] <- sa - sb
  }
  
  res <- data.table::data.table(
    dataset = ds,
    tr_a = tr_a, N_a = nrow(da), mean_r_a = ma,
    tr_b = tr_b, N_b = nrow(db), mean_r_b = mb,
    delta_mean = delta,
    delta_q025 = as.numeric(stats::quantile(deltas, 0.025)),
    delta_q975 = as.numeric(stats::quantile(deltas, 0.975))
  )
  
  f_out <- file.path(out_dir, paste0("mpl_r_group_compare_", tr_a, "_vs_", tr_b, ".csv"))
  if (!should_skip(f_out, cfg, type = "output", label = paste0("MPL group r compare (", ds, ")"))) {
    data.table::fwrite(res, f_out)
    msg("Saved: ", f_out)
  }
  
  res
}

compare_group_r_by_treat(cfg, tr_a = "m25", tr_b = "m19", B = 5000L)