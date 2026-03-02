# scripts/analysis/02_score_lotr.R
source(here::here("scripts", "00_setup.R"))

library(data.table)
library(psych)

score_lotr <- function(cfg) {
  
  dataset <- as.character(cfg$run$dataset)
  
  # =================================================
  # DESIGN PARAMETERS (from design cfg)
  # =================================================
  lotr_cfg <- cfg$design$lotr
  
  LOTR_SCALE_MIN <- as.integer(lotr_cfg$scale_min)
  LOTR_SCALE_MAX <- as.integer(lotr_cfg$scale_max)
  
  item_cols    <- as.character(lotr_cfg$item_cols)
  scored_items <- as.character(lotr_cfg$scored_items)
  rev_items    <- as.character(lotr_cfg$rev_items)
  
  stopifnot(
    length(LOTR_SCALE_MIN) == 1L,
    length(LOTR_SCALE_MAX) == 1L,
    LOTR_SCALE_MAX > LOTR_SCALE_MIN,
    length(item_cols) > 0L,
    length(scored_items) > 0L,
    all(rev_items %in% scored_items),
    all(scored_items %in% item_cols)
  )
  
  # =================================================
  # INPUT
  # =================================================
  infile <- file.path(path_clean_ds(dataset), "lotr.csv")
  stopifnot(file.exists(infile))
  
  dt <- fread(infile, encoding = "UTF-8")
  stopifnot(all(c("pid", item_cols) %in% names(dt)))
  
  # =================================================
  # OUTPUT PATHS
  # =================================================
  out_dir <- path_out_ds(dataset)
  mod_dir <- path_mod_ds(dataset)
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(mod_dir, showWarnings = FALSE, recursive = TRUE)
  
  outfile_scored    <- file.path(path_clean_ds(dataset), "lotr_scored.csv")
  outfile_alpha_rds <- file.path(mod_dir, "lotr_alpha.rds")
  outfile_alpha_csv <- file.path(out_dir, "lotr_alpha.csv")
  
  # =================================================
  # should_skip guards (scored output)
  # =================================================
  if (should_skip(
    paths = outfile_scored,
    cfg   = cfg,
    type  = "output",
    label = paste0("LOT-R scored (", dataset, ")")
  )) {
    return(invisible(NULL))
  }
  
  # =================================================
  # CLEANING
  # =================================================
  dt[, pid := as.character(pid)]
  dt[, (item_cols) := lapply(.SD, as.numeric), .SDcols = item_cols]
  
  # Guard: no missing in scored items (forced responses)
  if (!dt[, all(complete.cases(.SD)), .SDcols = scored_items]) {
    bad <- dt[!complete.cases(dt[, ..scored_items]), .(pid)]
    stop("LOT-R: missing values in scored items for pid(s): ", paste(bad$pid, collapse = ", "))
  }
  
  # =================================================
  # SCORING
  # =================================================
  dt[, (rev_items) := lapply(.SD, function(x) LOTR_SCALE_MAX + LOTR_SCALE_MIN - x),
     .SDcols = rev_items]
  
  dt[, lotr_score := rowSums(.SD),  .SDcols = scored_items]
  dt[, lotr_mean  := rowMeans(.SD), .SDcols = scored_items]
  dt[, lotr_z     := as.numeric(scale(lotr_score))]
  
  out_scored <- dt[, .(pid, lotr_score, lotr_mean, lotr_z)]
  fwrite(out_scored, outfile_scored)
  msg("Saved: ", outfile_scored, " | rows: ", nrow(out_scored))
  
  # =================================================
  # Append lotr_score into master_sequences.csv (ONLY lotr_score)
  # =================================================
  f_master <- file.path(path_clean_ds(dataset), "master_sequences.csv")
  stopifnot(file.exists(f_master))
  
  # Guard overwrite behavior for master update
  if (!should_skip(
    paths = f_master,
    cfg   = cfg,
    type  = "output",
    label = paste0("Append lotr_score to master (", dataset, ")")
  )) {
    
    master <- fread(f_master)
    stopifnot("pid" %in% names(master))
    
    master[, pid := as.character(pid)]
    out_scored[, pid := as.character(pid)]
    
    # Remove existing lotr_score if re-running
    if ("lotr_score" %in% names(master)) {
      master[, lotr_score := NULL]
    }
    
    master <- merge(
      master,
      out_scored[, .(pid, lotr_score)],
      by = "pid",
      all.x = TRUE
    )
    
    fwrite(master, f_master)
    msg("Saved: ", f_master, " (updated with lotr_score)")
  }
  
  # =================================================
  # RELIABILITY: Cronbach alpha
  # =================================================
  alpha_result <- psych::alpha(dt[, ..scored_items])
  alpha_val    <- alpha_result$total$raw_alpha
  
  msg("Cronbach alpha (LOT-R scored items): ", round(alpha_val, 3))
  
  alpha_interpretation <- if (alpha_val >= 0.9) {
    "Excellent internal consistency (may be inflated in small samples)."
  } else if (alpha_val >= 0.8) {
    "Good internal consistency."
  } else if (alpha_val >= 0.7) {
    "Acceptable internal consistency."
  } else if (alpha_val >= 0.6) {
    "Questionable internal consistency."
  } else if (alpha_val >= 0.5) {
    "Poor internal consistency."
  } else {
    "Unacceptable internal consistency."
  }
  
  # Save RDS (model artifact)
  if (!should_skip(
    paths = outfile_alpha_rds,
    cfg   = cfg,
    type  = "model",
    label = paste0("LOT-R alpha rds (", dataset, ")")
  )) {
    saveRDS(alpha_result, outfile_alpha_rds)
    msg("Saved: ", outfile_alpha_rds)
  }
  
  alpha_table <- data.table(
    dataset = dataset,
    n_participants = nrow(dt),
    n_items = length(scored_items),
    raw_alpha = round(alpha_val, 4),
    standardized_alpha = round(alpha_result$total$std.alpha, 4),
    average_inter_item_correlation = round(alpha_result$total$average_r, 4),
    interpretation = alpha_interpretation
  )
  
  # Save CSV (output artifact)
  if (!should_skip(
    paths = outfile_alpha_csv,
    cfg   = cfg,
    type  = "output",
    label = paste0("LOT-R alpha csv (", dataset, ")")
  )) {
    fwrite(alpha_table, outfile_alpha_csv)
    msg("Saved: ", outfile_alpha_csv)
  }
  
  invisible(list(scored = out_scored, alpha = alpha_result, alpha_table = alpha_table))
}