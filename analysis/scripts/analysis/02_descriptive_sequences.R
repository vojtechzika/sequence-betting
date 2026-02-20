# scripts/analysis/12_descriptive_sequences.R
source(here::here("scripts", "00_setup.R"))

library(data.table)

descriptive_sequences <- function(dataset = "pilot") {
  
  infile <- file.path(path_clean_ds(dataset), "master_sequences.csv")
  stopifnot(file.exists(infile))
  
  dt <- fread(infile, encoding = "UTF-8")
  stopifnot(all(c("pid", "uid", "seq", "treat", "block", "btn_order", "sex", "stake", "side", "screen_ms") %in% names(dt)))
  
  dt[, sex := as.character(sex)]
  dt[, treat := as.character(treat)]
  dt[, btn_order := as.character(btn_order)]
  
  # ------------------------------------------------------------
  # Definitions
  # ------------------------------------------------------------
  # Bet indicator (robust): stake > 0. If your no-bet is coded differently, change here.
  dt[, bet := !is.na(stake) & stake > 0]
  
  # Side-specific bet indicators
  dt[, bet_H := bet & side == "H"]
  dt[, bet_O := bet & side == "O"]
  
  # ------------------------------------------------------------
  # Helper: per-sequence descriptives (returns one table)
  # ------------------------------------------------------------
  seq_desc <- function(x) {
    
    out <- x[, .(
      n_pid = uniqueN(pid),
      
      bet_rate_mean = mean(bet, na.rm = TRUE),
      bet_rate_sd   = sd(bet, na.rm = TRUE),
      
      bet_rate_H_mean = mean(bet_H, na.rm = TRUE),
      bet_rate_H_sd   = sd(bet_H, na.rm = TRUE),
      
      bet_rate_O_mean = mean(bet_O, na.rm = TRUE),
      bet_rate_O_sd   = sd(bet_O, na.rm = TRUE),
      
      # Stakes: conditional on betting
      stake_mean = mean(stake[bet], na.rm = TRUE),
      stake_sd   = sd(stake[bet], na.rm = TRUE),
      
      # Stakes conditional on side (also conditional on betting)
      stake_H_mean = mean(stake[bet & side == "H"], na.rm = TRUE),
      stake_H_sd   = sd(stake[bet & side == "H"], na.rm = TRUE),
      
      stake_O_mean = mean(stake[bet & side == "O"], na.rm = TRUE),
      stake_O_sd   = sd(stake[bet & side == "O"], na.rm = TRUE),
      
      screen_ms_mean = mean(screen_ms, na.rm = TRUE),
      screen_ms_sd   = sd(screen_ms, na.rm = TRUE)
    ), by = .(seq)]
    
    # Alphabetical order by seq
    setorder(out, seq)
    out
  }
  
  # ------------------------------------------------------------
  # OUTPUT 1: all sequences (all treatments pooled)
  # ------------------------------------------------------------
  out_all <- seq_desc(dt)
  f_all <- file.path(path_out_ds(dataset), paste0("descriptive_sequences_", dataset, ".csv"))
  fwrite(out_all, f_all)
  msg("Saved:", f_all)
  
  # ------------------------------------------------------------
  # OUTPUT 2–3: by treatment (two files)
  # ------------------------------------------------------------
  for (tr in c("m25", "m19")) {
    x <- dt[treat == tr]
    if (nrow(x) == 0) {
      warning("No rows for treat = ", tr, " (dataset=", dataset, "). Skipping file.")
      next
    }
    out_tr <- seq_desc(x)
    f_tr <- file.path(path_out_ds(dataset), paste0("descriptive_sequences_", tr, "_", dataset, ".csv"))
    fwrite(out_tr, f_tr)
    msg("Saved:", f_tr)
  }
  
  # ------------------------------------------------------------
  # OUTPUT 4: grouped table (block, btn_order, sex) with 95% CI
  # IMPORTANT: CI computed from participant-level aggregation (avoids pseudo-replication)
  # Rows: block1..4, orderOH/orderHO, sexM/sexF
  # ------------------------------------------------------------
  
  mean_ci <- function(v) {
    v <- v[is.finite(v)]
    n <- length(v)
    if (n == 0) return(list(mean = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_, n = 0L))
    m <- mean(v)
    se <- sd(v) / sqrt(n)
    list(mean = m, ci_lo = m - 1.96 * se, ci_hi = m + 1.96 * se, n = n)
  }
  
  # Participant-level summary for a given grouping variable
  group_summary <- function(group_var, levels_keep, row_label_fun) {
    
    # aggregate within pid x group_var first
    pid_g <- dt[, .(
      bet_rate      = mean(bet, na.rm = TRUE),
      bet_rate_H    = mean(bet_H, na.rm = TRUE),
      bet_rate_O    = mean(bet_O, na.rm = TRUE),
      
      stake_mean    = mean(stake[bet], na.rm = TRUE),
      stake_H_mean  = mean(stake[bet & side == "H"], na.rm = TRUE),
      stake_O_mean  = mean(stake[bet & side == "O"], na.rm = TRUE),
      
      screen_ms_mean = mean(screen_ms, na.rm = TRUE)
    ), by = c("pid", group_var)]
    
    # keep only requested levels
    pid_g <- pid_g[get(group_var) %in% levels_keep]
    
    # CI across participants for each level
    out <- pid_g[, {
      b  <- mean_ci(bet_rate)
      bH <- mean_ci(bet_rate_H)
      bO <- mean_ci(bet_rate_O)
      
      s  <- mean_ci(stake_mean)
      sH <- mean_ci(stake_H_mean)
      sO <- mean_ci(stake_O_mean)
      
      sc <- mean_ci(screen_ms_mean)
      
      .(
        N_pid = b$n,
        
        bet_rate_mean = b$mean,
        bet_rate_ci_lo = b$ci_lo,
        bet_rate_ci_hi = b$ci_hi,
        
        bet_rate_H_mean = bH$mean,
        bet_rate_H_ci_lo = bH$ci_lo,
        bet_rate_H_ci_hi = bH$ci_hi,
        
        bet_rate_O_mean = bO$mean,
        bet_rate_O_ci_lo = bO$ci_lo,
        bet_rate_O_ci_hi = bO$ci_hi,
        
        stake_mean = s$mean,
        stake_ci_lo = s$ci_lo,
        stake_ci_hi = s$ci_hi,
        
        stake_H_mean = sH$mean,
        stake_H_ci_lo = sH$ci_lo,
        stake_H_ci_hi = sH$ci_hi,
        
        stake_O_mean = sO$mean,
        stake_O_ci_lo = sO$ci_lo,
        stake_O_ci_hi = sO$ci_hi,
        
        screen_ms_mean = sc$mean,
        screen_ms_ci_lo = sc$ci_lo,
        screen_ms_ci_hi = sc$ci_hi
      )
    }, by = group_var]
    
    setnames(out, group_var, "level")
    out[, row := row_label_fun(level)]
    out[, level := NULL]
    setcolorder(out, c("row", setdiff(names(out), "row")))
    out[]
  }
  
  # blocks
  blocks_present <- sort(unique(dt$block))
  block_levels <- intersect(blocks_present, c(1, 2, 3, 4))
  out_block <- group_summary(
    group_var = "block",
    levels_keep = block_levels,
    row_label_fun = function(z) paste0("block", z)
  )
  
  # button order
  # expected values: "OH" / "HO" (keep exactly what exists)
  ord_present <- sort(unique(na.omit(dt$btn_order)))
  ord_keep <- intersect(ord_present, c("OH", "HO"))
  out_order <- group_summary(
    group_var = "btn_order",
    levels_keep = ord_keep,
    row_label_fun = function(z) paste0("order", z)
  )
  
  # sex (requested: M/F only)
  sex_present <- sort(unique(na.omit(dt$sex)))
  sex_keep <- intersect(sex_present, c("M", "F"))
  out_sex <- group_summary(
    group_var = "sex",
    levels_keep = sex_keep,
    row_label_fun = function(z) paste0("sex", z)
  )
  
  out_groups <- rbindlist(list(out_block, out_order, out_sex), fill = TRUE)
  
  f_groups <- file.path(path_out_ds(dataset), paste0("descriptive_sequences_groups_", dataset, ".csv"))
  fwrite(out_groups, f_groups)
  msg("Saved:", f_groups)
  
  invisible(list(
    all = out_all,
    groups = out_groups
  ))
}

# Example:
# descriptive_sequences("pilot")
# descriptive_sequences("main")