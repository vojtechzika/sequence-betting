# scripts/analysis/02_descriptive_sequences.R
source(here::here("scripts", "00_setup.R"))

library(data.table)

descriptive_sequences <- function(cfg) {
  
  stopifnot(is.list(cfg), !is.null(cfg$run), !is.null(cfg$run$dataset))
  stopifnot(!is.null(cfg$plan), !is.null(cfg$plan$by))
  
  ds <- cfg$run$dataset
  
  infile <- file.path(path_clean_ds(ds), "master_sequences.csv")
  stopifnot(file.exists(infile))
  
  dt <- fread(infile, encoding = "UTF-8")
  stopifnot(all(c("pid", "uid", "seq", "treat", "block", "btn_order", "sex", "stake", "side", "screen_ms") %in% names(dt)))
  
  dt[, pid       := as.character(pid)]
  dt[, sex       := as.character(sex)]
  dt[, treat     := as.character(treat)]
  dt[, btn_order := as.character(btn_order)]
  dt[, side      := as.character(side)]
  dt[, seq       := as.character(seq)]
  
  # oTree convention: anything not exactly F/M -> O
  dt[, sex := fifelse(sex %in% c("F", "M"), sex, "O")]
  dt[is.na(sex) | !nzchar(sex), sex := "O"]
  
  # Bet indicator: stake > 0
  dt[, bet := !is.na(stake) & stake > 0]
  dt[, bet_H := bet & side == "H"]
  dt[, bet_O := bet & side == "O"]
  
  # ------------------------------------------------------------
  # Helpers
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
      
      stake_mean = mean(stake[bet], na.rm = TRUE),
      stake_sd   = sd(stake[bet], na.rm = TRUE),
      
      stake_H_mean = mean(stake[bet & side == "H"], na.rm = TRUE),
      stake_H_sd   = sd(stake[bet & side == "H"], na.rm = TRUE),
      
      stake_O_mean = mean(stake[bet & side == "O"], na.rm = TRUE),
      stake_O_sd   = sd(stake[bet & side == "O"], na.rm = TRUE),
      
      screen_ms_mean = mean(screen_ms, na.rm = TRUE),
      screen_ms_sd   = sd(screen_ms, na.rm = TRUE)
    ), by = .(seq)]
    
    setorder(out, seq)
    out
  }
  
  mean_ci <- function(v) {
    v <- v[is.finite(v)]
    n <- length(v)
    if (n == 0) return(list(mean = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_, n = 0L))
    m <- mean(v)
    se <- sd(v) / sqrt(n)
    list(mean = m, ci_lo = m - 1.96 * se, ci_hi = m + 1.96 * se, n = n)
  }
  
  group_summary <- function(x, group_var, levels_keep, row_label_fun) {
    
    pid_g <- x[, .(
      bet_rate       = mean(bet, na.rm = TRUE),
      bet_rate_H     = mean(bet_H, na.rm = TRUE),
      bet_rate_O     = mean(bet_O, na.rm = TRUE),
      
      stake_mean     = mean(stake[bet], na.rm = TRUE),
      stake_H_mean   = mean(stake[bet & side == "H"], na.rm = TRUE),
      stake_O_mean   = mean(stake[bet & side == "O"], na.rm = TRUE),
      
      screen_ms_mean = mean(screen_ms, na.rm = TRUE)
    ), by = c("pid", group_var)]
    
    pid_g <- pid_g[get(group_var) %in% levels_keep]
    
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
        
        bet_rate_mean   = b$mean,
        bet_rate_ci_lo  = b$ci_lo,
        bet_rate_ci_hi  = b$ci_hi,
        
        bet_rate_H_mean  = bH$mean,
        bet_rate_H_ci_lo = bH$ci_lo,
        bet_rate_H_ci_hi = bH$ci_hi,
        
        bet_rate_O_mean  = bO$mean,
        bet_rate_O_ci_lo = bO$ci_lo,
        bet_rate_O_ci_hi = bO$ci_hi,
        
        stake_mean   = s$mean,
        stake_ci_lo  = s$ci_lo,
        stake_ci_hi  = s$ci_hi,
        
        stake_H_mean  = sH$mean,
        stake_H_ci_lo = sH$ci_lo,
        stake_H_ci_hi = sH$ci_hi,
        
        stake_O_mean  = sO$mean,
        stake_O_ci_lo = sO$ci_lo,
        stake_O_ci_hi = sO$ci_hi,
        
        screen_ms_mean  = sc$mean,
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
  
  run_one <- function(x, tag) {
    out_seq <- seq_desc(x)
    f_seq <- file.path(path_out_ds(ds), paste0("descriptive_sequences_", tag, ".csv"))
    fwrite(out_seq, f_seq)
    msg("Saved:", f_seq)
    
    blocks_present <- sort(unique(x$block))
    block_levels <- intersect(blocks_present, c(1, 2, 3, 4))
    out_block <- group_summary(x, "block", block_levels, function(z) paste0("block", z))
    
    ord_present <- sort(unique(na.omit(x$btn_order)))
    ord_keep <- intersect(ord_present, c("OH", "HO"))
    out_order <- group_summary(x, "btn_order", ord_keep, function(z) paste0("order", z))
    
    sex_present <- sort(unique(na.omit(x$sex)))
    sex_keep <- intersect(sex_present, c("F", "M", "O"))
    out_sex <- group_summary(x, "sex", sex_keep, function(z) paste0("sex", z))
    
    out_groups <- rbindlist(list(out_block, out_order, out_sex), fill = TRUE)
    f_groups <- file.path(path_out_ds(ds), paste0("descriptive_sequences_groups_", tag, ".csv"))
    fwrite(out_groups, f_groups)
    msg("Saved:", f_groups)
    
    list(all = out_seq, groups = out_groups)
  }
  
  # ------------------------------------------------------------
  # Execute: always pooled + always by-treatment
  # ------------------------------------------------------------
  outputs <- list()
  
  outputs$pooled <- run_one(dt, tag = "pooled")
  
  for (tr in cfg$plan$by) {
    x <- dt[treat == tr]
    if (nrow(x) == 0) {
      warning("No rows for treat='", tr, "' in dataset='", ds, "'. Skipping.")
      next
    }
    outputs[[tr]] <- run_one(x, tag = tr)
  }
  
  invisible(outputs)
}