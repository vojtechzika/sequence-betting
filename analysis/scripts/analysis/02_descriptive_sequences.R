# ============================================================
# scripts/analysis/02_descriptive_sequences.R
# Descriptive statistics at sequence level
# ============================================================

source(here::here("scripts", "00_setup.R"))
library(data.table)

descriptive_sequences <- function(cfg) {
  
  ds   <- as.character(cfg$run$dataset)
  seed <- as.integer(cfg$run$seed)
  tr_vec <- unique(as.character(cfg$run$treatment))
  
  infile <- file.path(path_clean_ds(ds), "master_sequences.csv")
  stopifnot(file.exists(infile))
  
  dt <- fread(infile, encoding="UTF-8")
  
  required <- c("pid","seq","treat","block","btn_order",
                "sex","stake","side","screen_ms")
  stopifnot(all(required %in% names(dt)))
  
  dt[, pid   := as.character(pid)]
  dt[, treat := as.character(treat)]
  dt[, seq   := as.character(seq)]
  dt[, side  := as.character(side)]
  dt[, btn_order := as.character(btn_order)]
  
  dt[, sex := fifelse(sex %in% c("F","M"), sex, "O")]
  dt[is.na(sex) | !nzchar(sex), sex := "O"]
  
  dt[, bet   := !is.na(stake) & stake > 0]
  dt[, bet_H := bet & side=="H"]
  dt[, bet_O := bet & side=="O"]
  
  # --------------------------------------------------------
  seq_desc <- function(x) {
    
    out <- x[, .(
      n_pid = uniqueN(pid),
      
      bet_rate_mean = mean(bet, na.rm=TRUE),
      bet_rate_sd   = sd(bet, na.rm=TRUE),
      
      stake_mean = mean(stake[bet], na.rm=TRUE),
      stake_sd   = sd(stake[bet], na.rm=TRUE),
      
      screen_ms_mean = mean(screen_ms, na.rm=TRUE),
      screen_ms_sd   = sd(screen_ms, na.rm=TRUE)
    ), by=seq]
    
    setorder(out, seq)
    out
  }
  
  run_one <- function(x, tag) {
    
    out_seq <- seq_desc(x)
    
    f_seq <- file.path(path_out_ds(ds),
                       paste0("descriptive_sequences_", tag, ".csv"))
    fwrite(out_seq, f_seq)
    msg("Saved:", f_seq)
    
    out_seq
  }
  
  outputs <- list()
  
  outputs$pooled <- run_one(dt, "pooled")
  
  for (tr in cfg$plan$by) {
    x <- dt[treat==tr]
    if (nrow(x)==0) next
    outputs[[tr]] <- run_one(x, tr)
  }
  
  invisible(outputs)
}