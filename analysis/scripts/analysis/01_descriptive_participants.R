# ============================================================
# scripts/analysis/01_descriptive_participants.R
# Descriptive statistics at participant level
# ============================================================

source(here::here("scripts", "00_setup.R"))
library(data.table)

descriptive_participants <- function(cfg) {
  
  ds   <- as.character(cfg$run$dataset)
  seed <- as.integer(cfg$run$seed)
  tr_vec <- unique(as.character(cfg$run$treatment))
  
  # ----------------------------
  # Inputs
  # ----------------------------
  f_par  <- file.path(path_clean_ds(ds), "participants.csv")
  f_lotr <- file.path(path_clean_ds(ds), "lotr_scored.csv")
  f_mpl  <- file.path(path_out_ds(ds), "mpl_scored.csv") # would be better to read from models/ but is too much work to redo. This will do ...
  
  stopifnot(file.exists(f_par), file.exists(f_mpl), file.exists(f_lotr))
  
  par  <- fread(f_par,  encoding = "UTF-8")
  mpl  <- fread(f_mpl,  encoding = "UTF-8")
  lotr <- fread(f_lotr, encoding = "UTF-8")
  
  stopifnot(all(c("pid","age","sex","treat","payoff") %in% names(par)))
  stopifnot(all(c("pid","r_mean","inconsistent") %in% names(mpl)))
  stopifnot(all(c("pid","lotr_score") %in% names(lotr)))
  
  par[, pid := as.character(pid)]
  par[, sex := fifelse(sex %in% c("F","M"), sex, "O")]
  par[is.na(sex) | !nzchar(sex), sex := "O"]
  par[, treat := as.character(treat)]
  
  mpl[, pid := as.character(pid)]
  lotr[, pid := as.character(pid)]
  
  mpl[, inconsistent := as.integer(inconsistent)]
  
  stopifnot(mpl[, uniqueN(pid)] == nrow(mpl))
  stopifnot(lotr[, uniqueN(pid)] == nrow(lotr))
  
  dt <- merge(par[, .(pid, age, sex, treat, payoff)],
              mpl[, .(pid, r_mean, inconsistent)],
              by="pid", all.x=TRUE)
  
  dt <- merge(dt,
              lotr[, .(pid, lotr_score)],
              by="pid", all.x=TRUE)
  
  # ----------------------------
  # Summary helper
  # ----------------------------
  one_row <- function(x) {
    x[, .(
      N = uniqueN(pid),
      
      age_mean = mean(age, na.rm=TRUE),
      age_sd   = sd(age, na.rm=TRUE),
      
      sex_F_n = sum(sex=="F"),
      sex_M_n = sum(sex=="M"),
      sex_O_n = sum(sex=="O"),
      
      payoff_mean = mean(payoff, na.rm=TRUE),
      payoff_sd   = sd(payoff, na.rm=TRUE),
      
      lotr_mean = mean(lotr_score, na.rm=TRUE),
      lotr_sd   = sd(lotr_score, na.rm=TRUE),
      
      r_mean = mean(r_mean, na.rm=TRUE),
      r_sd   = sd(r_mean, na.rm=TRUE),
      
      hl_inconsistent_n    = sum(inconsistent==1L, na.rm=TRUE),
      hl_inconsistent_frac = mean(inconsistent==1L, na.rm=TRUE),
      
      r_mean_F = mean(r_mean[sex=="F"], na.rm=TRUE),
      r_mean_M = mean(r_mean[sex=="M"], na.rm=TRUE),
      r_mean_O = mean(r_mean[sex=="O"], na.rm=TRUE)
    )]
  }
  
  outputs <- list()
  
  # pooled
  out_pooled <- one_row(dt)
  f_pooled <- file.path(path_out_ds(ds), "descriptive_participants_pooled.csv")
  fwrite(out_pooled, f_pooled)
  msg("Saved:", f_pooled)
  outputs$pooled <- out_pooled
  
  # by treatment
  for (tr in cfg$plan$by) {
    x <- dt[treat == tr]
    if (nrow(x) == 0) next
    
    out_tr <- one_row(x)
    f_tr <- file.path(path_out_ds(ds),
                      paste0("descriptive_participants_", tr, ".csv"))
    fwrite(out_tr, f_tr)
    msg("Saved:", f_tr)
    
    outputs[[tr]] <- out_tr
  }
  
  invisible(outputs)
}