# scripts/analysis/11_descriptive_participants.R
source(here::here("scripts", "00_setup.R"))

library(data.table)

descriptive_participants <- function(dataset = "pilot") {
  
  # ----------------------------
  # Inputs
  # ----------------------------
  f_par  <- file.path(path_clean_ds(dataset), "participants.csv")
  f_mpl  <- file.path(path_clean_ds(dataset), "mpl_scored.csv")
  f_lotr <- file.path(path_clean_ds(dataset), "lotr_scored.csv")
  
  stopifnot(file.exists(f_par), file.exists(f_mpl), file.exists(f_lotr))
  
  par  <- fread(f_par,  encoding = "UTF-8")
  mpl  <- fread(f_mpl,  encoding = "UTF-8")
  lotr <- fread(f_lotr, encoding = "UTF-8")
  
  stopifnot(all(c("pid", "age", "sex", "treat", "payoff") %in% names(par)))
  stopifnot(all(c("pid", "r_mean", "inconsistent") %in% names(mpl)))
  stopifnot(all(c("pid", "lotr_score") %in% names(lotr)))
  
  par[, sex := as.character(sex)]
  par[, treat := as.character(treat)]
  
  # Ensure inconsistent is 0/1 (handles TRUE/FALSE or 0/1)
  mpl[, inconsistent := {
    x <- inconsistent
    if (is.logical(x)) as.integer(x) else as.integer(as.character(x))
  }]
  
  # One row per pid in mpl/lotr
  stopifnot(mpl[, uniqueN(pid)] == nrow(mpl))
  stopifnot(lotr[, uniqueN(pid)] == nrow(lotr))
  
  # Merge pid-level info for summaries
  dt <- merge(par[, .(pid, age, sex, treat, payoff)], mpl[, .(pid, r_mean, inconsistent)], by = "pid", all.x = TRUE)
  dt <- merge(dt, lotr[, .(pid, lotr_score)], by = "pid", all.x = TRUE)
  
  # ----------------------------
  # Overall summary (ONE ROW)
  # ----------------------------
   sex_levels <- c("F", "M", "other")
  
  out <- dt[, .(
    N = uniqueN(pid),
    
    sex_F_n     = sum(sex == "F", na.rm = TRUE),
    sex_M_n     = sum(sex == "M", na.rm = TRUE),
    sex_other_n = sum(sex == "other", na.rm = TRUE),
    
    treat_m19_n = sum(treat == "m19", na.rm = TRUE),
    treat_m25_n = sum(treat == "m25", na.rm = TRUE),
    
    payoff_mean = mean(payoff, na.rm = TRUE),
    payoff_sd   = sd(payoff, na.rm = TRUE),
    
    age_mean = mean(age, na.rm = TRUE),
    age_sd   = sd(age, na.rm = TRUE),
    
    lotr_mean = mean(lotr_score, na.rm = TRUE),
    lotr_sd   = sd(lotr_score, na.rm = TRUE),
    
    r_mean = mean(r_mean, na.rm = TRUE),
    r_sd   = sd(r_mean, na.rm = TRUE),
    
    hl_inconsistent_n    = sum(inconsistent == 1L, na.rm = TRUE),
    hl_inconsistent_frac = mean(inconsistent == 1L, na.rm = TRUE)
    
    
  )]
  
  # ----------------------------
  # r_mean by sex (sanity check)
  # ----------------------------
  r_by_sex <- dt[, .(
    r_mean_by_sex = mean(r_mean, na.rm = TRUE)
  ), by = sex]
  
  # Put r_mean_by_sex into wide columns: r_mean_F, r_mean_M, r_mean_other
  r_wide <- dcast(r_by_sex, . ~ sex, value.var = "r_mean_by_sex")
  r_wide[, "." := NULL]
  
  # Ensure columns exist even if a category is absent
  for (s in sex_levels) {
    if (!s %in% names(r_wide)) r_wide[, (s) := NA_real_]
  }
  
  # Rename to requested style
  setnames(r_wide,
           old = sex_levels,
           new = paste0("r_mean_", sex_levels))
  
  # Final single-row output
  out <- cbind(out, r_wide)
  
  # ----------------------------
  # Save
  # ----------------------------
  outfile <- file.path(path_out_ds(dataset), paste0("descriptive_participants_", dataset, ".csv"))
  fwrite(out, outfile)
  
  msg("Participants descriptives saved:", outfile)
  
  invisible(out)
}

# Example:
# participants_descriptives("pilot")
# participants_descriptives("main")