# ============================================================
# 02_stan_r_from_mpl.R
#
# PURPOSE
#   Estimates individual CRRA risk parameters (r_i) from Holt-Laury
#   MPL choices using a Bayesian Stan model. Fits separately for
#   pooled sample, each treatment, and (optionally) consistent-only
#   subsets when cfg$run$consistent_only == TRUE.
#
# INPUT
#   path_src/mpl.csv
#   path_src/master_sequences.csv  -- for treatment assignment
#   path_src/participants.csv      -- for sex variable in summary
#
# OUTPUT
#   path_mod/mpl_fit_<tag>.rds         -- Stan fit objects
#   path_mod/mpl_r_draws_<tag>.rds     -- posterior draws of r_i
#   path_out/mpl_scored_<tag>.csv      -- r_mean, r_median, inconsistent per pid
#   path_out/mpl_scored_summary.csv    -- summary statistics across groups
#
# TAGS
#   pooled, by treatment
#   pooled_consistent, by treatment _consistent (if consistent_only = TRUE)
#
# NOTES
#   - Join at analysis time from mpl_scored_<tr>.csv
#   - Prior parameters from cfg$model$stan$mpl$prior (Bland 2023)
# ============================================================

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

stan_r_from_mpl <- function(cfg) {
  
  seed             <- as.integer(cfg$run$seed)
  consistent_only  <- isTRUE(cfg$run$consistent_only)
  design           <- cfg$design
  model            <- cfg$model
  
  label_safe  <- toupper(trimws(as.character(design$mpl$label_safe)))
  label_risky <- toupper(trimws(as.character(design$mpl$label_risky)))
  
  stopifnot(
    length(label_safe) == 1L, length(label_risky) == 1L,
    nzchar(label_safe), nzchar(label_risky),
    label_safe != label_risky
  )
  
  stan_file <- here::here("stan", "r_from_mpl.stan")
  f_mpl     <- file.path(path_src, "mpl.csv")
  f_master  <- file.path(path_src, "master_sequences.csv")
  
  stopifnot(file.exists(stan_file), file.exists(f_mpl), file.exists(f_master))
  
  mpl    <- fread(f_mpl,    encoding = "UTF-8")
  master <- fread(f_master, encoding = "UTF-8")
  
  stopifnot("pid" %in% names(mpl), "pid" %in% names(master), "treat" %in% names(master))
  
  mpl[, pid := as.character(pid)]
  master[, pid   := as.character(pid)]
  master[, treat := as.character(treat)]
  
  pid_treat <- unique(master[, .(pid, treat)])
  mpl <- merge(mpl, pid_treat, by = "pid", all.x = TRUE)
  
  if (anyNA(mpl$treat)) {
    bad <- mpl[is.na(treat), unique(pid)]
    stop("Some MPL participants have no treatment match in master_sequences.csv. pid(s): ",
         paste(head(bad, 5), collapse = ", "))
  }
  
  K_cfg  <- as.integer(design$mpl$K)
  A_high <- as.numeric(design$mpl$A_high)
  A_low  <- as.numeric(design$mpl$A_low)
  B_high <- as.numeric(design$mpl$B_high)
  B_low  <- as.numeric(design$mpl$B_low)
  
  stopifnot(K_cfg > 0L, is.finite(A_high), is.finite(A_low),
            is.finite(B_high), is.finite(B_low))
  
  p_sched <- (1:K_cfg) / K_cfg
  stopifnot(length(p_sched) == K_cfg, all(is.finite(p_sched)),
            all(p_sched > 0), all(p_sched <= 1))
  
  choice_cols <- grep("^c\\d+$", names(mpl), value = TRUE)
  choice_cols <- choice_cols[order(as.integer(sub("^c", "", choice_cols)))]
  stopifnot(length(choice_cols) == K_cfg)
  
  st              <- model$stan$mpl
  iter_val        <- as.integer(st$iter)
  warmup_val      <- as.integer(st$warmup)
  chains_val      <- as.integer(st$chains)
  adapt_delta_val <- as.numeric(st$adapt_delta)
  treedepth_val   <- as.integer(st$treedepth)
  prior_r_mean    <- as.numeric(st$prior$r_mean)
  prior_r_sd      <- as.numeric(st$prior$r_sd)
  prior_lambda_mean <- as.numeric(st$prior$lambda_mean)
  prior_lambda_sd   <- as.numeric(st$prior$lambda_sd)

  
  stopifnot(iter_val > 0L, warmup_val >= 0L, chains_val > 0L,
            is.finite(adapt_delta_val), adapt_delta_val > 0, adapt_delta_val < 1,
            treedepth_val > 0L, length(prior_r_mean) == 1L, is.finite(prior_r_mean),
            length(prior_r_sd) == 1L, is.finite(prior_r_sd), prior_r_sd > 0,
            length(prior_lambda_mean) == 1L, is.finite(prior_lambda_mean), prior_lambda_mean > 0,
            length(prior_lambda_sd)   == 1L, is.finite(prior_lambda_sd),   prior_lambda_sd   > 0
            )
  
  writeLines(readLines(stan_file, warn = FALSE), stan_file)
  sm <- rstan::stan_model(stan_file)
  
  # ---- Inconsistency flag ----
  hl_inconsistent_flag <- function(long_dt) {
    setorder(long_dt, pid, row)
    long_dt[, {
      yy       <- y
      n_switch <- sum(yy[-1] != yy[-length(yy)])
      reversal <- any(diff(yy) > 0)
      .(inconsistent = as.integer(n_switch > 1 || reversal))
    }, by = pid]
  }
  
  # ---- Melt and code choices ----
  melt_and_code <- function(mpl_subset) {
    long <- melt(
      mpl_subset,
      id.vars       = c("pid", "treat"),
      measure.vars  = choice_cols,
      variable.name = "row",
      value.name    = "choice"
    )
    long[, row    := as.integer(sub("^c", "", row))]
    long[, choice := toupper(trimws(as.character(choice)))]
    
    bad_labels <- setdiff(unique(na.omit(long$choice)), c(label_safe, label_risky))
    if (length(bad_labels) > 0) {
      warning("Unexpected MPL choice label(s): ", paste(bad_labels, collapse = ", "))
    }
    
    long[, y := fifelse(
      choice == label_safe,  1L,
      fifelse(choice == label_risky, 0L, NA_integer_)
    )]
    
    bad_rows <- long[!is.na(choice) & is.na(y)]
    if (nrow(bad_rows) > 0) {
      stop("Could not map some MPL choices. Example values: ",
           paste(head(unique(bad_rows$choice), 5), collapse = ", "))
    }
    long
  }
  
  # ---- Fit one subset ----
  fit_one <- function(mpl_subset, tag) {
    
    fit_file   <- file.path(path_mod, paste0("mpl_fit_",     tag, ".rds"))
    draw_file  <- file.path(path_mod, paste0("mpl_r_draws_", tag, ".rds"))
    scored_csv <- file.path(path_out, paste0("mpl_scored_",  tag, ".csv"))
    
    skip_model  <- should_skip(c(fit_file, draw_file), cfg, "model",
                               paste0("MPL Stan (", tag, ")"))
    skip_scored <- should_skip(scored_csv, cfg, "output",
                               paste0("MPL scored (", tag, ")"))
    
    if (skip_model && skip_scored) return(invisible(NULL))
    
    long       <- melt_and_code(mpl_subset)
    pid_levels <- sort(unique(long$pid))
    long[, uid := match(pid, pid_levels)]
    Tobs       <- nrow(long)
    inc_tbl    <- hl_inconsistent_flag(long[!is.na(y), .(pid, row, y)])
    
    fit  <- NULL
    post <- NULL
    
    if (!skip_model) {
      data_list <- list(
        N            = length(pid_levels),
        T            = Tobs,
        uid          = long$uid,
        p            = p_sched[long$row],
        A1           = rep(A_high, Tobs),
        A2           = rep(A_low,  Tobs),
        B1           = rep(B_high, Tobs),
        B2           = rep(B_low,  Tobs),
        y            = long$y,
        prior_r_mean = prior_r_mean,
        prior_r_sd   = prior_r_sd,
        prior_lambda_mean = as.numeric(st$prior$lambda_mean),
        prior_lambda_sd   = as.numeric(st$prior$lambda_sd)
      )
      
      msg("MPL Stan: fitting tag=", tag,
          " | N=", data_list$N, " | T=", data_list$T,
          " | prior r ~ Normal(", prior_r_mean, ", ", prior_r_sd, ")")
      
      fit <- rstan::sampling(
        sm,
        data    = data_list,
        iter    = iter_val,
        warmup  = warmup_val,
        chains  = chains_val,
        seed    = seed,
        control = list(adapt_delta = adapt_delta_val, max_treedepth = treedepth_val)
      )
      
      saveRDS(fit, fit_file)
      msg("Saved: ", fit_file)
      
      post <- rstan::extract(fit)
      saveRDS(list(pid = pid_levels, r_draws = post$r), draw_file)
      msg("Saved: ", draw_file)
      
    } else {
      stopifnot(file.exists(fit_file))
      fit  <- readRDS(fit_file)
      post <- rstan::extract(fit)
      stopifnot(!is.null(post$r))
    }
    
    if (!skip_scored) {
      scored <- data.table(
        pid      = pid_levels,
        tag      = tag,
        r_mean   = apply(post$r, 2, mean),
        r_median = apply(post$r, 2, median)
      )
      scored <- merge(scored, inc_tbl, by = "pid", all.x = TRUE)
      scored[is.na(inconsistent), inconsistent := 0L]
      fwrite(scored, scored_csv)
      msg("Saved: ", scored_csv)
    }
    
    invisible(list(fit_file = fit_file, draw_file = draw_file, scored_csv = scored_csv))
  }
  
  # ---- Compute inconsistency flags on full data ----
  long_all        <- melt_and_code(mpl)
  inc_all         <- hl_inconsistent_flag(long_all[!is.na(y), .(pid, row, y)])
  consistent_pids <- inc_all[inconsistent == 0L, pid]
  
  msg("Consistent participants: ", length(consistent_pids),
      " of ", nrow(inc_all), " total")
  
  # ---- Fits ----
  # Pooled
  fit_one(mpl, "pooled")
  
  # Per treatment
  for (tr in names(design$seq$treatments)) {
    mpl_tr <- mpl[treat == tr]
    if (nrow(mpl_tr) == 0) {
      warning("No MPL rows for treat='", tr, "'. Skipping.")
      next
    }
    fit_one(mpl_tr, tr)
  }
  
  # Consistent-only (if requested)
  if (consistent_only) {
    msg("Fitting consistent-only subsets (cfg$run$consistent_only = TRUE)")
    
    fit_one(mpl[pid %in% consistent_pids], "pooled_consistent")
    
    for (tr in names(design$seq$treatments)) {
      mpl_tr_consistent <- mpl[pid %in% consistent_pids & treat == tr]
      if (nrow(mpl_tr_consistent) == 0) {
        warning("No consistent MPL rows for treat='", tr, "'. Skipping.")
        next
      }
      fit_one(mpl_tr_consistent, paste0(tr, "_consistent"))
    }
  }
  
  # ---- Summary CSV ----
  f_participants <- file.path(path_src, "participants.csv")
  participants   <- if (file.exists(f_participants)) {
    p <- fread(f_participants, encoding = "UTF-8")
    p[, pid := as.character(pid)]
    p[, sex_label := fcase(
      toupper(trimws(sex)) == "F", "Women",
      toupper(trimws(sex)) == "M", "Men",
      default = NA_character_
    )]
    p[, .(pid, sex_label)]
  } else NULL
  
  attach_sex <- function(d) {
    if (!is.null(participants))
      merge(d, participants, by = "pid", all.x = TRUE)
    else
      d[, sex_label := NA_character_][]
  }
  
  summarise_scored <- function(d, label) {
    d <- d[!is.na(r_mean)]
    if (nrow(d) == 0) return(NULL)
    data.table(
      group            = label,
      n                = nrow(d),
      r_mean           = mean(d$r_mean),
      r_sd             = sd(d$r_mean),
      r_median         = median(d$r_mean),
      r_q25            = quantile(d$r_mean, 0.25),
      r_q75            = quantile(d$r_mean, 0.75),
      r_min            = min(d$r_mean),
      r_max            = max(d$r_mean),
      pct_inconsistent = mean(d$inconsistent == 1L, na.rm = TRUE)
    )
  }
  
  summary_list <- list()
  
  all_tags <- c("pooled", names(design$seq$treatments))
  if (consistent_only) {
    all_tags <- c(all_tags, "pooled_consistent",
                  paste0(names(design$seq$treatments), "_consistent"))
  }
  
  for (tag in all_tags) {
    f_tag <- file.path(path_out, paste0("mpl_scored_", tag, ".csv"))
    if (!file.exists(f_tag)) next
    
    d <- attach_sex(fread(f_tag))
    
    for (subset_label in c("all", "consistent", "inconsistent")) {
      d_sub <- switch(subset_label,
                      all          = d,
                      consistent   = d[inconsistent == 0L],
                      inconsistent = d[inconsistent == 1L]
      )
      grp <- paste0(tag, " / ", subset_label)
      summary_list[[grp]] <- summarise_scored(d_sub, grp)
      
      for (sx in c("Women", "Men")) {
        grp_sx <- paste0(tag, " / ", subset_label, " / ", sx)
        summary_list[[grp_sx]] <- summarise_scored(d_sub[sex_label == sx], grp_sx)
      }
    }
  }
  
  summary_dt <- rbindlist(summary_list, fill = TRUE)
  f_summary  <- file.path(path_out, "mpl_scored_summary.csv")
  fwrite(summary_dt, f_summary)
  msg("Saved: ", f_summary)
  
  invisible(TRUE)
}