# ============================================================
# scripts/analysis/54_ex1_2_stan_ghi.R
#
# EX1.2 Predictors of the GHI (per treatment run)
#
# Uses:
#   chi_i draws from EX1.1 participants RDS:
#       data/clean/<ds>/models/ex1_1_participants_<tr>.rds
#
# Covariates:
#   lotr_score  (from master_sequences.csv)
#   screen_ms   (reaction time -> log RT on betting trials)
#   r_median    (computed from mpl_r_draws_<tr>.rds)
#
# Uncertainty propagation:
#   integrates over chi-draws via MI likelihood in Stan
#
# Outputs:
#   data/clean/<ds>/models/ex1_2_fit_<tr>.rds
#   data/clean/<ds>/output/ex1_2_coeffs_<tr>.csv
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(rstan)
})

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

ex1_2_stan_ghi <- function(cfg) {
  
  stopifnot(is.list(cfg), !is.null(cfg$run), !is.null(cfg$design), !is.null(cfg$model))
  stopifnot(!is.null(cfg$run$dataset), !is.null(cfg$run$seed))
  
  ds   <- as.character(cfg$run$dataset)
  seed <- as.integer(cfg$run$seed)
  
  # ------------------------------------------------------------
  # Run controls
  # ------------------------------------------------------------
  
  Trep_cfg <- cfg$run$ex1_trep
  if (is.null(Trep_cfg)) Trep_cfg <- 300L
  Trep_cfg <- as.integer(Trep_cfg)
  stopifnot(length(Trep_cfg) == 1L, Trep_cfg >= 10L)
  
  rt_scope <- "bet"
  
  # ------------------------------------------------------------
  # Files
  # ------------------------------------------------------------
  
  mod_dir <- path_mod_ds(ds)
  out_dir <- path_out_ds(ds)
  
  dir.create(mod_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  f_master <- file.path(path_clean_ds(ds), "master_sequences.csv")
  stopifnot(file.exists(f_master))
  
  stan_file <- here::here("stan", "ex1_ghi.stan")
  stopifnot(file.exists(stan_file))
  writeLines(readLines(stan_file, warn = FALSE), stan_file)
  
  sm <- rstan::stan_model(stan_file)
  
  # ------------------------------------------------------------
  # Load master
  # ------------------------------------------------------------
  
  master <- fread(f_master, encoding = "UTF-8")
  
  req <- c("pid","treat","stake","screen_ms","lotr_score")
  miss <- setdiff(req, names(master))
  if (length(miss) > 0)
    stop("master_sequences.csv missing: ", paste(miss, collapse=", "))
  
  master[, pid        := as.character(pid)]
  master[, treat      := as.character(treat)]
  master[, stake      := as.numeric(stake)]
  master[, screen_ms  := as.numeric(screen_ms)]
  master[, lotr_score := as.numeric(lotr_score)]
  
  master[is.na(stake), stake := 0]
  
  # ------------------------------------------------------------
  # Treatments
  # ------------------------------------------------------------
  
  tr_vec <- cfg$run$treatment
  
  if (is.null(tr_vec) || length(tr_vec)==0 || all(!nzchar(as.character(tr_vec)))) {
    tr_vec <- sort(unique(master$treat))
  } else {
    tr_vec <- unique(as.character(tr_vec))
  }
  
  stopifnot(length(tr_vec)>0, all(nzchar(tr_vec)))
  
  # ------------------------------------------------------------
  # helper
  # ------------------------------------------------------------
  
  z_within <- function(x){
    m <- mean(x, na.rm=TRUE)
    s <- sd(x, na.rm=TRUE)
    if (!is.finite(s) || s<=0) stop("Cannot z-score: sd <= 0")
    (x-m)/s
  }
  
  outputs <- list()
  
  for (tr in tr_vec) {
    
    # ------------------------------------------------------------
    # output skip logic
    # ------------------------------------------------------------
    
    f_fit <- file.path(mod_dir, paste0("ex1_2_fit_", tr, ".rds"))
    f_csv <- file.path(out_dir, paste0("ex1_2_coeffs_", tr, ".csv"))
    
    skip_model <- should_skip(
      paths = f_fit,
      cfg   = cfg,
      type  = "model",
      label = paste0("EX1.2 GHI (", ds, "/", tr, ")")
    )
    
    skip_output <- should_skip(
      paths = f_csv,
      cfg   = cfg,
      type  = "output",
      label = paste0("EX1.2 coeffs (", ds, "/", tr, ")")
    )
    
    if (skip_model && skip_output) {
      next
    }
    
    # ------------------------------------------------------------
    # load chi draws
    # ------------------------------------------------------------
    
    f_chi <- file.path(mod_dir,paste0("ex1_1_",tr,"_participants.rds"))
    if (!file.exists(f_chi))
      stop("Missing EX1.1 participants file: ",f_chi)
    
    chi_obj <- readRDS(f_chi)
    
    pid_levels <- chi_obj$pid_levels
    if (is.null(pid_levels)) pid_levels <- chi_obj$pid
    
    chi_draws <- chi_obj$chi_draws
    if (is.null(chi_draws)) chi_draws <- chi_obj$chi_rep
    
    stopifnot(is.matrix(chi_draws))
    stopifnot(ncol(chi_draws)==length(pid_levels))
    
    K_all <- nrow(chi_draws)
    
    Trep <- min(Trep_cfg,K_all)
    
    set.seed(seed)
    t_idx <- sort(sample.int(K_all,Trep,FALSE))
    y_rep <- chi_draws[t_idx,,drop=FALSE]
    
    # ------------------------------------------------------------
    # load MPL r draws
    # ------------------------------------------------------------
    
    r_file <- file.path(mod_dir,paste0("mpl_r_draws_",tr,".rds"))
    if (!file.exists(r_file))
      stop("Missing MPL r file: ",r_file)
    
    r_obj <- readRDS(r_file)
    
    r_dt <- data.table(
      pid = as.character(r_obj$pid),
      r_mean = apply(r_obj$r_draws,2,mean)
    )
    
    # ------------------------------------------------------------
    # build covariates
    # ------------------------------------------------------------
    
    d <- master[treat==tr & pid %in% pid_levels]
    
    d_rt <- d[stake>0 & is.finite(screen_ms) & screen_ms>0]
    
    rt_pid <- d_rt[,.(rt_log_mean = mean(log(screen_ms))),by=pid]
    
    pid_cov <- d[,.(lotr_score = lotr_score[1]),by=pid]
    
    pid_cov <- merge(pid_cov,r_dt,by="pid",all.x=TRUE)
    pid_cov <- merge(pid_cov,rt_pid,by="pid",all.x=TRUE)
    
    pid_cov <- pid_cov[pid %in% pid_levels]
    
    setkey(pid_cov,pid)
    pid_cov <- pid_cov[.(pid_levels)]
    
    keep <- is.finite(pid_cov$lotr_score) &
      is.finite(pid_cov$r_mean) &
      is.finite(pid_cov$rt_log_mean)
    
    pid_keep <- pid_cov[keep,pid]
    
    pid_cov <- pid_cov[keep]
    y_rep   <- y_rep[,keep,drop=FALSE]
    
    N <- nrow(pid_cov)
    
    # ------------------------------------------------------------
    # standardize predictors
    # ------------------------------------------------------------
    
    Zopt <- z_within(pid_cov$lotr_score)
    Zrt  <- z_within(pid_cov$rt_log_mean)
    Zr   <- z_within(pid_cov$r_mean)
    
    # ------------------------------------------------------------
    # stan
    # ------------------------------------------------------------
    
    data_list <- list(
      N=N,
      Trep=Trep,
      y_rep=y_rep,
      Zopt=as.vector(Zopt),
      Zrt=as.vector(Zrt),
      Zr=as.vector(Zr)
    )
    
    iter_val   <- cfg$model$stan$ex1_2[[ds]]$iter
    warmup_val <- cfg$model$stan$ex1_2[[ds]]$warmup
    chains_val <- cfg$model$stan$ex1_2[[ds]]$chains
    adapt_delta_val <- cfg$model$stan$ex1_2[[ds]]$adapt_delta
    treedepth_val   <- cfg$model$stan$ex1_2[[ds]]$treedepth
    
    if (!skip_model) {
      fit <- rstan::sampling(
        sm,
        data=data_list,
        iter=iter_val,
        warmup=warmup_val,
        chains=chains_val,
        seed=seed,
        control=list(adapt_delta=adapt_delta_val,
                     max_treedepth=treedepth_val)
      )
      
      # ------------------------------------------------------------
      # save fit
      # ------------------------------------------------------------
      
      saveRDS(fit,f_fit)
    } else {
      fit <- readRDS(f_fit)
    }
    
    # ------------------------------------------------------------
    # summarize
    # ------------------------------------------------------------
    
    post <- rstan::extract(fit)
    
    summ <- function(x){
      q <- quantile(x, c(0.025, 0.975), names = FALSE)
      out <- c(
        median = median(x),
        mean   = mean(x),
        q025   = q[1],
        q975   = q[2],
        p_gt0  = mean(x > 0)
      )
      return(out)
    }
    
    coefs <- rbind(
      alpha    = summ(post$alpha),
      beta_opt = summ(post$b_opt),
      beta_rt  = summ(post$b_rt),
      beta_r   = summ(post$b_r)
    )
    
    tbl <- data.table(
      term   = rownames(coefs),
      median = coefs[,"median"],
      mean   = coefs[,"mean"],
      q025   = coefs[,"q025"],
      q975   = coefs[,"q975"],
      p_gt0  = coefs[,"p_gt0"]
    )
    
    tbl[,`:=`(
      dataset=ds,
      treatment=tr,
      N=N,
      Trep=Trep,
      rt_scope=rt_scope
    )]
    
    setcolorder(tbl,
                c("dataset","treatment","N","Trep","rt_scope",
                  "term","median","mean","q025","q975","p_gt0"))
    
    if (!skip_output) {
      fwrite(tbl,f_csv)
    }
    
    outputs[[tr]] <- list(
      fit_file=f_fit,
      coeffs_csv=f_csv,
      pid_levels=pid_keep
    )
  }
  
  invisible(outputs)
}

# Example
# ex1_2_stan_ghi(cfg)