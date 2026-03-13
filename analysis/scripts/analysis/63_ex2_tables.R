# ============================================================
# scripts/analysis/61_ex2_stan.R
# ============================================================

library(data.table)
library(rstan)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

ex2_stan <- function(cfg) {
  
  stopifnot(is.list(cfg), !is.null(cfg$run), !is.null(cfg$model))
  stopifnot(!is.null(cfg$run$dataset), !is.null(cfg$run$seed))
  
  ds   <- as.character(cfg$run$dataset)
  seed <- as.integer(cfg$run$seed)
  
  tr_vec <- as.character(cfg$run$treatment)
  
  mod_dir <- path_mod_ds(ds)
  out_dir <- path_out_ds(ds)
  
  dir.create(mod_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  f_master <- file.path(path_clean_ds(ds), "master_sequences.csv")
  stopifnot(file.exists(f_master))
  
  stan_file <- here::here("stan", "ex2_associations.stan")
  stopifnot(file.exists(stan_file))
  
  writeLines(readLines(stan_file, warn = FALSE), stan_file)
  
  sm <- rstan::stan_model(stan_file)
  
  st <- cfg$model$stan$ex2[[ds]]
  
  iter_val        <- as.integer(st$iter)
  warmup_val      <- as.integer(st$warmup)
  chains_val      <- as.integer(st$chains)
  adapt_delta_val <- as.numeric(st$adapt_delta)
  treedepth_val   <- as.integer(st$treedepth)
  
  ex2_trep <- as.integer(cfg$model$simulation$ex2_trep[[ds]])
  
  min_Nk <- if (ds == "pilot") 0L else as.integer(cfg$design$ex2$min_Nk)
  
  rq2_min_bets <- as.integer(cfg$design$rq2$min_bets)
  
  master <- fread(f_master, encoding = "UTF-8")
  
  req <- c("pid","treat","stake","screen_ms","lotr_score")
  miss <- setdiff(req,names(master))
  if(length(miss)>0)
    stop("master_sequences.csv missing columns: ",paste(miss,collapse=", "))
  
  master[, pid := as.character(pid)]
  master[, treat := as.character(treat)]
  master[, stake := as.numeric(stake)]
  master[, screen_ms := as.numeric(screen_ms)]
  master[, lotr_score := as.numeric(lotr_score)]
  
  master[is.na(stake), stake := 0]
  
  master[, lotr_z := (lotr_score - mean(lotr_score, na.rm = TRUE)) /
           sd(lotr_score, na.rm = TRUE)]
  
  z_strict <- function(x){
    m <- mean(x)
    s <- sd(x)
    if(!is.finite(s) || s<=0) stop("Cannot z-score.")
    (x-m)/s
  }
  
  z_within_draw <- function(mat){
    for(t in seq_len(nrow(mat))){
      m <- mean(mat[t,])
      s <- sd(mat[t,])
      if(!is.finite(s) || s<=0) stop("draw sd zero")
      mat[t,] <- (mat[t,]-m)/s
    }
    mat
  }
  
  outputs <- list()
  
  for(tr in tr_vec){
    
    msg("Running EX2 Stan model for treatment: ", tr)
    
    # ---------------------------------------------------------
    # load RQ fits
    # ---------------------------------------------------------
    
    fit1 <- readRDS(file.path(mod_dir,paste0("rq1_fit_sequences_",tr,"_full.rds")))
    fit2 <- readRDS(file.path(mod_dir,paste0("rq2_fit_sequences_",tr,"_full.rds")))
    fit3 <- readRDS(file.path(mod_dir,paste0("rq3_fit_sequences_",tr,"_full.rds")))
    fit4 <- readRDS(file.path(mod_dir,paste0("rq4_fit_sequences_",tr,"_full.rds")))
    
    pid1 <- readRDS(file.path(mod_dir,paste0("rq1_pid_levels_",tr,"_full.rds")))
    pid2 <- readRDS(file.path(mod_dir,paste0("rq2_pid_levels_",tr,"_full.rds")))
    pid3 <- readRDS(file.path(mod_dir,paste0("rq3_pid_levels_",tr,"_full.rds")))
    pid4 <- readRDS(file.path(mod_dir,paste0("rq4_pid_levels_",tr,"_full.rds")))
    
    post1 <- rstan::extract(fit1)
    post2 <- rstan::extract(fit2)
    post3 <- rstan::extract(fit3)
    post4 <- rstan::extract(fit4)
    
    mu_b <- post1$mu_b_i
    mu_a <- post2$mu_a_i
    mu_c <- post3$mu_c_i
    mu_h <- post4$mu_h_i
    
    dtr <- master[treat==tr]
    
    bet_n <- dtr[stake>0,.N,by=pid]
    setkey(bet_n,pid)
    
    build_predictors_all <- function() {
      
      pid_cov <- dtr[, .(lotr_z = lotr_z[1]), by = pid]
      
      d_rt <- dtr[stake > 0 & screen_ms > 0]
      rt_pid <- d_rt[, .(rt_log_mean = mean(log(screen_ms))), by = pid]
      
      pid_cov <- merge(pid_cov, rt_pid, by = "pid", all.x = TRUE)
      
      pid_cov
    }
    
    pid_cov_all <- build_predictors_all()
    
    # participants with valid EX2 predictors
    pid_pred <- pid_cov_all[
      is.finite(lotr_z) & is.finite(rt_log_mean),
      pid
    ]
    
    # common EX2 sample across all outcomes except RQ2's extra min-bets rule
    pid_common <- Reduce(intersect, list(
      pid_pred,
      as.character(pid1),
      as.character(pid3),
      as.character(pid4)
    ))
    
    # RQ2-specific common sample additionally requiring enough bets
    n_bets_all <- bet_n[.(pid_common), N]
    n_bets_all[is.na(n_bets_all)] <- 0
    pid_common_a <- pid_common[n_bets_all >= rq2_min_bets]
    
    # z-score once on the common treatment-level EX2 sample
    pid_cov_common <- pid_cov_all[pid %in% pid_common]
    setkey(pid_cov_common, pid)
    pid_cov_common <- pid_cov_common[.(sort(pid_common))]
    
    pid_cov_common[, Zopt := z_strict(lotr_z)]
    pid_cov_common[, Zrt  := z_strict(rt_log_mean)]
    
    outcome_blocks <- list(
      list(k=1,name="b",pid=pid1,y=mu_b),
      list(k=2,name="a",pid=pid2,y=mu_a),
      list(k=3,name="c",pid=pid3,y=mu_c),
      list(k=4,name="h",pid=pid4,y=mu_h)
    )
    
    blocks <- list()
    
    for(ob in outcome_blocks){
      
      pids <- as.character(ob$pid)
      y_full <- ob$y
      
      if (ob$name == "a") {
        keep_pid <- intersect(pids, pid_common_a)
      } else {
        keep_pid <- intersect(pids, pid_common)
      }
      
      keep <- pids %in% keep_pid
      
      if(sum(keep) < min_Nk) next
      
      y_rep <- y_full[,keep,drop=FALSE]
      
      K_all <- nrow(y_rep)
      Trep <- min(ex2_trep,K_all)
      
      set.seed(seed+ob$k)
      idx <- sort(sample.int(K_all,Trep))
      
      y_rep <- y_rep[idx,,drop=FALSE]
      
      y_rep <- z_within_draw(y_rep)
      
      pid_cov_k <- pid_cov_common[.(pids[keep])]
      
      blocks[[ob$name]] <- list(
        k = ob$k,
        pid = pids[keep],
        y = y_rep,
        Zopt = pid_cov_k$Zopt,
        Zrt  = pid_cov_k$Zrt
      )
    }
    
    if(length(blocks)==0) next
    
    Trep_joint <- min(sapply(blocks,function(x)nrow(x$y)))
    
    kid <- c()
    Zopt <- c()
    Zrt <- c()
    y_stack <- NULL
    
    for(b in blocks){
      
      Nk <- length(b$pid)
      
      kid  <- c(kid,rep(b$k,Nk))
      Zopt <- c(Zopt,b$Zopt)
      Zrt  <- c(Zrt,b$Zrt)
      
      yb <- b$y[1:Trep_joint,,drop=FALSE]
      
      y_stack <- if(is.null(y_stack)) yb else cbind(y_stack,yb)
    }
    
    data_list <- list(
      K=4L,
      Nobs=length(kid),
      Trep=Trep_joint,
      y_rep=y_stack,
      kid=as.integer(kid),
      Zopt=Zopt,
      Zrt=Zrt
    )
    
    msg("Sampling EX2 Stan model for treatment: ", tr)
    
    fit <- rstan::sampling(
      sm,
      data=data_list,
      iter=iter_val,
      warmup=warmup_val,
      chains=chains_val,
      seed=seed,
      control=list(
        adapt_delta=adapt_delta_val,
        max_treedepth=treedepth_val
      )
    )
    
    fit_file <- file.path(mod_dir,paste0("ex2_fit_",tr,".rds"))
    saveRDS(fit, fit_file)
    
    msg("Saved Stan fit: ", fit_file)
    
    post <- rstan::extract(fit)
    
    summ <- function(x){
      q <- quantile(x, c(0.025, 0.975), names = FALSE)
      c(
        median = median(x),
        mean   = mean(x),
        q025   = q[1],
        q975   = q[2],
        p_gt0  = mean(x > 0)
      )
    }
    
    pooled <- rbind(
      beta_opt_bar=summ(post$beta_opt_bar),
      beta_rt_bar=summ(post$beta_rt_bar)
    )
    
    tbl <- data.table(
      term=rownames(pooled),
      median=pooled[,"median"],
      mean=pooled[,"mean"],
      q025=pooled[,"q025"],
      q975=pooled[,"q975"],
      p_gt0=pooled[,"p_gt0"],
      dataset=ds,
      treatment=tr
    )
    
    coeff_file <- file.path(out_dir,paste0("ex2_coeffs_",tr,".csv"))
    fwrite(tbl, coeff_file)
    
    msg("Saved EX2 coefficient table: ", coeff_file)
    
    outputs[[tr]] <- tbl
  }
  
  invisible(outputs)
}