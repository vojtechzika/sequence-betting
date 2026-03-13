# ============================================================
# scripts/analysis/55_ex1_3_pure_approximation.R
#
# EX1.3 Pure-sequence approximation to the GHI
#
# Uses:
#   chi_i draws from:
#       data/clean/<ds>/models/ex1_1_<tr>_participants.rds
#
# Outputs:
#   data/clean/<ds>/output/ex1_3_rho_<tr>.csv
#   data/clean/<ds>/output/ex1_3_groups_<tr>.csv
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
})

ex1_3_ghi_vs_pure <- function(cfg) {
  
  stopifnot(is.list(cfg), !is.null(cfg$run), !is.null(cfg$design))
  stopifnot(!is.null(cfg$run$dataset))
  
  ds <- as.character(cfg$run$dataset)
  
  design <- cfg$design
  
  # ------------------------------------------------------------
  # treatments
  # ------------------------------------------------------------
  
  tr_vec <- cfg$run$treatment
  
  if (is.null(tr_vec) || length(tr_vec)==0 || all(!nzchar(as.character(tr_vec)))) {
    stop("EX1.3 requires cfg$run$treatment")
  }
  
  tr_vec <- unique(as.character(tr_vec))
  
  # ------------------------------------------------------------
  # labels
  # ------------------------------------------------------------
  
  lab_heads <- as.character(design$seq$side_labels$heads)
  lab_tails <- as.character(design$seq$side_labels$tails)
  
  pure_heads <- as.character(design$seq$anchor_labels$pure_heads)
  pure_tails <- as.character(design$seq$anchor_labels$pure_tails)
  
  # ------------------------------------------------------------
  # paths
  # ------------------------------------------------------------
  
  mod_dir <- path_mod_ds(ds)
  out_dir <- path_out_ds(ds)
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  f_master <- file.path(path_clean_ds(ds), "master_sequences.csv")
  stopifnot(file.exists(f_master))
  
  master <- fread(f_master, encoding="UTF-8")
  
  req <- c("pid","treat","seq","stake","side")
  miss <- setdiff(req,names(master))
  if(length(miss)>0)
    stop("master_sequences.csv missing: ",paste(miss,collapse=", "))
  
  master[, pid := as.character(pid)]
  master[, treat := as.character(treat)]
  master[, seq := as.character(seq)]
  master[, stake := as.numeric(stake)]
  master[, side := as.character(side)]
  
  master[is.na(stake), stake := 0]
  
  outputs <- list()
  
  for(tr in tr_vec){
    
    # ------------------------------------------------------------
    # load chi draws
    # ------------------------------------------------------------
    
    f_chi <- file.path(mod_dir,paste0("ex1_1_",tr,"_participants.rds"))
    
    if(!file.exists(f_chi))
      stop("EX1.3: Missing EX1.1 participants RDS: ",f_chi)
    
    chi_obj <- readRDS(f_chi)
    
    pid_levels <- chi_obj$pid_levels
    if(is.null(pid_levels)) pid_levels <- chi_obj$pid
    
    chi_draws <- chi_obj$chi_draws
    if(is.null(chi_draws)) chi_draws <- chi_obj$chi_rep
    
    stopifnot(is.matrix(chi_draws))
    stopifnot(ncol(chi_draws)==length(pid_levels))
    
    K <- nrow(chi_draws)
    
    # ------------------------------------------------------------
    # pure-sequence behavior
    # ------------------------------------------------------------
    
    d <- master[
      treat==tr &
        pid %in% pid_levels &
        seq %in% c(pure_heads,pure_tails)
    ]
    
    d <- d[stake>0]
    
    d <- d[side %in% c(lab_heads,lab_tails)]
    
    d[,h := as.integer(side==lab_heads)]
    
    chk <- d[, .N, by=.(pid,seq)]
    
    bad_chk <- chk[N!=1]
    
    if(nrow(bad_chk)>0){
      stop(
        "EX1.3 expected exactly one bet per participant on pure sequences.\n",
        paste(capture.output(print(head(bad_chk,10))),collapse="\n")
      )
    }
    
    pure_wide <- dcast(d,pid~seq,value.var="h")
    
    pure_wide <- pure_wide[
      is.finite(get(pure_heads)) &
        is.finite(get(pure_tails))
    ]
    
    pid_keep <- intersect(pid_levels,pure_wide$pid)
    
    if(length(pid_keep)==0)
      stop("EX1.3: no valid participants for treatment ",tr)
    
    idx_keep <- match(pid_keep,pid_levels)
    
    chi_draws_keep <- chi_draws[,idx_keep,drop=FALSE]
    
    setkey(pure_wide,pid)
    
    chi_pure <- pure_wide[.(pid_keep),
                          as.integer(get(pure_heads)==1L) -
                            as.integer(get(pure_tails)==1L)
    ]
    
    # ------------------------------------------------------------
    # safe correlation
    # ------------------------------------------------------------
    
    safe_cor <- function(x,y){
      
      ok <- is.finite(x) & is.finite(y)
      
      if(sum(ok)<3) return(NA_real_)
      
      sx <- sd(x[ok])
      sy <- sd(y[ok])
      
      if(!is.finite(sx) || !is.finite(sy) || sx<=0 || sy<=0)
        return(NA_real_)
      
      suppressWarnings(cor(x[ok],y[ok]))
    }
    
    rho_draws <- vapply(
      seq_len(K),
      function(k) safe_cor(chi_pure,chi_draws_keep[k,]),
      numeric(1)
    )
    
    rho_draws <- rho_draws[is.finite(rho_draws)]
    
    if (length(rho_draws) < 10) {
      warning("Too few finite rho draws for treatment ", tr)
      
      rho_tbl <- data.table(
        dataset = ds,
        treatment = tr,
        N = length(pid_keep),
        K = K,
        K_finite = length(rho_draws),
        rho_median = NA_real_,
        rho_mean   = NA_real_,
        rho_q025   = NA_real_,
        rho_q975   = NA_real_,
        p_rho_gt0  = NA_real_
      )
    } else {
      rho_tbl <- data.table(
        dataset = ds,
        treatment = tr,
        N = length(pid_keep),
        K = K,
        K_finite = length(rho_draws),
        rho_median = median(rho_draws),
        rho_mean   = mean(rho_draws),
        rho_q025   = quantile(rho_draws, 0.025),
        rho_q975   = quantile(rho_draws, 0.975),
        p_rho_gt0  = mean(rho_draws > 0)
      )
    }
    
    f_rho <- file.path(out_dir,paste0("ex1_3_rho_",tr,".csv"))
    
    fwrite(rho_tbl,f_rho)
    
    msg("Saved: ",f_rho)
    
    # ------------------------------------------------------------
    # grouped chi summaries
    # ------------------------------------------------------------
    
    groups <- sort(unique(chi_pure))
    
    grp_rows <- list()
    
    for(g in groups){
      
      idxg <- which(chi_pure==g)
      
      if(length(idxg)==0) next
      
      chi_g_draws <- rowMeans(
        chi_draws_keep[,idxg,drop=FALSE],
        na.rm=TRUE
      )
      
      grp_rows[[as.character(g)]] <- data.table(
        dataset = ds,
        treatment = tr,
        chi_pure = g,
        n = length(idxg),
        chi_median = median(chi_g_draws),
        chi_mean   = mean(chi_g_draws),
        chi_q025   = quantile(chi_g_draws,0.025),
        chi_q975   = quantile(chi_g_draws,0.975)
      )
      
    }
    
    grp_tbl <- rbindlist(grp_rows,use.names=TRUE,fill=TRUE)
    
    setorder(grp_tbl,chi_pure)
    
    f_grp <- file.path(out_dir,paste0("ex1_3_groups_",tr,".csv"))
    
    fwrite(grp_tbl,f_grp)
    
    msg("Saved: ",f_grp)
    
    outputs[[tr]] <- list(
      rho=rho_tbl,
      groups=grp_tbl
    )
    
  }
  
  invisible(outputs)
}

# Example
# ex1_3_ghi_vs_pure(cfg)