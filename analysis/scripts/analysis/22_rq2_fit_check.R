# ------------------------------------------------------------
# RQ2 Fit Check
#
# What this does (per dataset x treatment):
# 1) Load rq2_fit_sequences_<tr>.rds
# 2) Print convergence diagnostics (Rhat, ESS)
# 3) Report divergences and treedepth hits
# 4) Posterior predictive check for z (Normal model)
#
# Outputs:
#   data/clean/<ds>/output/rq2_fitcheck_<tr>.txt
#   data/clean/<ds>/output/rq2_ppc_<tr>.png
# ------------------------------------------------------------

library(data.table)
library(rstan)
library(ggplot2)

rq2_fit_check <- function(cfg) {
  
  ds   <- as.character(cfg$run$dataset)
  seed <- as.integer(cfg$run$seed)
  tr_vec <- unique(as.character(cfg$run$treatment))
  
  mod_dir <- path_mod_ds(ds)
  out_dir <- path_out_ds(ds)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  for (tr in tr_vec) {
    
    f_fit <- file.path(mod_dir, paste0("rq2_fit_sequences_", tr, ".rds"))
    stopifnot(file.exists(f_fit))
    
    f_txt  <- file.path(out_dir, paste0("rq2_fitcheck_", tr, ".txt"))
    f_ppc  <- file.path(out_dir, paste0("rq2_ppc_", tr, ".png"))
    
    if (should_skip(
      paths = c(f_txt, f_ppc),
      cfg   = cfg,
      type  = "output",
      label = paste0("RQ2 fit check (", ds, "/", tr, ")")
    )) next
    
    fit <- readRDS(f_fit)
    
    # ----------------------------
    # Convergence summary
    # ----------------------------
    summ <- summary(fit)$summary
    rhat_max <- max(summ[, "Rhat"], na.rm = TRUE)
    ess_min  <- min(summ[, "n_eff"], na.rm = TRUE)
    
    sampler_params <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
    
    divergences <- sum(sapply(sampler_params, function(x) sum(x[, "divergent__"])))
    treedepth_hits <- sum(sapply(sampler_params, function(x) sum(x[, "treedepth__"] == max(x[, "treedepth__"]))))
    
    sink(f_txt)
    cat("RQ2 Fit Diagnostics\n")
    cat("====================\n\n")
    cat("Dataset:", ds, "\n")
    cat("Treatment:", tr, "\n\n")
    
    cat("Max Rhat:", rhat_max, "\n")
    cat("Min n_eff:", ess_min, "\n\n")
    
    cat("Divergences:", divergences, "\n")
    cat("Max treedepth hits:", treedepth_hits, "\n\n")
    
    cat("Key parameters:\n")
    print(summary(fit, pars = c("alpha", "sigma", "sigma_u", "sigma_s"))$summary)
    
    sink()
    
    msg("Saved fit diagnostics: ", f_txt)
    
    # ----------------------------
    # Posterior predictive check
    # ----------------------------
    post <- rstan::extract(fit)
    
    # Observed z
    f_prep <- file.path(out_dir, paste0("rq2_prepared_", tr, ".csv"))
    stopifnot(file.exists(f_prep))
    d <- fread(f_prep)
    
    z_obs <- d$z
    
    # Simulate replicated z from posterior
    # z_rep = Normal(alpha + u_i + b_s, sigma)
    # For PPC, approximate using marginal sigma only
    alpha_draws <- post$alpha
    sigma_draws <- post$sigma
    
    n_rep <- min(1000, length(alpha_draws))
    idx <- seq_len(n_rep)
    
    z_rep <- rnorm(
      n = length(z_obs) * n_rep,
      mean = rep(alpha_draws[idx], each = length(z_obs)),
      sd   = rep(sigma_draws[idx], each = length(z_obs))
    )
    
    z_rep <- matrix(z_rep, ncol = n_rep)
    
    df_ppc <- data.table(
      z_obs = z_obs,
      z_rep_mean = rowMeans(z_rep)
    )
    
    p <- ggplot() +
      geom_density(aes(x = z_obs), data = df_ppc, color = "black") +
      geom_density(aes(x = z_rep_mean), data = df_ppc, color = "red") +
      labs(
        title = paste0("RQ2 PPC (", ds, "/", tr, ")"),
        x = "z",
        y = "Density"
      ) +
      theme_minimal()
    
    ggsave(f_ppc, p, width = 6, height = 4)
    
    msg("Saved PPC plot: ", f_ppc)
  }
  
  invisible(NULL)
}