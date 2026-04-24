# ============================================================
# 22_rq2_diagnostics.R
#
# PURPOSE
#   Convergence diagnostics and posterior predictive check for
#   RQ2 Gaussian likelihood. Uses Kolmogorov-Smirnov statistic
#   to assess distributional fit of standardized stake deviations.
#   Flags Beta-Binomial robustness if PPC is inadequate.
#
# ADEQUACY CRITERION:
#   PPC p-value based on KS statistic comparing observed z to
#   posterior predictive replications. Adequate if ppc_p in
#   [p_lo, p_hi] from cfg$model$ppc$rq1_interval.
#   Also reports pct_allin as a descriptive statistic.
#
# INPUT
#   path_mod/rq2_fit_sequences_<tr>_<tag>.rds
#   path_mod/rq2_prepared_<tr>_<tag>.rds
#
# OUTPUT
#   path_out/rq2_diagnostics.csv      -- adequacy flag per tr/tag
#   path_fig/rq2_ppc_<tr>_<tag>.png  -- PPC density plot
#
# CALL ORDER IN PIPELINE:
#   rq2_stan(cfg)                    -- primary fits
#   rq2_diagnostics(cfg)             -- adequacy check
#   rq2_stan(cfg, robustness = TRUE) -- BB if needed, no-op otherwise
# ============================================================

rq2_diagnostics <- function(cfg) {
  
  seed   <- as.integer(cfg$run$seed)
  tr_vec <- unique(as.character(cfg$run$treatment))
  tags   <- c("full", "confirmatory")
  e      <- as.numeric(cfg$design$seq$endowment)
  
  K_ppc <- as.integer(cfg$model$ppc$rq1_k)
  p_lo  <- as.numeric(cfg$model$ppc$rq1_interval[1])
  p_hi  <- as.numeric(cfg$model$ppc$rq1_interval[2])
  stopifnot(is.finite(p_lo), is.finite(p_hi), p_lo > 0, p_hi < 1, p_lo < p_hi)
  
  f_out <- file.path(path_out, "rq2_diagnostics.csv")
  
  if (should_skip(f_out, cfg, "output", "RQ2 diagnostics")) return(invisible(NULL))
  
  all_rows <- list()
  
  for (tr in tr_vec) {
    for (tag in tags) {
      
      if (tag == "confirmatory" && !isTRUE(cfg$design$a_flags$betting_normative[[tr]])) next
      
      f_fit  <- file.path(path_mod, paste0("rq2_fit_sequences_", tr, "_", tag, ".rds"))
      f_prep <- file.path(path_mod, paste0("rq2_prepared_",      tr, "_", tag, ".rds"))
      
      if (!file.exists(f_fit) || !file.exists(f_prep)) {
        warning("RQ2 diagnostics: missing artifacts for tr='", tr, "', tag='", tag, "'. Skipping.")
        next
      }
      
      fit  <- readRDS(f_fit)
      prep <- readRDS(f_prep)
      stopifnot(is.data.frame(prep), "stake" %in% names(prep), "z" %in% names(prep))
      
      # ---- Convergence diagnostics ----
      summ           <- summary(fit)$summary
      rhat_max       <- max(summ[, "Rhat"], na.rm = TRUE)
      ess_min        <- min(summ[, "n_eff"], na.rm = TRUE)
      sampler_params <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
      divergences    <- sum(sapply(sampler_params, function(x) sum(x[, "divergent__"])))
      treedepth_hits <- sum(sapply(sampler_params, function(x) {
        sum(x[, "treedepth__"] >= cfg$model$stan$rq2$treedepth)
      }))
      
      msg("RQ2 diagnostics (", tr, "/", tag, "):",
          " Rhat_max=", round(rhat_max, 4),
          " ESS_min=",  round(ess_min, 0),
          " divergences=", divergences,
          " treedepth_hits=", treedepth_hits)
      
      # ---- Descriptive: all-in proportion ----
      pct_allin <- mean(prep$stake >= e, na.rm = TRUE)
      
      # ---- PPC: KS statistic ----
      # Observed KS statistic vs standard normal
      post        <- rstan::extract(fit)
      alpha_draws <- as.numeric(post$alpha)
      sigma_draws <- as.numeric(post$sigma)
      z_obs       <- prep$z
      
      K_use <- min(K_ppc, length(alpha_draws))
      set.seed(seed)
      k_idx <- sort(sample.int(length(alpha_draws), K_use, replace = FALSE))
      
      D_obs <- ks.test(z_obs, "pnorm", mean = mean(z_obs), sd = sd(z_obs))$statistic
      
      D_rep <- sapply(k_idx, function(k) {
        z_rep <- rnorm(length(z_obs), mean = alpha_draws[k], sd = sigma_draws[k])
        ks.test(z_obs, z_rep)$statistic
      })
      
      ppc_p    <- mean(D_rep >= D_obs)
      adequate <- (ppc_p >= p_lo && ppc_p <= p_hi)
      
      if (!adequate) {
        warning("RQ2 PPC (", tr, "/", tag, "): INADEQUATE -- BB robustness will be fitted",
                " | ppc_p=", round(ppc_p, 4),
                " | KS_obs=", round(D_obs, 4),
                " | pct_allin=", round(pct_allin, 4))
      } else {
        msg("RQ2 PPC (", tr, "/", tag, "): adequate",
            " | ppc_p=", round(ppc_p, 4),
            " | KS_obs=", round(D_obs, 4),
            " | pct_allin=", round(pct_allin, 4))
      }
      
      # ---- PPC density plot ----
      f_ppc <- file.path(path_fig, paste0("rq2_ppc_", tr, "_", tag, ".png"))
      if (!should_skip(f_ppc, cfg, "output",
                       paste0("RQ2 PPC plot (", tr, "/", tag, ")"))) {
        
        set.seed(seed + 1L)
        k_plot <- sort(sample.int(length(alpha_draws), min(200L, K_use), replace = FALSE))
        
        z_rep_mat <- sapply(k_plot, function(k) {
          rnorm(length(z_obs), mean = alpha_draws[k], sd = sigma_draws[k])
        })
        
        df_obs <- data.table(z = z_obs,               source = "observed")
        df_rep <- data.table(z = as.vector(z_rep_mat), source = "replicated")
        df_ppc <- rbind(df_obs, df_rep)
        
        p <- ggplot(df_ppc, aes(x = z, colour = source)) +
          geom_density(linewidth = 0.7) +
          scale_colour_manual(values = c(observed = "black", replicated = "#C0392B"),
                              name = NULL) +
          annotate("text", x = Inf, y = Inf,
                   label = sprintf("KS p = %.3f\nAll-in: %.1f%%\nAdequate: %s",
                                   ppc_p, 100 * pct_allin,
                                   if (adequate) "Yes" else "No"),
                   hjust = 1.1, vjust = 1.5, size = 3, colour = "grey30") +
          labs(title = paste0("RQ2 PPC (", tr, "/", tag, ")"),
               x = "z (standardized stake deviation)", y = "Density") +
          theme_classic(base_size = 11) +
          theme(legend.position = "bottom")
        
        ggsave(f_ppc, p, width = 6, height = 4, dpi = 300)
        msg("Saved: ", f_ppc)
      }
      
      all_rows[[paste(tr, tag, sep = "_")]] <- data.table(
        treatment      = tr,
        tag            = tag,
        rhat_max       = round(rhat_max, 4),
        ess_min        = round(ess_min, 0),
        divergences    = divergences,
        treedepth_hits = treedepth_hits,
        pct_allin      = round(pct_allin, 4),
        KS_obs         = round(D_obs, 4),
        KS_rep_median  = round(median(D_rep), 4),
        ppc_p          = round(ppc_p, 4),
        adequate       = adequate,
        p_lo           = p_lo,
        p_hi           = p_hi
      )
    }
  }
  
  if (length(all_rows) > 0) {
    fwrite(rbindlist(all_rows), f_out)
    msg("Saved: ", f_out)
  }
  
  invisible(all_rows)
}