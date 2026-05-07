# ============================================================
# 32_rq3_diagnostics.R
#
# PURPOSE
#   Diagnostics for RQ3 welfare loss models. Checks convergence
#   of all three fits and determines selected model per tr/tag:
#     - hurdle probability negligible -> gamma_only
#     - hurdle probability substantial -> primary (hurdle-Gamma)
#     - fallback if primary/gamma_only fail convergence -> alternative (Gaussian)
#   Produces PPC density comparison plot.
#
# Per preregistration:
#   - Gaussian always fitted as diagnostic alternative
#   - If hurdle probability negligible AND PPC shows no point
#     mass at zero, Gamma-only specification used
#   - Substantive conclusions based on stability across specs
#   - Convergence criteria: Rhat <= 1.01, ESS >= 400
#   - If selected model fails convergence, fall back to alternative
#     if it meets convergence; otherwise keep selected and flag deviation
#
# INPUT
#   path_mod/rq3_fit_sequences_<tr>_<tag>.rds        -- primary
#   path_mod/rq3_fit_sequences_<tr>_<tag>_gamma.rds  -- Gamma-only
#   path_mod/rq3_fit_sequences_<tr>_<tag>_alt.rds    -- Gaussian
#   path_mod/rq3_pid_levels_<tr>_<tag>.rds
#   path_mod/rq3_seq_levels_<tr>_<tag>.rds
#
# OUTPUT
#   path_out/rq3_diagnostics.csv
#   path_fig/rq3_ppc_<tr>_<tag>.png
# ============================================================

rq3_diagnostics <- function(cfg) {
  
  seed   <- as.integer(cfg$run$seed)
  tr_vec <- unique(as.character(cfg$run$treatment))
  tags   <- c("full", "confirmatory")
  
  hurdle_negligible_threshold <- 0.05
  rhat_threshold              <- 1.01
  ess_threshold               <- 400
  
  f_out <- file.path(path_out, "rq3_diagnostics.csv")
  if (should_skip(f_out, cfg, "output", "RQ3 diagnostics")) return(invisible(NULL))
  
  # ---- Helper: convergence summary for one fit ----
  convergence_summary <- function(fit) {
    summ   <- summary(fit)$summary
    sp     <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
    list(
      rhat_max       = max(summ[, "Rhat"], na.rm = TRUE),
      ess_min        = min(summ[, "n_eff"], na.rm = TRUE),
      divergences    = sum(sapply(sp, function(x) sum(x[, "divergent__"]))),
      treedepth_hits = sum(sapply(sp, function(x) {
        sum(x[, "treedepth__"] >= cfg$model$stan$rq3$treedepth)
      }))
    )
  }
  
  # ---- Helper: meets preregistered convergence criteria ----
  meets_convergence <- function(conv) {
    !is.na(conv$rhat_max) &&
      conv$rhat_max  <= rhat_threshold &&
      conv$ess_min   >= ess_threshold  &&
      conv$divergences == 0
  }
  
  all_rows <- list()
  
  for (tr in tr_vec) {
    for (tag in tags) {
      
      if (tag == "confirmatory" && !isTRUE(cfg$design$a_flags$betting_normative[[tr]])) next
      
      f_fit_p   <- file.path(path_mod, paste0("rq3_fit_sequences_", tr, "_", tag, ".rds"))
      f_fit_g   <- file.path(path_mod, paste0("rq3_fit_sequences_", tr, "_", tag, "_gamma.rds"))
      f_fit_alt <- file.path(path_mod, paste0("rq3_fit_sequences_", tr, "_", tag, "_alt.rds"))
      f_seq     <- file.path(path_mod, paste0("rq3_seq_levels_",    tr, "_", tag, ".rds"))
      
      if (!file.exists(f_fit_p)) {
        warning("RQ3 diagnostics: missing primary fit for tr='", tr, "', tag='", tag, "'. Skipping.")
        next
      }
      
      fit_p  <- readRDS(f_fit_p)
      post_p <- rstan::extract(fit_p)
      conv_p <- convergence_summary(fit_p)
      
      msg("RQ3 diagnostics primary (", tr, "/", tag, "):",
          " Rhat_max=", round(conv_p$rhat_max, 4),
          " ESS_min=",  round(conv_p$ess_min, 0),
          " divergences=", conv_p$divergences,
          " treedepth_hits=", conv_p$treedepth_hits)
      
      # ---- Hurdle probability ----
      stopifnot(!is.null(post_p$a0))
      hurdle_mean       <- mean(plogis(post_p$a0))
      hurdle_negligible <- hurdle_mean < hurdle_negligible_threshold
      
      msg("RQ3 hurdle probability (", tr, "/", tag, "):",
          " mean_omega=", round(hurdle_mean, 4),
          if (hurdle_negligible) " [NEGLIGIBLE -- Gamma-only preferred]"
          else " [substantial -- hurdle-Gamma appropriate]")
      
      # ---- Gamma-only convergence ----
      conv_g       <- list(rhat_max = NA_real_, ess_min = NA_real_,
                           divergences = NA_integer_, treedepth_hits = NA_integer_)
      gamma_exists <- file.exists(f_fit_g)
      
      if (gamma_exists) {
        fit_g  <- readRDS(f_fit_g)
        conv_g <- convergence_summary(fit_g)
        msg("RQ3 diagnostics Gamma-only (", tr, "/", tag, "):",
            " Rhat_max=", round(conv_g$rhat_max, 4),
            " ESS_min=",  round(conv_g$ess_min, 0),
            " divergences=", conv_g$divergences)
      }
      
      # ---- Alternative convergence + Spearman stability ----
      conv_alt     <- list(rhat_max = NA_real_, ess_min = NA_real_,
                           divergences = NA_integer_, treedepth_hits = NA_integer_)
      spearman_seq <- NA_real_
      alt_exists   <- file.exists(f_fit_alt)
      
      if (alt_exists) {
        fit_alt  <- readRDS(f_fit_alt)
        post_alt <- rstan::extract(fit_alt)
        conv_alt <- convergence_summary(fit_alt)
        
        msg("RQ3 diagnostics alternative (", tr, "/", tag, "):",
            " Rhat_max=", round(conv_alt$rhat_max, 4),
            " ESS_min=",  round(conv_alt$ess_min, 0),
            " divergences=", conv_alt$divergences)
        
        # Spearman: compare selected non-Gaussian model vs Gaussian
        ref_post <- if (gamma_exists && hurdle_negligible) {
          rstan::extract(readRDS(f_fit_g))
        } else {
          post_p
        }
        
        if (!is.null(ref_post$mu_c) && !is.null(post_alt$mu_c)) {
          mu_c_ref     <- apply(ref_post$mu_c, 2, mean)
          mu_c_alt     <- apply(post_alt$mu_c, 2, mean)
          spearman_seq <- cor(mu_c_ref, mu_c_alt, method = "spearman")
          msg("RQ3 sequence ordering stability (", tr, "/", tag, "):",
              " Spearman rho=", round(spearman_seq, 4))
        }
      }
      
      # ---- Selected model ----
      selected_model <- if (hurdle_negligible && gamma_exists) {
        "gamma_only"
      } else {
        "primary"
      }
      
      adequate <- if (!is.na(spearman_seq)) spearman_seq > 0.90 else NA
      
      # ---- Convergence check on selected model ----
      conv_selected <- if (selected_model == "gamma_only") conv_g else conv_p
      selected_converged <- meets_convergence(conv_selected)
      
      # ---- Convergence-based fallback ----
      fallback_used <- FALSE
      
      if (!selected_converged && selected_model != "alternative") {
        
        if (alt_exists && meets_convergence(conv_alt)) {
          msg("RQ3 (", tr, "/", tag, "): selected model fails convergence criteria",
              " (Rhat=", round(conv_selected$rhat_max, 4),
              ", ESS=", round(conv_selected$ess_min, 0), ")",
              " -- falling back to alternative (Gaussian)")
          selected_model   <- "alternative"
          adequate         <- TRUE
          fallback_used    <- TRUE
        } else {
          msg("RQ3 (", tr, "/", tag, "): selected model fails convergence criteria",
              " (Rhat=", round(conv_selected$rhat_max, 4),
              ", ESS=", round(conv_selected$ess_min, 0), ")",
              " -- alternative also fails; keeping selected [PREREGISTRATION DEVIATION]")
          fallback_used <- FALSE
        }
      }
      
      msg("RQ3 (", tr, "/", tag, "): selected=", selected_model,
          " | converged=", selected_converged,
          " | fallback=", fallback_used,
          if (!is.na(adequate)) paste0(" | adequate=", adequate) else "")
      # ---- PPC density plot ----
      f_ppc <- file.path(path_fig, paste0("rq3_ppc_", tr, "_", tag, ".png"))
      if (alt_exists && !should_skip(f_ppc, cfg, "output",
                                     paste0("RQ3 PPC plot (", tr, "/", tag, ")"))) {
        
        if (selected_model == "alternative") {
          ref_post_plot <- post_alt
          ref_label     <- "Gaussian (selected)"
          df_ppc <- data.table(mu_c   = as.vector(ref_post_plot$mu_c),
                               source = ref_label)
        } else {
          ref_post_plot <- if (selected_model == "gamma_only") {
            rstan::extract(readRDS(f_fit_g))
          } else {
            post_p
          }
          ref_label <- switch(selected_model,
                              gamma_only = "Gamma-only (selected)",
                              primary    = "hurdle-Gamma (selected)")
          df_ppc <- rbind(
            data.table(mu_c = as.vector(ref_post_plot$mu_c), source = ref_label),
            data.table(mu_c = as.vector(post_alt$mu_c),      source = "Gaussian (alternative)")
          )
        }
        
        if (!is.null(df_ppc) && nrow(df_ppc) > 0) {
          
          n_sources  <- length(unique(df_ppc$source))
          col_vals   <- if (n_sources == 1) "#2166AC" else c("#2166AC", "#C0392B")
          fill_vals  <- if (n_sources == 1) "#2166AC" else c("#2166AC", "#C0392B")
          src_levels <- unique(df_ppc$source)
          df_ppc[, source := factor(source, levels = src_levels)]
          
          annot_label <- paste0(
            sprintf("Spearman \u03c1 = %.3f", ifelse(is.na(spearman_seq), 0, spearman_seq)),
            sprintf("\nHurdle mean = %.3f", hurdle_mean)
          )
          
          p <- ggplot(df_ppc, aes(x = mu_c, colour = source, fill = source)) +
            geom_density(alpha = 0.15, linewidth = 0.8) +
            scale_colour_manual(values = setNames(col_vals,  src_levels), name = NULL) +
            scale_fill_manual(  values = setNames(fill_vals, src_levels), name = NULL) +
            scale_x_continuous(
              name   = expression(hat(mu)[s]^c ~ "(sequence-level welfare loss)"),
              labels = scales::percent_format(accuracy = 1)
            ) +
            scale_y_continuous(name = "Density") +
            annotate("text", x = Inf, y = Inf,
                     label  = annot_label,
                     hjust  = 1.05, vjust = 1.5,
                     size   = 3, colour = "grey40",
                     family = "mono") +
            theme_classic(base_size = 11) +
            theme(
              legend.position    = "bottom",
              legend.text        = element_text(size = 9),
              panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
              axis.title         = element_text(size = 10),
              plot.margin        = margin(8, 12, 8, 8)
            )
          
          ggsave(f_ppc, p, width = 6, height = 4, dpi = 300)
          msg("Saved: ", f_ppc)
        }
      }
      
      all_rows[[paste(tr, tag, sep = "_")]] <- data.table(
        treatment          = tr,
        tag                = tag,
        hurdle_mean        = round(hurdle_mean, 4),
        hurdle_negligible  = hurdle_negligible,
        selected_model     = selected_model,
        fallback_used      = fallback_used,
        spearman_seq       = round(spearman_seq, 4),
        adequate           = adequate,
        # Primary
        rhat_max_primary   = round(conv_p$rhat_max, 4),
        ess_min_primary    = round(conv_p$ess_min, 0),
        divergences_p      = conv_p$divergences,
        treedepth_p        = conv_p$treedepth_hits,
        # Gamma-only
        gamma_exists       = gamma_exists,
        rhat_max_gamma     = round(conv_g$rhat_max, 4),
        ess_min_gamma      = round(conv_g$ess_min, 0),
        divergences_g      = conv_g$divergences,
        # Alternative
        alt_exists         = alt_exists,
        rhat_max_alt       = round(conv_alt$rhat_max, 4),
        ess_min_alt        = round(conv_alt$ess_min, 0),
        divergences_alt    = conv_alt$divergences
      )
    }
  }
  
  if (length(all_rows) > 0) {
    fwrite(rbindlist(all_rows, fill = TRUE), f_out)
    msg("Saved: ", f_out)
  }
  
  invisible(all_rows)
}