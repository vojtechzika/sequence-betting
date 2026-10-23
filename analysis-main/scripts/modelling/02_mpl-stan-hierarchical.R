# ============================================================
# 02_mpl-stan-hierarchical.R
#
# PURPOSE
#   Alternative to the preregistered individual-level Holt-Laury model
#   (01_mpl-stan-bernoulli.R). The choice model is the same (Bernoulli
#   with a logit link in lambda_i times the expected-utility difference
#   between lotteries A and B, CRRA utility), but the participant
#   parameters are partially pooled: r_i ~ Normal(mu_r, sigma_r) and
#   log(lambda_i) ~ Normal(mu_ll, sigma_ll), with hyperpriors centred on
#   the preregistered prior (cfg$model$stan$mpl_hierarchical$prior).
#   Pooling shrinks the estimates of participants whose choices are
#   weakly informative (inconsistent patterns) toward the population.
#
# WHAT IT DOES
#   The participant parameters are integrated out on a fixed grid (see
#   scripts/stan/mpl-hierarchical.stan), so Stan samples only the four
#   hyperparameters. The likelihood of each distinct choice pattern on the
#   grid is computed here and passed as data. Draws of r_i given the
#   hyperparameters are generated in Stan from the exact conditional
#   posterior on the r grid. This avoids the sampling problems of
#   participants who choose the risky lottery in every row (see
#   01_mpl-stan-bernoulli.R) and needs no replacement of draws.
#   As in the individual-level model, the model is fitted separately for
#   each treatment, so that a participant's estimate reflects only
#   within-treatment information (the Holt-Laury task follows the sequence
#   experiment and may be affected by the incentive frame). Because a
#   participant's estimate now depends on the other participants in the
#   fit, the Holt-Laury-consistent sample is estimated in its own fit, in
#   which the population is estimated from the consistent participants of
#   the treatment only.
#
# INPUT
#   path_src/holt-laury.csv               -- MPL choices c1..cK
#   path_src/participants.csv             -- treatment and sex
#   path_der/switching-point-r.csv        -- inconsistency flag per pid
#   scripts/stan/mpl-hierarchical.stan
#
# OUTPUT (model folder path_mod/mpl-hierarchical/)
#   fit-<tag>.rds             -- Stan fit objects
#   r-draws-<tag>.rds         -- list(pid, r_draws): posterior draws of r_i
#   scored-<tag>.csv          -- r_mean, r_median, inconsistent per pid
#   hyperparameters-<tag>.csv -- posterior mean, sd and 95% interval of
#                                mu_r, sigma_r, mu_ll, sigma_ll
#   diagnostics.csv           -- convergence checks per fit and the largest
#                                posterior mass of any participant in the
#                                outermost cells of the r grid
#   model-info.csv            -- model property read by the a* step
#                                (partial_pooling = TRUE)
#   inputs.csv, run-info.txt
#
# TAGS AND FILE NAMES
#   <treatment>               all participants of the treatment
#   <treatment>_consistent    Holt-Laury-consistent participants of the
#                             treatment (if consistent_only)
#   The file layout of draws and scored files is that of the
#   individual-level model.
#
# NOTES
#   - If the model is listed in cfg$run$prior_sensitivity, a second pass
#     refits it with the broadened priors (prior_broad in settings.R) and
#     writes a complete separate folder with the suffix -broad
#   - There is no fit that combines treatments
#   - The grid covers r in [-8, 4] (step 0.02) and log lambda in [-2, 8]
#     (step 0.2); sigma_r >= 0.05 and sigma_ll >= 0.2 keep the population
#     densities resolved by these steps
#   - The model is not used downstream unless it is set in selected-models.csv;
#     its model-info.csv records partial_pooling = TRUE, so the a* step
#     stops unless its consistent-sample path reads the consistent fit
#   - The model folder is built as mpl-hierarchical.tmp and renamed when
#     complete; an existing folder is never overwritten
#   - Convergence criteria are cfg$model$convergence
# ============================================================

options(mc.cores = parallel::detectCores())
rstan::rstan_options(auto_write = FALSE)

# Log-likelihood of one choice pattern on the (r, log lambda) grid
#   y: choices (1 = safe, 0 = risky) in the rows of the task
#   p: probability of the high payoff in each row
mpl_pattern_loglik <- function(y, p, A1, A2, B1, B2, r_grid, ll_grid) {
  crra <- function(x, r) ifelse(abs(r - 1) < 1e-8, log(x), x^(1 - r) / (1 - r))
  d1  <- crra(A1, r_grid) - crra(B1, r_grid)
  d2  <- crra(A2, r_grid) - crra(B2, r_grid)
  lam <- exp(ll_grid)
  ll  <- matrix(0, length(r_grid), length(ll_grid))
  for (k in seq_along(p)) {
    z     <- outer(p[k] * d1 + (1 - p[k]) * d2, lam)
    log_s <- -(pmax(-z, 0) + log1p(exp(-abs(z))))     # log logistic(z)
    ll    <- ll + (if (y[k] == 1L) log_s else log_s - z)
  }
  ll
}

mpl_stan_hierarchical <- function(cfg) {
  for (pass in prior_passes(cfg, "mpl_hierarchical")) mpl_stan_hierarchical_pass(cfg, pass)
  invisible(NULL)
}

mpl_stan_hierarchical_pass <- function(cfg, pass) {
  
  run <- model_begin(paste0("mpl-hierarchical", pass$suffix))
  if (is.null(run)) return(invisible(NULL))
  
  seed            <- as.integer(cfg$run$seed)
  edge_mass_max   <- 0.01     # largest tolerated posterior mass in the outermost r cells
  consistent_only <- isTRUE(cfg$run$consistent_only)
  design          <- cfg$design
  treatments      <- names(design$seq$treatments)
  
  label_safe  <- toupper(trimws(as.character(design$mpl$label_safe)))
  label_risky <- toupper(trimws(as.character(design$mpl$label_risky)))
  stopifnot(length(label_safe) == 1L, length(label_risky) == 1L,
            nzchar(label_safe), nzchar(label_risky), label_safe != label_risky)
  
  # ---- Design of the MPL ----
  K_cfg  <- as.integer(design$mpl$K)
  A_high <- as.numeric(design$mpl$A_high)
  A_low  <- as.numeric(design$mpl$A_low)
  B_high <- as.numeric(design$mpl$B_high)
  B_low  <- as.numeric(design$mpl$B_low)
  p_sched <- as.numeric(design$mpl$p)
  stopifnot(K_cfg > 0L, is.finite(A_high), is.finite(A_low),
            is.finite(B_high), is.finite(B_low),
            A_high > 0, A_low > 0, B_high > 0, B_low > 0,
            length(p_sched) == K_cfg, all(is.finite(p_sched)),
            all(p_sched > 0), all(p_sched < 1), !is.unsorted(p_sched))
  
  # ---- Sampler settings and hyperpriors ----
  st <- cfg$model$stan$mpl_hierarchical
  if (is.null(st)) stop("cfg$model$stan$mpl_hierarchical is missing.", call. = FALSE)
  iter_val        <- as.integer(st$iter)
  warmup_val      <- as.integer(st$warmup)
  chains_val      <- as.integer(st$chains)
  adapt_delta_val <- as.numeric(st$adapt_delta)
  treedepth_val   <- as.integer(st$treedepth)
  pr              <- pass$prior
  stopifnot(iter_val > 0L, warmup_val >= 0L, chains_val > 0L,
            is.finite(adapt_delta_val), adapt_delta_val > 0, adapt_delta_val < 1,
            treedepth_val > 0L,
            is.finite(pr$mu_r_mean), pr$mu_r_sd > 0, pr$sigma_r_sd > 0,
            pr$lambda_mean > 0, pr$mu_ll_sd > 0, pr$sigma_ll_sd > 0)
  
  # ---- Grids ----
  r_grid  <- seq(-8, 4, by = 0.02)
  ll_grid <- seq(-2, 8, by = 0.2)
  R <- length(r_grid); L <- length(ll_grid)
  dr <- 0.02; dl <- 0.2
  
  # ---- Inputs ----
  stan_file <- here::here("scripts", "stan", "mpl-hierarchical.stan")
  f_mpl     <- file.path(path_src, "holt-laury.csv")
  f_par     <- file.path(path_src, "participants.csv")
  f_sw      <- file.path(path_der, "switching-point-r.csv")
  stopifnot(file.exists(stan_file), file.exists(f_mpl),
            file.exists(f_par), file.exists(f_sw))
  
  mpl <- fread(f_mpl, encoding = "UTF-8")
  par <- fread(f_par, encoding = "UTF-8")
  sw  <- fread(f_sw,  encoding = "UTF-8")
  mpl[, pid := as.character(pid)]
  par[, pid := as.character(pid)]
  sw[,  pid := as.character(pid)]
  
  if (anyDuplicated(par$pid)) stop("Repeated pid in participants.csv.")
  if (!all(mpl$pid %in% par$pid)) stop("MPL participants missing from participants.csv.")
  par[, sex_label := fcase(toupper(trimws(sex)) == "F", "Women",
                           toupper(trimws(sex)) == "M", "Men",
                           default = NA_character_)]
  mpl <- merge(mpl, par[, .(pid, treat = as.character(treat), sex_label)], by = "pid")
  if (!all(mpl$treat %in% treatments)) stop("Unknown treatment in MPL participants.")
  if (!all(mpl$pid %in% sw$pid)) stop("MPL participants missing from switching-point-r.csv.")
  mpl <- merge(mpl, sw[, .(pid, inconsistent = as.integer(inconsistent))], by = "pid")
  setorder(mpl, pid)
  
  choice_cols <- grep("^c\\d+$", names(mpl), value = TRUE)
  choice_cols <- choice_cols[order(as.integer(sub("^c", "", choice_cols)))]
  stopifnot(length(choice_cols) == K_cfg)
  
  sm <- rstan::stan_model(stan_file)
  
  # ---- Fit one sample ----
  fit_one <- function(d, tag, treatment, sample) {
    
    if (nrow(d) == 0L) stop("No MPL participants for tag '", tag, "'.")
    
    # Choices coded 1 = safe, 0 = risky
    ch <- toupper(trimws(as.matrix(d[, ..choice_cols])))
    if (!all(ch %in% c(label_safe, label_risky))) {
      stop("MPL choices contain labels other than the configured safe/risky labels.")
    }
    ymat <- matrix(as.integer(ch == label_safe), nrow = nrow(d))
    
    # Distinct choice patterns and their likelihood on the grid, scaled to a
    # maximum of one within each pattern
    key  <- apply(ymat, 1L, paste, collapse = "")
    ukey <- unique(key)
    pat  <- match(key, ukey)
    ypat <- ymat[match(ukey, key), , drop = FALSE]
    P    <- length(ukey)
    lik  <- matrix(0, nrow = P * R, ncol = L)
    for (j in seq_len(P)) {
      ll <- mpl_pattern_loglik(ypat[j, ], p_sched, A_high, A_low, B_high, B_low,
                               r_grid, ll_grid)
      lik[((j - 1L) * R + 1L):(j * R), ] <- exp(ll - max(ll))
    }
    
    data_list <- list(
      N = nrow(ymat), P = P, R = R, L = L,
      pat = pat, n_pat = tabulate(pat, nbins = P),
      r_grid = r_grid, ll_grid = ll_grid, dr = dr, dl = dl, lik = lik,
      h_mu_r_mean = pr$mu_r_mean, h_mu_r_sd = pr$mu_r_sd,
      h_sigma_r_sd = pr$sigma_r_sd, h_lambda_mean = pr$lambda_mean,
      h_mu_ll_sd = pr$mu_ll_sd, h_sigma_ll_sd = pr$sigma_ll_sd
    )
    
    # Initial values near the preregistered prior
    init_fun <- function() list(
      mu_r = pr$mu_r_mean + rnorm(1, 0, 0.05),
      sigma_r = 0.8 * exp(rnorm(1, 0, 0.1)),
      mu_ll = log(pr$lambda_mean) + rnorm(1, 0, 0.05),
      sigma_ll = 0.8 * exp(rnorm(1, 0, 0.1))
    )
    
    fit <- rstan::sampling(
      sm, data = data_list,
      iter = iter_val, warmup = warmup_val, chains = chains_val,
      seed = seed, init = init_fun,
      control = list(adapt_delta = adapt_delta_val, max_treedepth = treedepth_val)
    )
    tag_file <- gsub("_", "-", tag)
    saveRDS(fit, file.path(run, paste0("fit-", tag_file, ".rds")))
    
    # Hyperparameters
    hp <- rstan::extract(fit, pars = c("mu_r", "sigma_r", "mu_ll", "sigma_ll"))
    hyper <- rbindlist(lapply(names(hp), function(nm) {
      x <- as.numeric(hp[[nm]])
      data.table(parameter = nm, mean = mean(x), sd = sd(x),
                 q025 = quantile(x, 0.025), q975 = quantile(x, 0.975))
    }))
    fwrite(hyper, file.path(run, paste0("hyperparameters-", tag_file, ".csv")))
    
    # Participant draws
    r_draws <- rstan::extract(fit, pars = "r_draw")$r_draw
    list(fit = fit, pid = d$pid, r_draws = r_draws, inconsistent = d$inconsistent,
         tag = tag, treatment = treatment, sample = sample)
  }
  
  # ---- Draw and scored files of one fit ----
  write_draws <- function(res) {
    tag_file <- gsub("_", "-", res$tag)
    saveRDS(list(pid = res$pid, r_draws = res$r_draws),
            file.path(run, paste0("r-draws-", tag_file, ".rds")))
    fwrite(data.table(pid = res$pid, tag = res$tag,
                      r_mean = colMeans(res$r_draws),
                      r_median = apply(res$r_draws, 2, median),
                      inconsistent = res$inconsistent),
           file.path(run, paste0("scored-", tag_file, ".csv")))
  }
  
  diag_of <- function(res) {
    edge <- max(apply(res$r_draws, 2, function(x) mean(x < r_grid[1] + dr / 2 | x > r_grid[R] - dr / 2)))
    rbind(
      stan_diagnostics(res$fit, tag = res$tag, treatment = res$treatment, sample = res$sample,
                       cfg = cfg, max_treedepth = treedepth_val),
      data.table(tag = res$tag, treatment = res$treatment, sample = res$sample,
                 check = "grid_edge_mass_max", value = edge,
                 threshold = edge_mass_max, pass = edge <= edge_mass_max)
    )
  }
  
  # ---- Fits: each treatment and, if requested, its consistent sample ----
  diag_all <- NULL
  for (tr in treatments) {
    res <- fit_one(mpl[treat == tr], tr, tr, "full")
    write_draws(res)
    diag_all <- rbind(diag_all, diag_of(res))
    
    if (consistent_only) {
      res_c <- fit_one(mpl[treat == tr & inconsistent == 0L], paste0(tr, "_consistent"),
                       tr, "consistent")
      write_draws(res_c)
      diag_all <- rbind(diag_all, diag_of(res_c))
    }
  }
  
  # ---- Model property read by the a* step ----
  fwrite(data.table(property = "partial_pooling", value = TRUE), file.path(run, "model-info.csv"))
  
  # ---- Diagnostics (written last: marks the run as complete) ----
  write_diagnostics(run, diag_all)
  
  settings_run <- st
  if (pass$name != "main")
    settings_run <- modifyList(settings_run, list(prior = pass$prior, prior_pass = pass$name))
  model_end(run, cfg, settings = settings_run,
            key = if (pass$name == "main") "mpl" else NULL)
}
