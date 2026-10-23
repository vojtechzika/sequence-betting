# ============================================================
# 01_mpl-stan-bernoulli.R
#
# PURPOSE
#   Estimates individual CRRA risk parameters r_i from Holt-Laury MPL
#   choices with a Bayesian Bernoulli (logit) choice model: the
#   probability of the safe choice is logistic in lambda_i times the
#   expected-utility difference between lotteries A and B, with a
#   normal prior on r_i and a lognormal prior on lambda_i
#   (cfg$model$stan$mpl_bernoulli$prior). The model is fitted to the pooled
#   sample (both treatments combined) and to each treatment. Because
#   participants are not pooled, estimates for the Holt-Laury-consistent
#   subset are obtained by filtering these fits (07_a-star-draws.R). Convergence
#   diagnostics are part of this script.
#
#   Each participant's posterior depends only on that participant's own
#   choices and has two parameters, so it can also be computed exactly on a
#   grid. The posterior mean of r_i from the grid is compared with the
#   sampler's (check grid_max_abs_diff in diagnostics.csv). The check
#   matters for participants who choose the risky lottery in every row: their
#   posterior has a narrow mode far in the lower tail of the prior and a flat
#   spurious mode at high r and low lambda, and a chain captured by the latter
#   is not detected by the usual convergence diagnostics. The sampler uses
#   Stan's default random initial values.
#
# DEPARTURE FROM THE PREREGISTERED ESTIMATION PROCEDURE
#   The preregistered model (likelihood, priors, estimand) is estimated in
#   Stan for every participant, with one exception that is handled by
#   exact computation. For a participant who chooses the risky lottery in
#   every row, the posterior of r_i is the prior truncated at a hard
#   boundary: the observed choice in the first row is probable only if r_i
#   lies below the value at which the first lottery pair is indifferent
#   (r = -2.48 under the design), and the log-likelihood falls by more than
#   100 within 0.001 of r_i above that value. The posterior mass lies in a
#   narrow region far in the lower tail of the prior. The sampler cannot
#   traverse this boundary: with initial values at the posterior mode it
#   saturated the maximum tree depth and produced divergent transitions, and
#   with random initial values all chains end in a distant local mode at
#   high r and low lambda that has negligible posterior mass, where the
#   convergence diagnostics are satisfied but the estimate of r_i has the
#   wrong sign.
#
#   Each participant's posterior has only two parameters (r_i, log lambda_i),
#   so it can be computed to numerical precision on a grid. The rule is:
#   for every participant whose sampler posterior mean of r_i differs from
#   the grid posterior mean by more than grid_tol, the draws of r_i are
#   replaced by an equal number of draws from the exact posterior, evaluated
#   on a grid that is refined around the region carrying the posterior mass.
#   The model is unchanged, no participant is excluded, and the rule uses
#   the Holt-Laury choices only. A replacement share above exact_share_max
#   fails a check, so that this step cannot mask a general sampling problem.
#   Replaced participants are listed in exact-draw-participants.csv. The
#   saved Stan fits (fit-<tag>.rds) contain the original sampler draws.
#   Participants for whom the posterior has this boundary are those with an
#   all-risky choice pattern; the data then identify that r_i lies below the
#   boundary, and the location of the posterior within that region reflects
#   the prior.
#
# INPUT
#   path_src/holt-laury.csv               -- MPL choices c1..cK
#   path_src/participants.csv             -- treatment and sex
#   path_der/switching-point-r.csv        -- inconsistency flag per pid
#   scripts/stan/mpl-bernoulli.stan
#
# OUTPUT (model folder path_mod/mpl-bernoulli/)
#   fit-<tag>.rds        -- Stan fit objects
#   r-draws-<tag>.rds    -- list(pid, r_draws): posterior draws of r_i
#                           (exact-posterior draws for replaced participants)
#   exact-draw-participants.csv -- participants whose draws were replaced:
#                           sampler, grid and exact posterior means of r_i
#   scored-<tag>.csv     -- r_mean, r_median, inconsistent per pid
#   summary.csv          -- summary of r_mean by tag, consistency, sex
#   diagnostics.csv      -- convergence checks per fit, and the largest
#                           absolute difference between the sampler's and the
#                           grid posterior mean of r_i (grid_max_abs_diff)
#   model-info.csv       -- model property read by the a* step
#                           (partial_pooling = FALSE: no pooling of r_i across
#                           participants, so subsets of participants can be
#                           obtained by filtering the full-sample fit)
#   inputs.csv, run-info.txt
#
# TAGS
#   pooled, <treatment>                        (all participants; "pooled" here
#                                               means both treatments combined)
#
# NOTES
#   - The model folder is built as mpl-bernoulli.tmp and renamed when
#     complete; an existing folder is never overwritten
#   - Convergence criteria are cfg$model$convergence; the tolerance of the
#     grid comparison is grid_tol below (in units of r)
#   - Unrecognised MPL choice labels or unknown treatments stop the script
# ============================================================

options(mc.cores = parallel::detectCores())
rstan::rstan_options(auto_write = FALSE)

# ------------------------------------------------------------
# Exact posterior of (r_i, log lambda_i) on a grid, per participant
#   ymat: participants x rows matrix of choices (1 = safe, 0 = risky)
#   p:    probability of the high payoff in each row
# Returns a matrix with the posterior mean of r_i and the location of the
# posterior mode (r, log lambda). The prior on lambda is lognormal, so on
# the log scale it is normal and no Jacobian term enters.
# ------------------------------------------------------------
mpl_grid_posterior <- function(ymat, p, A1, A2, B1, B2, prior,
                               r_grid  = seq(-8, 4, by = 0.025),
                               ll_grid = seq(-2, 8, by = 0.05)) {
  crra <- function(x, r) ifelse(abs(r - 1) < 1e-8, log(x), x^(1 - r) / (1 - r))
  d1  <- crra(A1, r_grid) - crra(B1, r_grid)   # utility difference, high payoffs
  d2  <- crra(A2, r_grid) - crra(B2, r_grid)   # utility difference, low payoffs
  lam <- exp(ll_grid)
  
  lprior <- outer(dnorm(r_grid,  prior$r_mean, prior$r_sd, log = TRUE),
                  dnorm(ll_grid, log(prior$lambda_mean), prior$lambda_sd, log = TRUE),
                  "+")
  
  # Log-probability of a safe (S) and of a risky (R) choice in each row
  K <- length(p)
  logS <- logR <- vector("list", K)
  for (k in seq_len(K)) {
    z <- outer(p[k] * d1 + (1 - p[k]) * d2, lam)
    logS[[k]] <- -(pmax(-z, 0) + log1p(exp(-abs(z))))   # log logistic(z)
    logR[[k]] <- logS[[k]] - z                          # log (1 - logistic(z))
  }
  
  t(vapply(seq_len(nrow(ymat)), function(i) {
    lp <- lprior
    for (k in seq_len(K)) lp <- lp + (if (ymat[i, k] == 1L) logS[[k]] else logR[[k]])
    w <- exp(lp - max(lp)); w <- w / sum(w)
    j <- arrayInd(which.max(lp), dim(lp))
    c(r_mean = sum(w * r_grid), r_mode = r_grid[j[1]], ll_mode = ll_grid[j[2]])
  }, numeric(3)))
}

# ------------------------------------------------------------
# Draws of r_i from the exact posterior of one participant
#   y:    choices of the participant (1 = safe, 0 = risky), one per row
#   n:    number of draws
# The posterior of (r_i, log lambda_i) is evaluated on a grid. A coarse pass
# locates the region carrying the posterior mass; the grid is then refined
# there (step r_step / refine, at most max_points points) and r_i is drawn
# from its marginal posterior, with a uniform draw within the grid cell.
# Returns the draws and the posterior mean of r_i on the refined grid.
# ------------------------------------------------------------
mpl_exact_draws <- function(y, p, A1, A2, B1, B2, prior, n, seed,
                            r_step = 0.025, refine = 250, max_points = 40000L,
                            mass_tol = 1e-10,
                            ll_grid = seq(-2, 8, by = 0.05)) {
  crra <- function(x, r) ifelse(abs(r - 1) < 1e-8, log(x), x^(1 - r) / (1 - r))
  
  # Log posterior (up to a constant): r in rows, log lambda in columns
  log_post <- function(r_grid) {
    d1  <- crra(A1, r_grid) - crra(B1, r_grid)
    d2  <- crra(A2, r_grid) - crra(B2, r_grid)
    lam <- exp(ll_grid)
    lp  <- outer(dnorm(r_grid,  prior$r_mean, prior$r_sd, log = TRUE),
                 dnorm(ll_grid, log(prior$lambda_mean), prior$lambda_sd, log = TRUE),
                 "+")
    for (k in seq_along(p)) {
      z    <- outer(p[k] * d1 + (1 - p[k]) * d2, lam)
      logS <- -(pmax(-z, 0) + log1p(exp(-abs(z))))   # log logistic(z)
      lp   <- lp + (if (y[k] == 1L) logS else logS - z)
    }
    lp
  }
  
  # Marginal posterior weights of r on a grid
  r_weights <- function(r_grid) {
    lp <- log_post(r_grid)
    w  <- rowSums(exp(lp - max(lp)))
    w / sum(w)
  }
  
  # Coarse pass: region carrying the posterior mass
  r_coarse <- seq(-8, 4, by = r_step)
  w0       <- r_weights(r_coarse)
  keep     <- which(w0 / max(w0) > mass_tol)
  lo       <- r_coarse[min(keep)] - r_step
  hi       <- r_coarse[max(keep)] + r_step
  
  # Refined pass
  step   <- max(r_step / refine, (hi - lo) / max_points)
  r_fine <- seq(lo, hi, by = step)
  w      <- r_weights(r_fine)
  
  set.seed(seed)
  idx <- sample.int(length(r_fine), size = n, replace = TRUE, prob = w)
  list(draws = r_fine[idx] + runif(n, -step / 2, step / 2),
       mean  = sum(w * r_fine))
}

mpl_stan_bernoulli <- function(cfg) {
  
  run <- model_begin("mpl-bernoulli")
  if (is.null(run)) return(invisible(NULL))
  
  seed            <- as.integer(cfg$run$seed)
  grid_tol        <- 0.10   # tolerated |sampler - grid| posterior mean of r_i
  exact_share_max <- 0.05   # largest share of participants with replaced draws
  design          <- cfg$design
  treatments      <- names(design$seq$treatments)
  
  label_safe  <- toupper(trimws(as.character(design$mpl$label_safe)))
  label_risky <- toupper(trimws(as.character(design$mpl$label_risky)))
  stopifnot(
    length(label_safe) == 1L, length(label_risky) == 1L,
    nzchar(label_safe), nzchar(label_risky),
    label_safe != label_risky
  )
  
  # ---- Design of the MPL ----
  K_cfg  <- as.integer(design$mpl$K)
  A_high <- as.numeric(design$mpl$A_high)
  A_low  <- as.numeric(design$mpl$A_low)
  B_high <- as.numeric(design$mpl$B_high)
  B_low  <- as.numeric(design$mpl$B_low)
  stopifnot(K_cfg > 0L, is.finite(A_high), is.finite(A_low),
            is.finite(B_high), is.finite(B_low),
            A_high > 0, A_low > 0, B_high > 0, B_low > 0)
  
  p_sched <- as.numeric(design$mpl$p)
  stopifnot(length(p_sched) == K_cfg, all(is.finite(p_sched)),
            all(p_sched > 0), all(p_sched < 1), !is.unsorted(p_sched))
  
  # ---- Sampler settings and priors ----
  st                <- cfg$model$stan$mpl_bernoulli
  # This script fits the individual-level model only (fixed priors, no pooling)
  if (!identical(st$partial_pooling, FALSE)) {
    stop("mpl_bernoulli$partial_pooling must be FALSE: this script fits the individual-level ",
         "model. A pooled model needs its own script.", call. = FALSE)
  }
  iter_val          <- as.integer(st$iter)
  warmup_val        <- as.integer(st$warmup)
  chains_val        <- as.integer(st$chains)
  adapt_delta_val   <- as.numeric(st$adapt_delta)
  treedepth_val     <- as.integer(st$treedepth)
  prior_r_mean      <- as.numeric(st$prior$r_mean)
  prior_r_sd        <- as.numeric(st$prior$r_sd)
  prior_lambda_mean <- as.numeric(st$prior$lambda_mean)
  prior_lambda_sd   <- as.numeric(st$prior$lambda_sd)
  
  stopifnot(iter_val > 0L, warmup_val >= 0L, chains_val > 0L,
            is.finite(adapt_delta_val), adapt_delta_val > 0, adapt_delta_val < 1,
            treedepth_val > 0L,
            length(prior_r_mean) == 1L, is.finite(prior_r_mean),
            length(prior_r_sd) == 1L, is.finite(prior_r_sd), prior_r_sd > 0,
            length(prior_lambda_mean) == 1L, is.finite(prior_lambda_mean),
            prior_lambda_mean > 0,
            length(prior_lambda_sd) == 1L, is.finite(prior_lambda_sd),
            prior_lambda_sd > 0)
  
  # ---- Inputs ----
  stan_file <- here::here("scripts", "stan", "mpl-bernoulli.stan")
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
  
  # Treatment and sex from participants.csv (one row per pid)
  if (anyDuplicated(par$pid)) stop("Repeated pid in participants.csv.")
  if (!all(mpl$pid %in% par$pid)) stop("MPL participants missing from participants.csv.")
  par[, sex_label := fcase(toupper(trimws(sex)) == "F", "Women",
                           toupper(trimws(sex)) == "M", "Men",
                           default = NA_character_)]
  mpl <- merge(mpl, par[, .(pid, treat = as.character(treat), sex_label)], by = "pid")
  if (!all(mpl$treat %in% treatments)) stop("Unknown treatment in MPL participants.")
  
  # Holt-Laury inconsistency flag, as defined in 22_score-switching-point-r.R
  if (!all(mpl$pid %in% sw$pid)) stop("MPL participants missing from switching-point-r.csv.")
  mpl <- merge(mpl, sw[, .(pid, inconsistent = as.integer(inconsistent))], by = "pid")
  
  choice_cols <- grep("^c\\d+$", names(mpl), value = TRUE)
  choice_cols <- choice_cols[order(as.integer(sub("^c", "", choice_cols)))]
  stopifnot(length(choice_cols) == K_cfg)
  
  sm <- rstan::stan_model(stan_file)
  
  # ---- Long format with coded choices (1 = safe, 0 = risky) ----
  melt_and_code <- function(d) {
    long <- melt(d, id.vars = c("pid", "treat"), measure.vars = choice_cols,
                 variable.name = "row", value.name = "choice")
    long[, row    := as.integer(sub("^c", "", row))]
    long[, choice := toupper(trimws(as.character(choice)))]
    if (!all(long$choice %in% c(label_safe, label_risky))) {
      stop("MPL choices contain labels other than the configured safe/risky labels.")
    }
    long[, y := fifelse(choice == label_safe, 1L, 0L)]
    long
  }
  
  # ---- Fit one subset ----
  fit_one <- function(d, tag, treatment, sample) {
    
    if (nrow(d) == 0L) stop("No MPL participants for tag '", tag, "'.")
    tag_file <- gsub("_", "-", tag)    # tags are written with hyphens in file names
    
    long       <- melt_and_code(d)
    pid_levels <- sort(unique(long$pid))
    long[, uid := match(pid, pid_levels)]
    Tobs       <- nrow(long)
    
    data_list <- list(
      N = length(pid_levels), T = Tobs, uid = long$uid,
      p  = p_sched[long$row],
      A1 = A_high, A2 = A_low, B1 = B_high, B2 = B_low,
      y  = long$y,
      prior_r_mean = prior_r_mean, prior_r_sd = prior_r_sd,
      prior_lambda_mean = prior_lambda_mean, prior_lambda_sd = prior_lambda_sd
    )
    
    # Exact grid posterior per participant, used as a reference for the sampler
    ymat <- dcast(long[, .(pid, row, y)], pid ~ row, value.var = "y")
    ymat <- as.matrix(ymat[match(pid_levels, pid), -"pid"])
    grid <- mpl_grid_posterior(ymat, p_sched, A_high, A_low, B_high, B_low,
                               prior = list(r_mean = prior_r_mean, r_sd = prior_r_sd,
                                            lambda_mean = prior_lambda_mean,
                                            lambda_sd = prior_lambda_sd))
    fit <- rstan::sampling(
      sm, data = data_list,
      iter = iter_val, warmup = warmup_val, chains = chains_val,
      seed = seed,
      control = list(adapt_delta = adapt_delta_val, max_treedepth = treedepth_val)
    )
    
    saveRDS(fit, file.path(run, paste0("fit-", tag_file, ".rds")))
    
    post <- rstan::extract(fit)
    
    # Participants whose sampler posterior mean differs from the exact grid
    # posterior mean get draws from the exact posterior (see header)
    sampler_mean <- colMeans(post$r)
    off          <- which(abs(sampler_mean - grid[, "r_mean"]) > grid_tol)
    exact_log    <- data.table(pid = character(), tag = character(),
                               sampler_mean = numeric(), grid_mean = numeric(),
                               exact_mean = numeric())
    for (j in off) {
      ex <- mpl_exact_draws(ymat[j, ], p_sched, A_high, A_low, B_high, B_low,
                            prior = list(r_mean = prior_r_mean, r_sd = prior_r_sd,
                                         lambda_mean = prior_lambda_mean,
                                         lambda_sd = prior_lambda_sd),
                            n = nrow(post$r), seed = seed + j)
      post$r[, j] <- ex$draws
      exact_log <- rbind(exact_log,
                         data.table(pid = pid_levels[j], tag = tag,
                                    sampler_mean = sampler_mean[j],
                                    grid_mean = grid[j, "r_mean"],
                                    exact_mean = ex$mean))
    }
    
    saveRDS(list(pid = pid_levels, r_draws = post$r),
            file.path(run, paste0("r-draws-", tag_file, ".rds")))
    
    scored <- data.table(
      pid          = pid_levels,
      tag          = tag,
      r_mean       = apply(post$r, 2, mean),
      r_median     = apply(post$r, 2, median)
    )
    scored <- merge(scored, mpl[, .(pid, inconsistent, sex_label)], by = "pid")
    fwrite(scored[, .(pid, tag, r_mean, r_median, inconsistent)],
           file.path(run, paste0("scored-", tag_file, ".csv")))
    
    # Agreement between the draws used downstream and the grid posterior mean
    # of r_i, and the share of participants whose draws were replaced
    grid_diff <- max(abs(colMeans(post$r) - grid[, "r_mean"]))
    grid_row  <- data.table(tag = tag, treatment = treatment, sample = sample,
                            check = "grid_max_abs_diff", value = grid_diff,
                            threshold = grid_tol, pass = grid_diff <= grid_tol)
    exact_share <- length(off) / length(pid_levels)
    exact_row   <- data.table(tag = tag, treatment = treatment, sample = sample,
                              check = "exact_draw_share", value = exact_share,
                              threshold = exact_share_max,
                              pass = exact_share <= exact_share_max)
    
    list(
      scored = scored,
      exact  = exact_log,
      diag   = rbind(stan_diagnostics(fit, tag = tag, treatment = treatment,
                                      sample = sample, cfg = cfg,
                                      max_treedepth = treedepth_val),
                     grid_row, exact_row)
    )
  }
  
  # ---- Fits: pooled and by treatment ----
  results <- list()
  results[["pooled"]] <- fit_one(mpl, "pooled", "pooled", "full")
  for (tr in treatments) {
    results[[tr]] <- fit_one(mpl[treat == tr], tr, tr, "full")
  }
  
  # ---- Summary of r_mean by tag, consistency and sex ----
  summarise_scored <- function(d, label) {
    if (nrow(d) == 0L) return(NULL)
    data.table(
      group = label, n = nrow(d),
      r_mean = mean(d$r_mean), r_sd = sd(d$r_mean), r_median = median(d$r_mean),
      r_q25 = quantile(d$r_mean, 0.25), r_q75 = quantile(d$r_mean, 0.75),
      r_min = min(d$r_mean), r_max = max(d$r_mean),
      pct_inconsistent = mean(d$inconsistent == 1L)
    )
  }
  
  summary_list <- list()
  for (tag in names(results)) {
    d <- results[[tag]]$scored
    for (sub in c("all", "consistent", "inconsistent")) {
      d_sub <- switch(sub, all = d,
                      consistent   = d[inconsistent == 0L],
                      inconsistent = d[inconsistent == 1L])
      grp <- paste0(tag, " / ", sub)
      summary_list[[grp]] <- summarise_scored(d_sub, grp)
      for (sx in c("Women", "Men")) {
        grp_sx <- paste0(grp, " / ", sx)
        summary_list[[grp_sx]] <- summarise_scored(d_sub[sex_label == sx], grp_sx)
      }
    }
  }
  fwrite(rbindlist(summary_list), file.path(run, "summary.csv"))
  
  # ---- Participants whose draws come from the exact posterior ----
  fwrite(rbindlist(lapply(results, `[[`, "exact")),
         file.path(run, "exact-draw-participants.csv"))
  
  # ---- Diagnostics (written last: marks the run as complete) ----
  # Model property read by the a* step: whether participants are partially pooled
  fwrite(data.table(property = "partial_pooling", value = FALSE), file.path(run, "model-info.csv"))
  
  write_diagnostics(run, rbindlist(lapply(results, `[[`, "diag")))
  
  model_end(run, cfg, settings = st, key = "mpl")
}
