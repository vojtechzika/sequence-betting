# ============================================================
# 03_eu-functions-test.R
#
# PURPOSE
#   Validates the expected-utility functions in 02_eu-functions.R at the
#   design values (endowment, minimum payoff, coin probability and
#   multipliers from 01_design.R) and writes a report stating what each
#   function does and whether it passed its checks.
#
#   Checks
#     - CRRA utility: the logarithmic case at r = 1, the linear case at
#       r = 0, and agreement of the vector and scalar versions
#     - expected utility at a zero stake equals the utility of the endowment
#     - optimal stake: agreement with brute-force maximisation of the
#       expected utility, with the closed-form interior optimum, with the
#       corner solution under risk neutrality, and non-increase in r
#     - certainty equivalent: u(CE) reproduces EU on a grid of risk
#       parameters (including values around r = 1) and stakes, in relative
#       error, and the vector wrapper agrees with the EU-to-CE pipeline
#
# USAGE
#   source() this file, then run eu_functions_test(). It does not need cfg.
#   All checks are run and reported before the function stops on a failure.
#
# INPUT
#   scripts/config/01_design.R       -- design values (design_cfg())
#   scripts/config/02_eu-functions.R -- functions under test
#
# OUTPUT
#   path_rep/eu-functions-test.txt        -- report: functions, checks, result
#   path_der/eu-functions-test-grid.csv   -- certainty-equivalent inversion on
#                                            the grid (internal; EU, CE, u(CE),
#                                            absolute and relative error)
# ============================================================

eu_functions_test <- function() {
  
  if (!exists("design_cfg")) source(here::here("scripts", "config", "01_design.R"))
  if (!exists("crra_u"))     source(here::here("scripts", "config", "02_eu-functions.R"))
  if (!exists("msg"))        msg <- function(...) cat(..., "\n")
  
  dir_rep <- if (exists("path_rep")) path_rep else here::here("data", "output", "reports")
  dir_der <- if (exists("path_der")) path_der else here::here("data", "derived")
  dir.create(dir_rep, showWarnings = FALSE, recursive = TRUE)
  dir.create(dir_der, showWarnings = FALSE, recursive = TRUE)
  
  # ---- Design values ----
  design <- design_cfg()
  e      <- design$seq$endowment
  xmin   <- design$seq$xmin
  p_win  <- design$seq$coin_prob
  mults  <- unlist(design$seq$treatments)
  
  tol     <- 1e-8   # identities that should hold to numerical precision
  tol_rel <- 1e-9   # relative error of the certainty-equivalent inversion
  
  # ---- Result collector ----
  res <- list()
  record <- function(group, check, statistic, value, threshold) {
    res[[length(res) + 1L]] <<- data.frame(
      group = group, check = check, statistic = statistic,
      value = value, threshold = threshold,
      pass = is.finite(value) && value <= threshold,
      stringsAsFactors = FALSE
    )
  }
  
  # ---- 1. CRRA utility ----
  x <- c(1, 2, 10)
  record("utility", "u(x; r = 1) equals log(x)", "max abs error",
         max(abs(crra_u(x, r = 1, xmin = xmin) - log(x))), tol)
  record("utility", "u(x; r = 0) equals x", "max abs error",
         max(abs(crra_u(x, r = 0, xmin = xmin) - x)), tol)
  
  r_vec <- c(0, 0.5, 1)
  x_vec <- c(1, 2, 10)
  u_ref <- vapply(seq_along(r_vec), function(i) crra_u(x_vec[i], r_vec[i], xmin = xmin),
                  numeric(1))
  record("utility", "vector version matches scalar version", "max abs error",
         max(abs(crra_u_vec(x_vec, r_vec, xmin = xmin) - u_ref)), tol)
  
  # ---- 2. Expected utility at a zero stake equals u(e) ----
  r_set <- c(-1, 0, 0.5, 1, 2)
  for (m in mults) {
    err <- max(vapply(r_set, function(r) {
      abs(eu_stake(a = 0, r = r, m = m, e = e, p_win = p_win, xmin = xmin) -
            crra_u(e, r = r, xmin = xmin))
    }, numeric(1)))
    record("expected utility", paste0("EU at a = 0 equals u(e), m = ", m),
           "max abs error", err, tol)
  }
  
  # ---- 3. Optimal stake ----
  r_opt <- c(-1, -0.2, 0, 0.25, 0.5, 1, 2, 5)
  grid  <- 0:e
  a_opt <- function(r, m) a_star_discrete(r = r, m = m, e = e, p_win = p_win, xmin = xmin)
  
  record("optimal stake", "returns a single integer", "failures",
         as.numeric(!(is.integer(a_opt(0, mults[1])) && length(a_opt(0, mults[1])) == 1L)), 0)
  
  for (m in mults) {
    
    # brute force over the integer grid
    n_bad <- sum(vapply(r_opt, function(r) {
      eu_grid <- vapply(grid, function(a) {
        eu_stake(a, r = r, m = m, e = e, p_win = p_win, xmin = xmin)
      }, numeric(1))
      !identical(a_opt(r, m), as.integer(grid[which.max(eu_grid)]))
    }, logical(1)))
    record("optimal stake",
           paste0("matches brute-force maximisation, m = ", m, " (", length(r_opt), " values of r)"),
           "mismatches", n_bad, 0)
    
    # closed-form interior optimum, r > 0 (first-order condition of
    # p u(e + (m - 1)a) + (1 - p) u(e - a) with u(x) = x^(1 - r)/(1 - r)):
    # (W/L)^r = p (m - 1)/(1 - p), W = e + (m - 1)a, L = e - a.
    # The discrete optimum must be one of the two integers bracketing it.
    r_pos <- c(0.25, 0.5, 1, 2, 5)
    n_out <- sum(vapply(r_pos, function(r) {
      k      <- (p_win * (m - 1) / (1 - p_win))^(1 / r)
      a_cont <- if (k > 1) e * (k - 1) / (m - 1 + k) else 0
      !(a_opt(r, m) %in% c(floor(a_cont), ceiling(a_cont)))
    }, logical(1)))
    record("optimal stake",
           paste0("is an integer bracketing the closed-form optimum, m = ", m, " (r > 0)"),
           "violations", n_out, 0)
    
    # risk neutrality: corner solution
    a_rn <- a_opt(0, m)
    a_rn_expected <- if (p_win * m > 1) e else 0
    record("optimal stake",
           paste0("risk-neutral corner solution, m = ", m, " (expected ", a_rn_expected, ")"),
           "abs deviation (ECU)", abs(a_rn - a_rn_expected), 0)
    
    # stake does not increase with risk aversion
    a_path <- vapply(sort(r_opt), function(r) a_opt(r, m), integer(1))
    record("optimal stake",
           paste0("non-increasing in r, m = ", m),
           "violations", sum(diff(a_path) > 0), 0)
  }
  
  # ---- 4. Certainty-equivalent inversion on a grid ----
  r_grid <- c(-2, -1, -0.5, -0.000001, 0, 0.000001, 0.2, 0.5, 0.9, 0.99, 0.999999, 1, 1.000001, 1.01, 2, 5)
  a_grid <- c(0, 1, 5, 10, 25, 50, 75, 99, 100)
  
  out <- do.call(rbind, lapply(mults, function(m) {
    g <- expand.grid(r = r_grid, a = a_grid)
    g$m  <- m
    g$EU <- mapply(function(a, r) eu_stake(a, r, m = m, e = e, p_win = p_win, xmin = xmin),
                   g$a, g$r)
    g$CE      <- ce_from_eu_vec(g$EU, g$r, xmin = xmin)
    g$u_CE    <- mapply(function(c, r) crra_u(c, r, xmin = xmin), g$CE, g$r)
    g$abs_err <- abs(g$u_CE - g$EU)
    g$rel_err <- g$abs_err / pmax(abs(g$EU), .Machine$double.xmin)
    g
  }))
  
  for (m in mults) {
    d <- out[out$m == m, ]
    record("certainty equivalent",
           paste0("CE is finite and positive, m = ", m, " (", nrow(d), " cells)"),
           "failures", sum(!is.finite(d$CE) | d$CE <= 0), 0)
    record("certainty equivalent",
           paste0("u(CE) reproduces EU, m = ", m),
           "max rel error", max(d$rel_err, na.rm = TRUE), tol_rel)
  }
  
  d1 <- out[out$m == mults[1], ]
  ce_wrap <- ce_stake_vec_r(a = d1$a, r = d1$r, m = mults[1], e = e, p_win = p_win, xmin = xmin)
  record("certainty equivalent", "ce_stake_vec_r matches the EU to CE pipeline",
         "max abs error", max(abs(ce_wrap - d1$CE)), tol)
  
  report <- do.call(rbind, res)
  all_pass <- all(report$pass)
  
  # ---- Grid file (internal) ----
  f_grid <- file.path(dir_der, "eu-functions-test-grid.csv")
  data.table::fwrite(out[, c("m", "r", "a", "EU", "CE", "u_CE", "abs_err", "rel_err")], f_grid)
  msg("Saved: ", f_grid)
  
  # ---- Text report ----
  fun_desc <- c(
    crra_u         = "CRRA utility u(x; r) = x^(1 - r)/(1 - r), and ln(x) when |r - 1| < 1e-8; payoffs are floored at the minimum payoff.",
    crra_u_vec     = "Vectorised utility for a vector of payoffs and a vector of risk parameters.",
    eu_stake       = "Expected utility of staking a out of the endowment e on a coin win with probability p: p u(e - a + m a) + (1 - p) u(e - a).",
    eu_stake_vec_r = "Vectorised expected utility for vectors of stakes and risk parameters.",
    a_star_discrete = "Optimal stake: the integer in 0..e maximising expected utility; ties are resolved toward the smaller stake.",
    ce_from_eu_vec  = "Certainty equivalent c solving u(c; r) = EU; closed form when available, root finding otherwise.",
    ce_stake_vec_r  = "Certainty equivalent of a stake: eu_stake_vec_r followed by ce_from_eu_vec."
  )
  
  fmt <- function(v) formatC(v, format = "e", digits = 2)
  rows <- sprintf("  %s  [%s] %s\n        %s = %s (limit %s)",
                  ifelse(report$pass, "PASS", "FAIL"), report$group, report$check,
                  report$statistic, fmt(report$value), fmt(report$threshold))
  
  lines <- c(
    "TEST OF THE EXPECTED-UTILITY FUNCTIONS",
    paste0("Written: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste0("R:       ", R.version.string),
    "",
    paste0("RESULT: ", if (all_pass) "all " else "", sum(report$pass), " of ", nrow(report),
           " checks passed", if (!all_pass) " -- FAILED" else ""),
    "",
    "DESIGN VALUES USED (scripts/config/01_design.R)",
    paste0("  endowment e = ", e, " ECU; minimum payoff = ", xmin,
           "; win probability p = ", p_win, "; multipliers m = ",
           paste(mults, collapse = ", ")),
    "",
    "FUNCTIONS (scripts/config/02_eu-functions.R)",
    paste0("  ", names(fun_desc), ": ", fun_desc),
    "",
    "CHECKS",
    rows,
    "",
    paste0("The grid of certainty-equivalent inversions (", nrow(out),
           " cells: r from -2 to 5, including values around 0 and 1, stakes from 0 to ", e,
           " ECU) is saved as derived/eu-functions-test-grid.csv.")
  )
  
  f_rep <- file.path(dir_rep, "eu-functions-test.txt")
  writeLines(lines, f_rep)
  msg("Saved: ", f_rep)
  
  if (!all_pass) {
    stop("EU function tests failed: ", sum(!report$pass), " of ", nrow(report),
         " checks. See ", f_rep, call. = FALSE)
  }
  invisible(report)
}
