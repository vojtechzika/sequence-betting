# ============================================================
# 03_eu_functions.R
# Pure expected utility helpers (no config)
# ============================================================

R_TOL_DEFAULT <- 1e-8

crra_u <- function(x, r, xmin, r_tol = R_TOL_DEFAULT) {
  x <- pmax(x, xmin)
  if (abs(r - 1) < r_tol) log(x) else x^(1 - r) / (1 - r)
}

eu_stake <- function(a, r, m, e, p_win, xmin, r_tol = R_TOL_DEFAULT) {
  w_win  <- (e - a) + m * a
  w_lose <- (e - a)
  p_win * crra_u(w_win, r, xmin, r_tol) +
    (1 - p_win) * crra_u(w_lose, r, xmin, r_tol)
}

a_star_discrete <- function(r, m, e, p_win, xmin,
                            r_tol = R_TOL_DEFAULT, grid = NULL) {
  if (is.null(grid)) grid <- 0:e
  EU <- eu_stake(grid, r, m, e, p_win, xmin, r_tol)
  grid[which.max(EU)]
}

# ----------------------------
# Vectorized CE helpers (handle vector r)
# ----------------------------

crra_u_vec <- function(x, r, xmin, r_tol = R_TOL_DEFAULT) {
  x <- pmax(x, xmin)
  is_log <- abs(r - 1) < r_tol
  out <- numeric(length(r))
  out[is_log]  <- log(x[is_log])
  out[!is_log] <- x[!is_log]^(1 - r[!is_log]) / (1 - r[!is_log])
  out
}

eu_stake_vec_r <- function(a, r, m, e, p_win, xmin, r_tol = R_TOL_DEFAULT) {
  stopifnot(length(a) == length(r))
  w_win  <- (e - a) + m * a
  w_lose <- (e - a)
  p_win * crra_u_vec(w_win,  r, xmin = xmin, r_tol = r_tol) +
    (1 - p_win) * crra_u_vec(w_lose, r, xmin = xmin, r_tol = r_tol)
}

ce_from_eu_vec <- function(EU, r, xmin, r_tol = R_TOL_DEFAULT) {
  stopifnot(length(EU) == length(r))
  is_log <- abs(r - 1) < r_tol
  
  ce <- numeric(length(r))
  ce[is_log] <- exp(EU[is_log])
  
  val <- (1 - r[!is_log]) * EU[!is_log]
  ok  <- is.finite(val) & val > 0
  ce_tmp <- rep(xmin, length(val))
  ce_tmp[ok] <- val[ok]^(1 / (1 - r[!is_log][ok]))
  ce[!is_log] <- ce_tmp
  
  pmax(ce, xmin)
}

ce_stake_vec_r <- function(a, r, m, e, p_win, xmin, r_tol = R_TOL_DEFAULT) {
  EU <- eu_stake_vec_r(a = a, r = r, m = m, e = e, p_win = p_win, xmin = xmin, r_tol = r_tol)
  ce_from_eu_vec(EU = EU, r = r, xmin = xmin, r_tol = r_tol)
}