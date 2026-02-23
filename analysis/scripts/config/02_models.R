model_cfg <- function() {
  list(
    stan = list(
      mpl = list(
        pilot = list(iter=1000, warmup=500,  chains=2, adapt_delta=0.90,  treedepth=10),
        main  = list(iter=4000, warmup=2000, chains=4, adapt_delta=0.995, treedepth=15)
      ),
      rq1 = list(
        pilot = list(iter=1500, warmup=750, chains=2, adapt_delta=0.9, treedepth=12),
        main  = list(iter=4000, warmup=2000, chains=4, adapt_delta=0.99, treedepth=15)
      ),
      rq2 = list(
        pilot = list(iter=1500, warmup=750, chains=2, adapt_delta=0.90, treedepth=12),
        main  = list(iter=4000, warmup=2000, chains=4, adapt_delta=0.99, treedepth=15)
      )
    )
  )
}

R_TOL_DEFAULT <- 1e-8

# ----------------------------
# EU helpers (pure functions)
# ----------------------------

crra_u <- function(x, r, xmin, r_tol = R_TOL_DEFAULT) {
  x <- pmax(x, xmin)
  
  if (abs(r - 1) < r_tol) {
    log(x)
  } else {
    x^(1 - r) / (1 - r)
  }
}

# Expected utility for stake a (vectorized in a)
eu_stake <- function(a, r, m, e, p_win, xmin, r_tol = R_TOL_DEFAULT) {
  stopifnot(is.finite(p_win), p_win >= 0, p_win <= 1)
  stopifnot(is.finite(e), e > 0)
  stopifnot(is.finite(xmin), xmin > 0)
  
  w_win  <- (e - a) + m * a
  w_lose <- (e - a)
  
  p_win * crra_u(w_win,  r, xmin = xmin, r_tol = r_tol) +
    (1 - p_win) * crra_u(w_lose, r, xmin = xmin, r_tol = r_tol)
}

# Discrete EU maximization: a in {0,...,e}, ties -> smaller a
a_star_discrete <- function(r, m, e, p_win, xmin, r_tol = R_TOL_DEFAULT, grid = NULL) {
  stopifnot(is.finite(m), m > 1)
  
  if (is.null(grid)) grid <- 0:e
  
  EU <- eu_stake(grid, r, m, e = e, p_win = p_win, xmin = xmin, r_tol = r_tol)
  grid[which.max(EU)]  # tie -> smaller stake (first max)
}