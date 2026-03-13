# ============================================================
# 03_eu_functions.R
# Pure expected utility helpers (no config)
# ============================================================

R_TOL_DEFAULT <- 1e-8

crra_u <- function(x, r, xmin, r_tol = R_TOL_DEFAULT) {
  x <- pmax(x, xmin)
  if (abs(r - 1) < r_tol) {
    log(x)
  } else {
    x^(1 - r) / (1 - r)
  }
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
  
  # keep prereg tie-break toward smaller stake
  grid <- sort(unique(as.integer(grid)))
  
  EU <- vapply(
    grid,
    function(a) eu_stake(a, r, m, e, p_win, xmin, r_tol),
    numeric(1)
  )
  
  as.integer(grid[which.max(EU)])
}

# ----------------------------
# Vectorized CE helpers (handle vector r)
# ----------------------------

crra_u_vec <- function(x, r, xmin, r_tol = R_TOL_DEFAULT) {
  stopifnot(length(x) == length(r))
  
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
  
  p_win * crra_u_vec(w_win, r, xmin = xmin, r_tol = r_tol) +
    (1 - p_win) * crra_u_vec(w_lose, r, xmin = xmin, r_tol = r_tol)
}

ce_from_eu_vec <- function(EU, r, xmin, r_tol = R_TOL_DEFAULT) {
  stopifnot(length(EU) == length(r))
  
  n <- length(r)
  ce <- rep(NA_real_, n)
  
  for (i in seq_len(n)) {
    
    EU_i <- EU[i]
    r_i  <- r[i]
    
    if (!is.finite(EU_i) || !is.finite(r_i)) {
      ce[i] <- NA_real_
      next
    }
    
    # log case: exact inverse
    if (abs(r_i - 1) < r_tol) {
      ce[i] <- exp(EU_i)
      ce[i] <- pmax(ce[i], xmin)
      next
    }
    
    # try closed form first when numerically admissible
    val <- (1 - r_i) * EU_i
    
    use_closed_form <- is.finite(val) && (val > 0)
    
    if (use_closed_form) {
      ce_cf <- val^(1 / (1 - r_i))
      
      if (is.finite(ce_cf) && ce_cf > 0) {
        ce[i] <- pmax(ce_cf, xmin)
        next
      }
    }
    
    # numerical inversion of u(c; r) = EU, c > 0
    f_root <- function(c) {
      crra_u(c, r_i, xmin = xmin, r_tol = r_tol) - EU_i
    }
    
    lower <- xmin
    f_lower <- f_root(lower)
    
    if (!is.finite(f_lower)) {
      ce[i] <- NA_real_
      next
    }
    
    # expand upper bound until the root is bracketed
    upper <- max(1, xmin * 2)
    f_upper <- f_root(upper)
    iter <- 0L
    max_iter <- 200L
    
    while (is.finite(f_upper) && (f_lower * f_upper > 0) && iter < max_iter) {
      upper <- upper * 2
      f_upper <- f_root(upper)
      iter <- iter + 1L
    }
    
    if (!is.finite(f_upper) || (f_lower * f_upper > 0)) {
      ce[i] <- NA_real_
      next
    }
    
    root <- tryCatch(
      uniroot(f_root, lower = lower, upper = upper, tol = 1e-10)$root,
      error = function(e) NA_real_
    )
    
    ce[i] <- if (is.finite(root)) pmax(root, xmin) else NA_real_
  }
  
  ce
}

ce_stake_vec_r <- function(a, r, m, e, p_win, xmin, r_tol = R_TOL_DEFAULT) {
  EU <- eu_stake_vec_r(
    a = a,
    r = r,
    m = m,
    e = e,
    p_win = p_win,
    xmin = xmin,
    r_tol = r_tol
  )
  
  ce_from_eu_vec(EU = EU, r = r, xmin = xmin, r_tol = r_tol)
}