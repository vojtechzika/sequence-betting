# ============================================================
# test_eu_functions.R
# Minimal validation for 03_eu_functions.R
# ============================================================

source("../analysis/scripts/config/03_eu_functions.R")

tol <- 1e-8

check <- function(cond, msg) {
  if (!isTRUE(cond)) stop(msg, call. = FALSE)
}

cat("Running EU helper tests...\n")

# ------------------------------------------------------------
# 1. Basic CRRA utility sanity checks
# ------------------------------------------------------------

x <- c(1, 2, 10)

u_log <- crra_u(x, r = 1, xmin = 0.01)
check(all.equal(u_log, log(x), tolerance = tol) == TRUE,
      "crra_u failed at r = 1")

u_r0 <- crra_u(x, r = 0, xmin = 0.01)
check(all.equal(u_r0, x, tolerance = tol) == TRUE,
      "crra_u failed at r = 0")

# vector version matches scalar version
r_vec <- c(0, 0.5, 1)
x_vec <- c(1, 2, 10)
u_vec <- crra_u_vec(x_vec, r_vec, xmin = 0.01)
u_ref <- vapply(seq_along(r_vec), function(i) {
  crra_u(x_vec[i], r_vec[i], xmin = 0.01)
}, numeric(1))
check(all.equal(u_vec, u_ref, tolerance = tol) == TRUE,
      "crra_u_vec does not match scalar crra_u")

cat("1 OK: utility functions\n")

# ------------------------------------------------------------
# 2. EU sanity checks
# ------------------------------------------------------------

# At a = 0, EU should equal u(e)
e <- 100
m <- 2.5
p <- 0.5
xmin <- 0.01

for (r in c(-1, 0, 0.5, 1, 2)) {
  eu0 <- eu_stake(a = 0, r = r, m = m, e = e, p_win = p, xmin = xmin)
  u0  <- crra_u(e, r = r, xmin = xmin)
  check(all.equal(eu0, u0, tolerance = tol) == TRUE,
        paste("eu_stake failed at a=0 for r =", r))
}

cat("2 OK: EU baseline identity\n")

# ------------------------------------------------------------
# 3. a_star_discrete sanity checks
# ------------------------------------------------------------

# return type
a1 <- a_star_discrete(r = 0, m = 2.5, e = 100, p_win = 0.5, xmin = 0.01)
check(is.integer(a1) && length(a1) == 1L,
      "a_star_discrete does not return integer(1)")

# compare against brute-force maximization
grid <- 0:100
for (r in c(-1, -0.2, 0, 0.5, 1, 2, 5)) {
  EU_grid <- vapply(grid, function(a) {
    eu_stake(a, r = r, m = 2.5, e = 100, p_win = 0.5, xmin = 0.01)
  }, numeric(1))
  
  a_ref <- as.integer(grid[which.max(EU_grid)])
  a_fun <- a_star_discrete(r = r, m = 2.5, e = 100, p_win = 0.5, xmin = 0.01)
  
  check(identical(a_fun, a_ref),
        paste("a_star_discrete mismatch for r =", r))
}

cat("3 OK: a_star_discrete\n")

# ------------------------------------------------------------
# 4. CE inversion check: u(CE) == EU
# ------------------------------------------------------------

r_test <- c(-1, -0.2, 0, 0.5, 0.9, 0.999999, 1, 1.000001, 2, 5)
a_test <- c(0, 1, 10, 50, 99, 100, 30, 70, 5, 80)

EU_test <- eu_stake_vec_r(
  a = a_test,
  r = r_test,
  m = 2.5,
  e = 100,
  p_win = 0.5,
  xmin = 0.01
)

CE_test <- ce_from_eu_vec(EU_test, r_test, xmin = 0.01)

# CE should be finite and positive
check(all(is.finite(CE_test)), "ce_from_eu_vec returned non-finite values")
check(all(CE_test > 0), "ce_from_eu_vec returned non-positive CE")

# plugging CE back into utility should reproduce EU
EU_back <- vapply(seq_along(r_test), function(i) {
  crra_u(CE_test[i], r_test[i], xmin = 0.01)
}, numeric(1))

err <- max(abs(EU_back - EU_test))
print(data.frame(r = r_test, a = a_test, EU = EU_test, CE = CE_test, EU_back = EU_back,
                 abs_err = abs(EU_back - EU_test)))
check(err < 1e-6, paste("CE inversion failed; max abs error =", err))

cat("4 OK: CE inversion\n")

# ------------------------------------------------------------
# 5. Vector CE wrapper check
# ------------------------------------------------------------

CE_wrap <- ce_stake_vec_r(
  a = a_test,
  r = r_test,
  m = 2.5,
  e = 100,
  p_win = 0.5,
  xmin = 0.01
)

check(all.equal(CE_wrap, CE_test, tolerance = tol) == TRUE,
      "ce_stake_vec_r does not match EU -> CE pipeline")

cat("5 OK: ce_stake_vec_r\n")

# ------------------------------------------------------------
# 6. Stress grid
# ------------------------------------------------------------

r_grid <- c(-2, -1, -0.5, 0, 0.2, 0.5, 0.9, 0.99, 0.999999, 1, 1.000001, 1.01, 2, 5)
a_grid <- c(0, 1, 5, 10, 25, 50, 75, 99, 100)

stress <- expand.grid(r = r_grid, a = a_grid)

stress$EU <- mapply(function(a, r) {
  eu_stake(a, r, m = 2.5, e = 100, p_win = 0.5, xmin = 0.01)
}, stress$a, stress$r)

stress$CE <- ce_from_eu_vec(stress$EU, stress$r, xmin = 0.01)
stress$EU_back <- mapply(function(c, r) {
  crra_u(c, r, xmin = 0.01)
}, stress$CE, stress$r)

stress$abs_err <- abs(stress$EU_back - stress$EU)

print(summary(stress$abs_err))

bad <- subset(stress, !is.finite(CE) | abs_err > 1e-5)
if (nrow(bad) > 0) {
  print(bad)
  stop("Stress test failed")
}

cat("6 OK: stress grid\n")
cat("All EU helper tests passed.\n")