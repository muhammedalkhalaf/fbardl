## Reference values from the Stata modules fbardl 1.3.0 and ardl
## (ardlbounds, Kripfganz and Schneider 2020) on the data below.
.ts_data <- function() {
  set.seed(20260928)
  for (i in 1:10) {
    x1 <- cumsum(rnorm(50)); x2 <- cumsum(rnorm(50))
    u <- as.numeric(arima.sim(list(ar = 0.5), 50))
  }
  T <- 80
  x1 <- cumsum(rnorm(T)); x2 <- cumsum(rnorm(T))
  brk <- c(rep(0, 40), rep(2, 40))
  y <- 1 + brk + 0.8 * x1 + 0.4 * x2 + as.numeric(arima.sim(list(ar = 0.3), T))
  data.frame(y = y, x1 = x1, x2 = x2)
}

test_that("Kripfganz-Schneider bounds reproduce ardlbounds", {
  r <- fbardl:::.ks_bounds("F", 3, 2, 80, 4, value = 3.9)
  expect_equal(unname(r$cv["I0", ]), c(3.2043, 3.8812, 5.4116), tolerance = 1e-4)
  expect_equal(unname(r$cv["I1", ]), c(4.2082, 4.9866, 6.7118), tolerance = 1e-4)
  expect_equal(unname(r$pvalue), c(0.04903, 0.13032), tolerance = 1e-4)
  r <- fbardl:::.ks_bounds("t", 1, 1, NULL, value = -2.2)
  expect_equal(unname(r$pvalue), c(0.02666, 0.11369), tolerance = 1e-4)
})

test_that("fardl statistics reproduce Stata fbardl for PSS cases 2 to 5", {
  d <- .ts_data()
  ref <- list(`2` = c(5.7073, -4.6603, 8.5421), `3` = c(7.3772, -4.6603, 8.5421),
              `4` = c(7.1391, -5.2939, 11.8777), `5` = c(9.4648, -5.2939, 11.8777))
  for (cs in 2:5) {
    r <- suppressMessages(fbardl(y ~ x1 + x2, data = d, type = "fardl",
                                 maxlag = 3, maxk = 2, case = cs))
    expect_equal(c(r$F.overall, r$t.dependent, r$F.independent),
                 ref[[as.character(cs)]], tolerance = 1e-4)
  }
  r <- suppressMessages(fbardl(y ~ x1 + x2, data = d, type = "fardl",
                               maxlag = 3, maxk = 2, case = 3))
  expect_equal(unname(r$cointegration$Fov.pval), c(0.0013, 0.0059), tolerance = 0.05)
  expect_equal(r$model.fit$r.squared, 0.6224, tolerance = 1e-3)
  expect_equal(r$diagnostics$bg4.stat, 7.7943, tolerance = 1e-4)
  expect_equal(r$diagnostics$bp.stat, 0.0780, tolerance = 1e-3)
})

test_that("bootstrap recursion reproduces the data with the observed residuals", {
  d <- .ts_data()
  X <- cbind(x1 = d$x1, x2 = d$x2); T <- nrow(d)
  fs <- sin(2 * pi * (1:T) / T); fc <- cos(2 * pi * (1:T) / T)
  for (tp in c("fbardl_mcnown", "fbardl_bvz")) {
    dev <- fbardl:::.bootstrap_ardl_test(d$y, X, 1, c(3, 0), fs, fc, FALSE, 0L, 3,
                                         tp, 10, 1, 1, 1, c("L.y", "L.x1", "L.x2"),
                                         dgpcheck = TRUE)
    expect_true(all(dev < 1e-8))
  }
})
