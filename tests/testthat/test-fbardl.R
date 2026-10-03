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
                                 maxlag = 3, maxk = 2, case = cs,
                                 kgrid = "fractional"))
    expect_equal(c(r$F.overall, r$t.dependent, r$F.independent),
                 ref[[as.character(cs)]], tolerance = 1e-4)
  }
  r <- suppressMessages(fbardl(y ~ x1 + x2, data = d, type = "fardl",
                               maxlag = 3, maxk = 2, case = 3,
                               kgrid = "fractional"))
  expect_equal(r$model.fit$r.squared, 0.6224, tolerance = 1e-3)
  expect_equal(r$diagnostics$bg4.stat, 7.7943, tolerance = 1e-4)
  expect_equal(r$diagnostics$bp.stat, 0.0780, tolerance = 1e-3)
})

test_that("fardl gives no verdict with Fourier terms and the KS decision without", {
  d <- .ts_data()
  r <- suppressMessages(fbardl(y ~ x1 + x2, data = d, type = "fardl",
                               maxlag = 3, maxk = 2, kgrid = "fractional"))
  expect_true(r$best.kstar > 0)
  co <- r$cointegration
  expect_false(co$valid)
  expect_identical(co$decision.code, "NOT_AVAILABLE")
  expect_true(all(is.na(c(co$Fov.pval, co$t.pval))))
  # bounds of the model without Fourier terms: sr excludes sin and cos
  expect_identical(co$sr, r$nparams - 1L - 3L - 2L)
  ks <- fbardl:::.ks_bounds("F", 3, 2, r$nobs, co$sr)
  expect_equal(co$F.cv, ks$cv)
  out <- capture.output(print(r))
  expect_true(any(grepl("without Fourier terms; not valid with Fourier terms", out)))
  expect_true(any(grepl("NOT AVAILABLE", out)))

  r0 <- suppressMessages(fbardl(y ~ x1 + x2, data = d, type = "fardl",
                                maxlag = 3, fourier = FALSE))
  co0 <- r0$cointegration
  expect_true(co0$valid)
  expect_identical(co0$sr, r0$nparams - 1L - 3L)
  ks <- fbardl:::.ks_bounds("F", 3, 2, r0$nobs, co0$sr, value = r0$F.overall)
  expect_equal(co0$Fov.pval, ks$pvalue)
  expect_true(co0$decision.code %in% c("COINTEGRATION", "NO_COINTEGRATION", "INCONCLUSIVE"))
})

test_that("fast selection designs equal .build_ardl_data_flex", {
  set.seed(3)
  n <- 40
  X <- cbind(a = cumsum(rnorm(n)), b = cumsum(rnorm(n)))
  y <- cumsum(rnorm(n))
  for (j0 in 0:1) for (p in 0:3) for (q1 in 0:3) for (tr in c(FALSE, TRUE)) {
    f <- fbardl:::.fbardl_fourier(n, 1.3)
    A <- fbardl:::.fbardl_allcols(y, X, 3, j0)
    d1 <- fbardl:::.build_ardl_data_flex(y, X, p, c(q1, 3 - q1), f$sin, f$cos, tr, j0)
    d2 <- fbardl:::.fbardl_subdesign(A, p, c(q1, 3 - q1),
                                     fbardl:::.fbardl_detcols(n, f, tr), j0)
    expect_identical(d1$Xmat, d2$Xmat)
    expect_identical(d1$Y, d2$Y)
  }
})

test_that("k* grids: integer by default, fractional on request", {
  d <- .ts_data()
  r <- suppressMessages(fbardl(y ~ x1 + x2, data = d, type = "fardl",
                               maxlag = 2, maxk = 3))
  expect_identical(r$ssr.by.k$k, 1:3)
  expect_identical(r$kgrid, "integer")
  # hand check: k* minimises the SSR of the model with every lag at maxlag
  ssr <- sapply(1:3, function(k) {
    f <- fbardl:::.fbardl_fourier(nrow(d), k)
    dd <- fbardl:::.build_ardl_data_flex(d$y, cbind(x1 = d$x1, x2 = d$x2), 2,
                                         c(2, 2), f$sin, f$cos)
    sum(lm.fit(dd$Xmat, dd$Y)$residuals^2)
  })
  expect_equal(r$ssr.by.k$ssr, ssr, tolerance = 1e-10)
  expect_equal(r$best.kstar, which.min(ssr))
  rf <- suppressMessages(fbardl(y ~ x1 + x2, data = d, type = "fardl",
                                maxlag = 2, maxk = 3, kgrid = "fractional"))
  expect_equal(rf$ssr.by.k$k, seq(0.1, 3, by = 0.1))
  rk <- suppressMessages(fbardl(y ~ x1 + x2, data = d, type = "fardl",
                                maxlag = 2, kstar = 1.7, lags = list(p = 1, q = 0)))
  expect_identical(rk$best.kstar, 1.7)
  expect_identical(unname(rk$best.q), c(0L, 0L))
  expect_error(fbardl(y ~ x1 + x2, data = d, kstar = 1, fourier = FALSE), "kstar")
  expect_error(fbardl(y ~ x1 + x2, data = d, maxk = 2.5), "maxk")
})
