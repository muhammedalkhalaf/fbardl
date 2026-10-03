## Monte Carlo size check (n = 100, y and x independent random walks,
## MC = 200, B = 199, 5% level). Slow: runs only when the environment
## variable FBARDL_SIZE_MC is "true" and never on CRAN. The values reported
## in NEWS.md were produced by this design (default call: Fov 0.075,
## t 0.135, Find 0.110, AND decision 0.040; fbardl 1.1.0, which held k* and
## the lags fixed in the bootstrap: Fov 0.375, t 0.420, Find 0.270).
test_that("bootstrap size with k* and lags selected again is far below the 1.1.0 size", {
  skip_on_cran()
  skip_if_not(identical(Sys.getenv("FBARDL_SIZE_MC"), "true"),
              "set FBARDL_SIZE_MC=true to run the size Monte Carlo")
  MC <- 200; B <- 199; n <- 100
  rej <- t(sapply(seq_len(MC), function(r) {
    set.seed(1000 + r)
    d <- data.frame(y = cumsum(rnorm(n)), x = cumsum(rnorm(n)))
    f <- suppressWarnings(suppressMessages(
      fbardl(y ~ x, d, maxlag = 1, reps = B, seed = r)))
    co <- f$cointegration
    c(co$Fov.pval < 0.05, co$t.pval < 0.05, co$Find.pval < 0.05,
      co$decision.code == "COINTEGRATION")
  }))
  size <- colMeans(rej)
  # Monte Carlo standard error about 0.02; t over-rejects (about 0.13)
  expect_lt(size[4], 0.075)
  expect_true(all(size[1:3] < 0.18))
})
