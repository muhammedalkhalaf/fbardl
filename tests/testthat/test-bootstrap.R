.boot_data <- function(seed = 11, n = 70) {
  set.seed(seed)
  x1 <- cumsum(rnorm(n)); x2 <- cumsum(rnorm(n))
  y <- 0.5 * x1 - 0.3 * x2 + c(rep(0, n %/% 2), rep(1, n - n %/% 2)) + rnorm(n)
  data.frame(y = y, x1 = x1, x2 = x2)
}

test_that("shared engine copy is version 1.1.0 and unchanged", {
  expect_identical(fbardl:::.abe_engine_version, "1.1.0")
  cand <- c(test_path("..", "..", "R", "ardl_boot_engine.R"),
            test_path("..", "..", "00_pkg_src", "fbardl", "R", "ardl_boot_engine.R"))
  cand <- cand[file.exists(cand)]
  skip_if(length(cand) == 0, "engine source file not reachable")
  expect_identical(unname(tools::md5sum(cand[1])), "bc08ba0a908a6cc92fb1ae491eced973")
})

test_that("bootstrap statistics equal the reported ones and the recursion reproduces the data", {
  d <- .boot_data()
  for (tp in c("fbardl_bvz", "fbardl_mcnown")) for (unc in c(FALSE, TRUE)) {
    r <- suppressMessages(fbardl(y ~ x1 + x2, data = d, type = tp, maxlag = 2,
                                 maxk = 2, reps = 9, seed = 1,
                                 unconditional = unc))
    co <- r$cointegration
    expect_true(all(co$dgpcheck < 1e-8))
    expect_equal(as.numeric(co$settings$j0), as.numeric(unc))
    expect_lt(co$design_check, 1e-8)
    expect_identical(co$engine.version, "1.1.0")
    expect_equal(sum(co$nvalid), 27)
    # p-values are shares of bootstrap statistics at least as extreme
    expect_equal(co$Fov.pval, mean(co$Fov.boot >= r$F.overall))
    expect_equal(co$t.pval, mean(co$t.boot <= r$t.dependent))
    expect_equal(co$Find.pval, mean(co$Find.boot >= r$F.independent))
    expect_identical(co$settings$nulls, if (tp == "fbardl_bvz") "separate" else "joint")
    expect_identical(co$settings$init, if (tp == "fbardl_bvz") "block" else "observed")
  }
  # the engine's canonical ECM gives the reported statistics (conditional model)
  r <- suppressMessages(fbardl(y ~ x1 + x2, data = d, maxlag = 2, maxk = 2,
                               reps = 5, seed = 2))
  X <- cbind(x1 = d$x1, x2 = d$x2)
  f <- fbardl:::.fbardl_fourier(nrow(d), r$best.kstar)
  D <- fbardl:::.abe_design(d$y, X, r$best.p, r$best.q, max(r$best.p, r$best.q) + 2,
                            3, cbind(f$sin, f$cos))
  expect_equal(unname(fbardl:::.abe_stats(D, 3)),
               c(r$F.overall, r$t.dependent, r$F.independent), tolerance = 1e-8)
})

test_that("what is selected again in each replication follows the arguments", {
  d <- .boot_data()
  sel <- function(...) suppressMessages(fbardl(y ~ x1 + x2, data = d, maxlag = 1,
                                               maxk = 2, reps = 3, seed = 1, ...))$cointegration$reselect
  expect_identical(sel(), c(kstar = TRUE, lags = TRUE))
  expect_identical(sel(kstar = 1), c(kstar = FALSE, lags = TRUE))
  expect_identical(sel(fourier = FALSE), c(kstar = FALSE, lags = TRUE))
  expect_identical(sel(kstar = 1, lags = list(p = 1, q = 1)), c(kstar = FALSE, lags = FALSE))
})

test_that("seed gives reproducible results and restores the session RNG", {
  d <- .boot_data()
  set.seed(99); u0 <- runif(1); set.seed(99)
  a <- suppressMessages(fbardl(y ~ x1 + x2, data = d, maxlag = 1, maxk = 2,
                               reps = 7, seed = 5))
  expect_identical(runif(1), u0)
  b <- suppressMessages(fbardl(y ~ x1 + x2, data = d, maxlag = 1, maxk = 2,
                               reps = 7, seed = 5))
  expect_identical(a$cointegration$Fov.boot, b$cointegration$Fov.boot)
})

test_that("decision uses the AND rule of the engine", {
  d <- .boot_data()
  r <- suppressMessages(fbardl(y ~ x1 + x2, data = d, maxlag = 1, maxk = 2,
                               reps = 19, seed = 3))
  co <- r$cointegration
  rej <- c(r$F.overall > co$F.cv["5%"], r$t.dependent < co$t.cv["5%"],
           r$F.independent > co$Find.cv["5%"])
  code <- if (!rej[1]) "NO_COINTEGRATION" else if (!rej[2]) "DEGENERATE_1" else
    if (!rej[3]) "DEGENERATE_2" else "COINTEGRATION"
  expect_identical(co$decision.code, code)
  out <- capture.output(print(r))
  expect_true(any(grepl("Selected again in each replication: k\\* and lags", out)))
})

test_that("re-selection on bootstrap samples equals an independent selection", {
  # independent design with data.frame and lm(); k* by minimum SSR with
  # every lag at maxlag, then p and q by AIC (AIC() differs from the
  # package's criterion by a constant, so the argmin is the same)
  hdesign <- function(y, X, p, q, k) {
    n <- length(y); dy <- c(NA, diff(y)); dX <- rbind(NA, diff(X))
    r <- (max(p, q) + 2):n
    df <- data.frame(dy = dy[r], Ly = y[r - 1])
    for (j in seq_len(ncol(X))) df[[paste0("Lx", j)]] <- X[r - 1, j]
    for (i in seq_len(p)) df[[paste0("dy", i)]] <- dy[r - i]
    for (j in seq_len(ncol(X))) for (l in 0:q[j]) df[[paste0("dx", j, "_", l)]] <- dX[r - l, j]
    if (k > 0) { df$s <- sin(2 * pi * k * r / n); df$c <- cos(2 * pi * k * r / n) }
    df
  }
  hselect <- function(y, X, maxlag, kv) {
    K <- ncol(X)
    ssr <- sapply(kv, function(k) sum(resid(lm(dy ~ ., hdesign(y, X, maxlag, rep(maxlag, K), k)))^2))
    k <- kv[which.min(ssr)]
    qg <- as.matrix(expand.grid(rep(list(0:maxlag), K)))
    best <- Inf
    for (p in seq_len(maxlag)) for (i in seq_len(nrow(qg))) {
      v <- AIC(lm(dy ~ ., hdesign(y, X, p, qg[i, ], k)))
      if (v < best - 1e-10) { best <- v; out <- list(p = p, q = unname(qg[i, ]), kstar = k) }
    }
    out
  }
  d <- .boot_data(seed = 5, n = 60)
  X <- cbind(x1 = d$x1, x2 = d$x2)
  setup <- list(maxlag = 2, hastrend = FALSE, j0 = 0L, case = 3, ic = "aic",
                fourier = TRUE, kstar = NULL, lags = NULL, kvalues = 1:3)
  sel <- fbardl:::.fbardl_select(d$y, X, setup)
  expect_equal(sel[c("p", "q", "kstar")], hselect(d$y, X, 2, 1:3))
  f <- fbardl:::.fbardl_fourier(nrow(d), sel$kstar)
  samples <- list()
  selfun <- function(yy, xx) {
    xx <- as.matrix(xx)
    colnames(xx) <- colnames(X)
    samples[[length(samples) + 1L]] <<- list(y = yy, x = xx)
    fbardl:::.fbardl_select(yy, xx, setup)[c("p", "q", "kstar")]
  }
  for (tp in c("separate", "joint")) {
    fbardl:::.ardl_boot_engine(d$y, X, p = sel$p, q = sel$q, case = 3,
                               det = cbind(f$sin, f$cos), nulls = tp,
                               xmodel = if (tp == "separate") "vecm" else "var",
                               B = 4, select = selfun, seed = 1)
  }
  expect_length(samples, 16)
  for (s in samples) {
    expect_false(isTRUE(all.equal(s$y, d$y)))
    expect_equal(fbardl:::.fbardl_select(s$y, s$x, setup)[c("p", "q", "kstar")],
                 hselect(s$y, s$x, 2, 1:3))
  }
})

test_that("KS bounds through fbardl with fourier = FALSE match Stata ardlbounds", {
  # n = 83 with lags p = 2, q = (0, 0) gives 80 observations, k = 2 and
  # sr = 4, the configuration of the Stata ardlbounds reference values
  d <- .boot_data(seed = 8, n = 83)
  r <- suppressMessages(fbardl(y ~ x1 + x2, data = d, type = "fardl", fourier = FALSE,
                               lags = list(p = 2, q = c(0, 0))))
  co <- r$cointegration
  expect_identical(r$nobs, 80L)
  expect_identical(co$sr, 4L)
  expect_equal(unname(co$F.cv["I0", ]), c(3.2043, 3.8812, 5.4116), tolerance = 1e-4)
  expect_equal(unname(co$F.cv["I1", ]), c(4.2082, 4.9866, 6.7118), tolerance = 1e-4)
  # p-values are consistent with the 5% bounds
  expect_identical(unname(co$Fov.pval[1] < 0.05), unname(r$F.overall > co$F.cv["I0", "5%"]))
  expect_identical(unname(co$Fov.pval[2] < 0.05), unname(r$F.overall > co$F.cv["I1", "5%"]))
  expect_identical(unname(co$t.pval[2] < 0.05), unname(r$t.dependent < co$t.cv["I1", "5%"]))
})

test_that("summary() stars the diagnostics when p is small", {
  d <- .boot_data()
  r <- suppressMessages(fbardl(y ~ x1 + x2, data = d, type = "fardl", fourier = FALSE,
                               maxlag = 1))
  out <- capture.output(summary(r))
  line <- grep("Breusch-Pagan", out, value = TRUE)
  expect_identical(grepl("\\*", line), r$diagnostics$bp.pval < 0.10)
  expect_identical(fbardl:::.stars(0.9177), "")
  expect_identical(fbardl:::.stars(0.03), " **")
})
