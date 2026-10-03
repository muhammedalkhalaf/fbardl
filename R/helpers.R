#' Build the ARDL equilibrium-correction design
#'
#' Rows are t = max(p, q) + 2, ..., n. Regressors: L.y, L.x, L1..Lp D.y,
#' D.x lags j0..q (per regressor), Fourier terms, trend (PSS cases 4 and 5),
#' constant.
#' @noRd
.build_ardl_data_flex <- function(y, X, p, q_vec, fourier_sin = NULL,
                                  fourier_cos = NULL, trend = FALSE, j0 = 0L) {
  n <- length(y)
  nindep <- ncol(X)
  if (length(q_vec) != nindep)
    stop("q_vec must have length equal to number of independent variables")
  start_idx <- max(p, max(q_vec)) + 2L
  if (start_idx > n) return(NULL)
  rows <- start_idx:n
  dy <- c(NA_real_, diff(y))
  dX <- rbind(NA_real_, diff(X))
  lagv <- function(v, j) v[rows - j]

  regressors <- list(y[rows - 1L])
  coef_names <- "L.y"
  for (j in seq_len(nindep)) {
    regressors[[length(regressors) + 1L]] <- X[rows - 1L, j]
    coef_names <- c(coef_names, paste0("L.", colnames(X)[j]))
  }
  for (lag in seq_len(p)) {
    regressors[[length(regressors) + 1L]] <- lagv(dy, lag)
    coef_names <- c(coef_names, paste0("L", lag, ".D.y"))
  }
  for (j in seq_len(nindep)) {
    if (q_vec[j] >= j0) for (lag in j0:q_vec[j]) {
      regressors[[length(regressors) + 1L]] <- lagv(dX[, j], lag)
      coef_names <- c(coef_names, if (lag == 0) paste0("D.", colnames(X)[j])
                      else paste0("L", lag, ".D.", colnames(X)[j]))
    }
  }
  if (!is.null(fourier_sin)) {
    regressors[[length(regressors) + 1L]] <- fourier_sin[rows]
    regressors[[length(regressors) + 1L]] <- fourier_cos[rows]
    coef_names <- c(coef_names, "sin", "cos")
  }
  if (trend) {
    regressors[[length(regressors) + 1L]] <- as.numeric(rows)
    coef_names <- c(coef_names, "trend")
  }
  regressors[[length(regressors) + 1L]] <- rep(1, length(rows))
  coef_names <- c(coef_names, "constant")

  Xmat <- do.call(cbind, regressors)
  colnames(Xmat) <- coef_names
  list(Y = dy[rows], Xmat = Xmat, coef_names = coef_names, rows = rows)
}


#' All candidate columns of .build_ardl_data_flex on every row t = 1..n
#' (NA where a lag is not available), for fast repeated selection
#' @noRd
.fbardl_allcols <- function(y, X, maxlag, j0 = 0L) {
  n <- length(y)
  K <- ncol(X)
  dy <- c(NA_real_, diff(y))
  dX <- rbind(NA_real_, diff(X))
  lagv <- function(v, j) c(rep(NA_real_, j), v[seq_len(n - j)])
  xn <- colnames(X)
  cols <- list(lagv(y, 1L))
  nm <- "L.y"
  for (j in seq_len(K)) {
    cols[[length(cols) + 1L]] <- lagv(X[, j], 1L)
    nm <- c(nm, paste0("L.", xn[j]))
  }
  for (lag in seq_len(maxlag)) {
    cols[[length(cols) + 1L]] <- lagv(dy, lag)
    nm <- c(nm, paste0("L", lag, ".D.y"))
  }
  # dX column of regressor j at lag l (l = j0..maxlag)
  idx_dx <- matrix(NA_integer_, K, maxlag + 1L)
  for (j in seq_len(K)) for (lag in j0:maxlag) {
    cols[[length(cols) + 1L]] <- lagv(dX[, j], lag)
    nm <- c(nm, if (lag == 0) paste0("D.", xn[j]) else paste0("L", lag, ".D.", xn[j]))
    idx_dx[j, lag + 1L] <- length(cols)
  }
  M <- do.call(cbind, cols)
  colnames(M) <- nm
  list(M = M, dy = dy, xn = xn, n = n, K = K, idx_dx = idx_dx)
}


#' Design of .build_ardl_data_flex taken from .fbardl_allcols (identical
#' values and column order). 'D' holds the deterministic columns that follow
#' the stochastic ones (Fourier sin and cos, trend, constant) on rows 1..n.
#' @noRd
.fbardl_subdesign <- function(A, p, q_vec, D, j0 = 0L, names = TRUE) {
  n <- A$n
  K <- A$K
  start_idx <- max(p, max(q_vec)) + 2L
  if (start_idx > n) return(NULL)
  rows <- start_idx:n
  idx <- c(seq_len(1L + K), 1L + K + seq_len(p))
  for (j in seq_len(K)) if (q_vec[j] >= j0) idx <- c(idx, A$idx_dx[j, (j0:q_vec[j]) + 1L])
  if (is.null(D)) {
    # A$M already holds the deterministic columns after the stochastic ones
    Xmat <- A$M[rows, c(idx, A$idx_det), drop = FALSE]
  } else {
    Xmat <- cbind(A$M[rows, idx, drop = FALSE], D[rows, , drop = FALSE])
  }
  nm <- if (names) colnames(Xmat) else NULL
  list(Y = A$dy[rows], Xmat = Xmat, coef_names = nm, rows = rows)
}


#' Deterministic columns of the design (sin, cos, trend, constant) on rows 1..n
#' @noRd
.fbardl_detcols <- function(n, fs, trend) {
  D <- cbind(sin = fs$sin, cos = fs$cos)
  if (trend) D <- cbind(D, trend = as.numeric(seq_len(n)))
  cbind(D, constant = rep(1, n))
}


#' Wald F-test
#' @noRd
.wald_f_test <- function(coefs, vcov_mat, R, r, df_resid) {
  # Compute (R*beta - r)'(R*V*R')^-1(R*beta - r) / q
  q <- nrow(R)
  Rb <- R %*% coefs - r
  RVR <- R %*% vcov_mat %*% t(R)
  F_stat <- as.numeric(t(Rb) %*% solve(RVR) %*% Rb / q)
  return(F_stat)
}


#' PSS bounds test with Kripfganz and Schneider (2020) critical values
#'
#' Finite-sample critical values and approximate p-values from the response
#' surface regressions of Kripfganz and Schneider (2020). The F_ind test has
#' no tabulated distribution; it is reported only with the bootstrap types.
#' With Fourier terms in the model (valid = FALSE) no valid bounds exist: the
#' bounds of the model without Fourier terms are returned for reference, with
#' no p-values and no decision.
#' @noRd
.pss_bounds_test <- function(Fov_stat, t_stat, k, case, nobs, sr, valid = TRUE) {
  Fb <- .ks_bounds("F", case, k, nobs, sr, value = if (valid) Fov_stat)
  tb <- .ks_bounds("t", case, k, nobs, sr, value = if (valid) t_stat)
  pF <- if (valid) Fb$pvalue else c(NA_real_, NA_real_)
  pt <- if (valid) tb$pvalue else c(NA_real_, NA_real_)
  dec <- function(stat, lo, hi, upper) {
    if (anyNA(c(lo, hi))) return(NA_character_)
    if (upper) {
      if (stat > hi) "reject" else if (stat < lo) "do not reject" else "inconclusive"
    } else {
      if (stat < hi) "reject" else if (stat > lo) "do not reject" else "inconclusive"
    }
  }
  if (!valid) {
    code <- "NOT_AVAILABLE"
    decision <- paste("NOT AVAILABLE: no valid bounds exist with Fourier terms;",
                      "use type = \"fbardl_bvz\" or \"fbardl_mcnown\"")
  } else {
    F_dec <- dec(Fov_stat, Fb$cv["I0", "5%"], Fb$cv["I1", "5%"], TRUE)
    t_dec <- dec(t_stat, tb$cv["I0", "5%"], tb$cv["I1", "5%"], FALSE)
    if (anyNA(c(F_dec, t_dec))) {
      code <- "UNAVAILABLE"
      decision <- "Critical values unavailable (fewer than twice as many observations as coefficients)"
    } else if (F_dec == "reject" && t_dec == "reject") {
      code <- "COINTEGRATION"
      decision <- "COINTEGRATION: F and t beyond the I(1) bounds at 5%"
    } else if (F_dec == "do not reject" || t_dec == "do not reject") {
      code <- "NO_COINTEGRATION"
      decision <- "NO COINTEGRATION: F or t within the I(0) bound at 5%"
    } else {
      code <- "INCONCLUSIVE"
      decision <- "INCONCLUSIVE at 5% (statistic between the bounds)"
    }
  }
  list(
    source = if (valid) "Kripfganz and Schneider (2020)" else
      paste("Kripfganz and Schneider (2020): bounds for the model without",
            "Fourier terms; not valid with Fourier terms"),
    valid = valid,
    sr = sr,
    F.cv = Fb$cv, t.cv = tb$cv,
    Fov.pval = pF, t.pval = pt, Find.pval = NA_real_,
    F.cv05.I0 = unname(Fb$cv["I0", "5%"]), F.cv05.I1 = unname(Fb$cv["I1", "5%"]),
    t.cv05.I0 = unname(tb$cv["I0", "5%"]), t.cv05.I1 = unname(tb$cv["I1", "5%"]),
    Find.cv05 = NA_real_,
    decision = decision, decision.code = code)
}


#' Fourier terms sin(2 pi k t / n), cos(2 pi k t / n), t = 1..n
#' @noRd
.fbardl_fourier <- function(n, k) {
  tt <- 1:n
  if (k > 0) list(sin = sin(2 * pi * k * tt / n), cos = cos(2 * pi * k * tt / n))
  else list(sin = NULL, cos = NULL)
}


#' Selection of k* (minimum SSR, every lag at maxlag) and of (p, q) (AIC or
#' BIC, k* fixed), as in fbardl 1.1.0. A k* or lags given by the user are
#' used as they are. The same function is used on the data and on every
#' bootstrap sample.
#' @noRd
.fbardl_select <- function(y, X, setup) {
  n <- length(y)
  nindep <- ncol(X)
  maxlag <- setup$maxlag
  hastrend <- setup$hastrend
  j0 <- setup$j0
  A <- .fbardl_allcols(y, X, max(maxlag, setup$lags$p, setup$lags$q), j0)
  ssr_of <- function(d) {
    if (is.null(d) || length(d$Y) < 10) return(NULL)
    if (!all(is.finite(d$Xmat)) || !all(is.finite(d$Y))) return(NULL)
    f <- stats::.lm.fit(d$Xmat, d$Y)
    list(rss = sum(f$residuals^2), nobs = length(f$residuals), k = ncol(d$Xmat))
  }
  ## Step 1: k*
  ssr_by_k <- NULL
  best_ssr_k <- NA_real_
  if (setup$fourier && !is.null(setup$kstar)) {
    best_kstar <- setup$kstar
  } else {
    kv <- if (setup$fourier) setup$kvalues else 0
    best_kstar <- 0
    best_ssr_k <- Inf
    ssr <- numeric(length(kv))
    for (i in seq_along(kv)) {
      fs <- .fbardl_fourier(n, kv[i])
      s <- ssr_of(.fbardl_subdesign(A, maxlag, rep(maxlag, nindep),
                                    .fbardl_detcols(n, fs, hastrend), j0, FALSE))
      ssr[i] <- if (is.null(s)) Inf else s$rss
      if (ssr[i] < best_ssr_k) {
        best_ssr_k <- ssr[i]
        best_kstar <- kv[i]
      }
    }
    if (!setup$fourier) best_kstar <- 0
    ssr_by_k <- data.frame(k = kv, ssr = ssr)
  }
  ## Step 2: lags
  fs <- .fbardl_fourier(n, best_kstar)
  AD <- A
  AD$M <- cbind(A$M, .fbardl_detcols(n, fs, hastrend))
  AD$idx_det <- (ncol(A$M) + 1L):ncol(AD$M)
  icval <- function(s) {
    ll <- -s$nobs/2 * (log(2 * pi) + log(s$rss/s$nobs) + 1)
    if (setup$ic == "aic") -2 * ll + 2 * s$k else -2 * ll + s$k * log(s$nobs)
  }
  if (!is.null(setup$lags)) {
    s <- ssr_of(.fbardl_subdesign(AD, setup$lags$p, setup$lags$q, NULL, j0, FALSE))
    if (is.null(s)) stop("the model with the lags given in 'lags' cannot be estimated")
    return(list(p = setup$lags$p, q = setup$lags$q, kstar = best_kstar,
                ic = icval(s), ssr_by_k = ssr_by_k, best_ssr_k = best_ssr_k,
                total_specs = 1L))
  }
  best_ic_val <- Inf
  best_p <- 1
  best_q <- rep(0, nindep)
  total_specs <- 0
  q_combinations <- as.matrix(expand.grid(replicate(nindep, 0:maxlag, simplify = FALSE)))
  for (p in 1:maxlag) {
    for (qidx in seq_len(nrow(q_combinations))) {
      total_specs <- total_specs + 1
      q_vec <- as.integer(q_combinations[qidx, ])
      s <- ssr_of(.fbardl_subdesign(AD, p, q_vec, NULL, j0, FALSE))
      if (is.null(s)) next
      ic_tmp <- icval(s)
      if (ic_tmp < best_ic_val) {
        best_ic_val <- ic_tmp
        best_p <- p
        best_q <- q_vec
      }
    }
  }
  list(p = best_p, q = best_q, kstar = best_kstar, ic = best_ic_val,
       ssr_by_k = ssr_by_k, best_ssr_k = best_ssr_k, total_specs = total_specs)
}


#' Fov, t and Find on a design built by .build_ardl_data_flex (Wald F with
#' the OLS covariance, as for the reported statistics)
#' @noRd
.fbardl_stats <- function(d, fov_names, ind_names) {
  if (!all(is.finite(d$Xmat)) || !all(is.finite(d$Y))) stop("non-finite design")
  qx <- qr(d$Xmat)
  if (qx$rank < ncol(d$Xmat)) stop("rank-deficient design")
  b <- qr.coef(qx, d$Y)
  e <- qr.resid(qx, d$Y)
  piv <- order(qx$pivot)
  V <- sum(e^2) / (length(e) - ncol(d$Xmat)) *
    chol2inv(qr.R(qx))[piv, piv, drop = FALSE]
  Ftest <- function(nms) {
    i <- match(nms, d$coef_names)
    as.numeric(t(b[i]) %*% solve(V[i, i, drop = FALSE]) %*% b[i]) / length(i)
  }
  iy <- match("L.y", d$coef_names)
  c(Fov = Ftest(fov_names), t = unname(b[iy] / sqrt(V[iy, iy])),
    Find = Ftest(ind_names))
}


#' Bootstrap ARDL test through the shared engine (R/ardl_boot_engine.R)
#'
#' "fbardl_bvz": separate nulls for Fov, t and Find (Bertelli, Vacca and
#' Zoia, 2022), marginal VECM for Delta x, residuals recentred after each
#' draw, random block of initial values. "fbardl_mcnown": the Fov null for
#' all statistics (McNown, Sam and Goh, 2018, Steps 1-8), applied to the
#' conditional ECM unless unconditional = TRUE (j0 = 1, MSG eq. 12);
#' unrestricted Delta x equation with y_{t-1}, residuals recentred once by
#' their mean, observed initial values.
#' The restricted equations are estimated once with the data-selected k*, p
#' and q; the Fourier terms of the data-generating process are held fixed.
#' On every bootstrap sample k* and the lags are selected again unless the
#' user fixed them, and the statistics are computed exactly as on the data.
#' @noRd
.fbardl_bootstrap <- function(y, X, sel, setup, type, reps, fov_names,
                              ind_names, seed = NULL) {
  n <- length(y)
  xn <- colnames(X)
  fs <- .fbardl_fourier(n, sel$kstar)
  det <- if (sel$kstar > 0) cbind(sin = fs$sin, cos = fs$cos) else NULL
  stat_fun <- function(yy, xx, sp) {
    xx <- as.matrix(xx)
    colnames(xx) <- xn
    f <- .fbardl_fourier(length(yy), sp$kstar)
    d <- .build_ardl_data_flex(yy, xx, sp$p, sp$q, f$sin, f$cos,
                               setup$hastrend, setup$j0)
    .fbardl_stats(d, fov_names, ind_names)
  }
  need_k <- setup$fourier && is.null(setup$kstar)
  need_l <- is.null(setup$lags)
  selfun <- NULL
  if (need_k || need_l) {
    fix <- setup
    selfun <- function(yy, xx) {
      xx <- as.matrix(xx)
      colnames(xx) <- xn
      .fbardl_select(yy, xx, fix)[c("p", "q", "kstar")]
    }
  }
  bvz <- type == "fbardl_bvz"
  run <- function() .ardl_boot_engine(
    y, X, p = sel$p, q = sel$q, case = setup$case, det = det, j0 = setup$j0,
    nulls = if (bvz) "separate" else "joint",
    xmodel = if (bvz) "vecm" else "var",
    B = reps,
    init = if (bvz) "block" else "observed",
    recentre = if (bvz) "draw" else "once",
    stat_fun = stat_fun, select = selfun, spec = sel[c("p", "q", "kstar")],
    use_find = TRUE, level = 0.05, seed = seed)
  eng <- run()

  nm <- c("10%", "5%", "2.5%", "1%")
  F_cv <- stats::setNames(eng$cv[, "Fov"], nm)
  t_cv <- stats::setNames(eng$cv[, "t"], nm)
  I_cv <- stats::setNames(eng$cv[, "Find"], nm)
  list(
    source = "bootstrap",
    scheme = if (bvz) "Bertelli, Vacca and Zoia (2022)" else "McNown, Sam and Goh (2018)",
    Fov.pval = unname(eng$p_value["Fov"]), t.pval = unname(eng$p_value["t"]),
    Find.pval = unname(eng$p_value["Find"]),
    F.cv = F_cv, t.cv = t_cv, Find.cv = I_cv,
    F.cv05.I0 = NA_real_, F.cv05.I1 = unname(F_cv["5%"]),
    t.cv05.I0 = NA_real_, t.cv05.I1 = unname(t_cv["5%"]),
    Find.cv05 = unname(I_cv["5%"]),
    nvalid = eng$n_valid, nfail = eng$n_fail,
    decision = eng$label, decision.code = eng$decision,
    reselect = c(kstar = need_k, lags = need_l),
    dgpcheck = eng$dgpcheck, design_check = eng$design_check,
    engine.version = eng$engine_version,
    settings = eng$settings,
    Fov.boot = unname(eng$boot[, "Fov"]), t.boot = unname(eng$boot[, "t"]),
    Find.boot = unname(eng$boot[, "Find"]))
}


#' Compute long-run coefficients via delta method
#' @noRd
.compute_long_run <- function(coefs, vcov_mat, coef_names, indepvars,
                               ecm_idx, indep_level_idx, df_resid, level) {
  nindep <- length(indepvars)

  lr_estimates <- numeric(nindep)
  lr_se <- numeric(nindep)
  lr_z <- numeric(nindep)
  lr_pval <- numeric(nindep)
  lr_lo <- numeric(nindep)
  lr_hi <- numeric(nindep)

  alpha <- coefs[ecm_idx]

  for (j in 1:nindep) {
    idx <- indep_level_idx[j]
    beta_x <- coefs[idx]

    # Long-run coefficient: -beta_x / alpha
    lr_estimates[j] <- -beta_x / alpha

    # Delta method for standard error
    # Gradient: d/d(alpha) = beta_x / alpha^2
    #           d/d(beta_x) = -1 / alpha
    grad <- numeric(length(coefs))
    grad[ecm_idx] <- beta_x / alpha^2
    grad[idx] <- -1 / alpha

    lr_var <- t(grad) %*% vcov_mat %*% grad
    lr_se[j] <- sqrt(as.numeric(lr_var))

    lr_z[j] <- lr_estimates[j] / lr_se[j]
    lr_pval[j] <- 2 * pnorm(-abs(lr_z[j]))

    z_crit <- qnorm(1 - (1 - level) / 2)
    lr_lo[j] <- lr_estimates[j] - z_crit * lr_se[j]
    lr_hi[j] <- lr_estimates[j] + z_crit * lr_se[j]
  }

  result <- data.frame(
    estimate = lr_estimates,
    std.error = lr_se,
    z.value = lr_z,
    p.value = lr_pval,
    ci.lower = lr_lo,
    ci.upper = lr_hi
  )
  rownames(result) <- indepvars

  return(result)
}


#' Extract short-run coefficients
#' @noRd
.extract_short_run <- function(coefs, se, tvals, pvals, coef_names,
                                depvar, indepvars, p, q_vec) {

  # Find short-run coefficient indices (lagged differences)
  sr_pattern <- paste0("^L[0-9]+\\.D\\.|^D\\.")
  sr_idx <- grep(sr_pattern, coef_names)

  # Also include Fourier terms and constant
  other_idx <- which(coef_names %in% c("sin", "cos", "constant"))
  all_sr_idx <- c(sr_idx, other_idx)

  result <- data.frame(
    estimate = coefs[all_sr_idx],
    std.error = se[all_sr_idx],
    t.value = tvals[all_sr_idx],
    p.value = pvals[all_sr_idx]
  )
  rownames(result) <- coef_names[all_sr_idx]

  return(result)
}


#' Diagnostic tests
#' @noRd
.run_diagnostics <- function(residuals, nobs, nparams, Xmat, fitted) {
  n <- length(residuals)

  # Jarque-Bera normality test
  m3 <- mean((residuals - mean(residuals))^3)
  m2 <- mean((residuals - mean(residuals))^2)
  skew <- m3 / m2^1.5
  m4 <- mean((residuals - mean(residuals))^4)
  kurt <- m4 / m2^2 - 3
  jb_stat <- n * (skew^2 / 6 + kurt^2 / 24)
  jb_pval <- pchisq(jb_stat, 2, lower.tail = FALSE)

  # Breusch-Godfrey LM test (as Stata's estat bgodfrey): regression of the
  # residuals on the original regressors and l lagged residuals (pre-sample
  # values set to zero); statistic n R^2 ~ chi2(l)
  bgtest <- function(l) tryCatch({
    elag <- sapply(seq_len(l), function(j) c(rep(0, j), residuals[seq_len(n - j)]))
    Z <- cbind(Xmat, elag)
    u <- stats::lm.fit(Z, residuals)$residuals
    r2 <- 1 - sum(u^2) / sum((residuals - mean(residuals))^2)
    stat <- n * r2
    list(stat = stat, pval = pchisq(stat, l, lower.tail = FALSE))
  }, error = function(e) list(stat = NA, pval = NA))
  bg1 <- bgtest(1)
  bg4 <- if (n > 8) bgtest(4) else list(stat = NA, pval = NA)

  # Breusch-Pagan / Cook-Weisberg test on the fitted values (as Stata's
  # estat hettest): ESS / 2 from the regression of e^2 / (RSS / n) on yhat
  bp <- tryCatch({
    g <- residuals^2 / (sum(residuals^2) / n)
    fb <- stats::lm.fit(cbind(1, fitted), g)
    ess <- sum((g - mean(g))^2) - sum(fb$residuals^2)
    stat <- ess / 2
    list(stat = stat, pval = pchisq(stat, 1, lower.tail = FALSE))
  }, error = function(e) list(stat = NA, pval = NA))

  # ARCH(1) LM test (as Stata's estat archlm): (n - 1) R^2 from the
  # regression of e_t^2 on e_{t-1}^2
  arch1 <- tryCatch({
    e2 <- residuals^2
    fa <- stats::lm.fit(cbind(1, e2[-n]), e2[-1])
    r2 <- 1 - sum(fa$residuals^2) / sum((e2[-1] - mean(e2[-1]))^2)
    stat <- (n - 1) * r2
    list(stat = stat, pval = pchisq(stat, 1, lower.tail = FALSE))
  }, error = function(e) list(stat = NA, pval = NA))

  return(list(
    jb.stat = jb_stat,
    jb.pval = jb_pval,
    bg1.stat = bg1$stat,
    bg1.pval = bg1$pval,
    bg4.stat = bg4$stat,
    bg4.pval = bg4$pval,
    bp.stat = bp$stat,
    bp.pval = bp$pval,
    arch1.stat = arch1$stat,
    arch1.pval = arch1$pval
  ))
}
