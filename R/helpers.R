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
#' @noRd
.pss_bounds_test <- function(Fov_stat, t_stat, k, case, nobs, sr) {
  Fb <- .ks_bounds("F", case, k, nobs, sr, value = Fov_stat)
  tb <- .ks_bounds("t", case, k, nobs, sr, value = t_stat)
  dec <- function(stat, lo, hi, upper) {
    if (anyNA(c(lo, hi))) return(NA_character_)
    if (upper) {
      if (stat > hi) "reject" else if (stat < lo) "do not reject" else "inconclusive"
    } else {
      if (stat < hi) "reject" else if (stat > lo) "do not reject" else "inconclusive"
    }
  }
  F_dec <- dec(Fov_stat, Fb$cv["I0", "5%"], Fb$cv["I1", "5%"], TRUE)
  t_dec <- dec(t_stat, tb$cv["I0", "5%"], tb$cv["I1", "5%"], FALSE)
  decision <- if (anyNA(c(F_dec, t_dec))) {
    "Critical values unavailable (fewer than twice as many observations as coefficients)"
  } else if (F_dec == "reject" && t_dec == "reject") {
    "COINTEGRATION: F and t beyond the I(1) bounds at 5%"
  } else if (F_dec == "do not reject" || t_dec == "do not reject") {
    "NO COINTEGRATION: F or t within the I(0) bound at 5%"
  } else {
    "INCONCLUSIVE at 5% (statistic between the bounds)"
  }
  list(
    source = "Kripfganz and Schneider (2020)",
    sr = sr,
    F.cv = Fb$cv, t.cv = tb$cv,
    Fov.pval = Fb$pvalue, t.pval = tb$pvalue, Find.pval = NA_real_,
    F.cv05.I0 = unname(Fb$cv["I0", "5%"]), F.cv05.I1 = unname(Fb$cv["I1", "5%"]),
    t.cv05.I0 = unname(tb$cv["I0", "5%"]), t.cv05.I1 = unname(tb$cv["I1", "5%"]),
    Find.cv05 = NA_real_,
    decision = decision)
}


#' Bootstrap ARDL test (McNown, Sam and Goh 2018; Bertelli, Vacca and Zoia 2022)
#'
#' Port of _fbardl_bootstrap.ado (Stata fbardl 1.3.0). Bootstrap data are
#' generated recursively from the restricted y equation (one null for the
#' McNown et al. version, one per statistic for Bertelli et al.) and the
#' equations for Delta x, resampling the residual pairs; the full model is
#' re-estimated on each bootstrap sample.
#' @noRd
.bootstrap_ardl_test <- function(y, X, best_p, best_q, fsin, fcos, trend, j0,
                                 case, type, reps, Fov_stat, t_stat, Find_stat,
                                 fov_names, dgpcheck = FALSE) {
  T <- length(y); K <- ncol(X); p <- best_p
  xn <- colnames(X)
  des <- function(yy, XX) .build_ardl_data_flex(yy, XX, p, best_q, fsin, fcos, trend, j0)
  full <- des(y, X)
  cn <- full$coef_names
  ind_names <- paste0("L.", xn)

  nnull <- if (type == "fbardl_mcnown") 1L else 3L
  nulluse <- if (type == "fbardl_mcnown") c(1L, 1L, 1L) else 1:3
  drops <- list(fov_names, "L.y", ind_names)

  ## restricted y equations
  eqY <- lapply(seq_len(nnull), function(h) {
    keep <- setdiff(cn, drops[[h]])
    Z <- full$Xmat[, keep, drop = FALSE]
    b <- stats::lm.fit(Z, full$Y)$coefficients
    b[is.na(b)] <- 0
    coef <- stats::setNames(rep(0, length(cn)), cn)
    coef[keep] <- b
    xb <- rep(NA_real_, T); r <- rep(NA_real_, T)
    xb[full$rows] <- as.numeric(Z %*% b)
    r[full$rows] <- full$Y - xb[full$rows]
    list(coef = coef, xb = xb, r = r)
  })

  ## equations for Delta x
  dy0 <- c(0, diff(y)); dX0 <- rbind(0, diff(X))
  rowsx <- (p + 2L):T
  Zx <- cbind(if (type == "fbardl_mcnown") y[rowsx - 1L],
              X[rowsx - 1L, , drop = FALSE],
              do.call(cbind, lapply(seq_len(p), function(j)
                cbind(dy0[rowsx - j], dX0[rowsx - j, , drop = FALSE]))),
              if (!is.null(fsin)) cbind(fsin[rowsx], fcos[rowsx]),
              if (trend) as.numeric(rowsx), 1)
  eqX <- lapply(seq_len(K), function(m) {
    b <- stats::lm.fit(Zx, dX0[rowsx, m])$coefficients
    b[is.na(b)] <- 0
    xb <- rep(NA_real_, T); r <- rep(NA_real_, T)
    xb[rowsx] <- as.numeric(Zx %*% b)
    r[rowsx] <- dX0[rowsx, m] - xb[rowsx]
    off <- if (type == "fbardl_mcnown") 1L else 0L
    list(by = if (off) b[1] else 0, bx = b[off + seq_len(K)],
         phi = b[off + K + (seq_len(p) - 1L) * (K + 1L) + 1L],
         th = matrix(b[off + K + outer(seq_len(K) + 1L, (seq_len(p) - 1L) * (K + 1L), `+`)],
                     K, p),
         xb = xb, r = r)
  })

  ## history parts evaluated on any series
  fity <- function(Ys, Xs, dYs, dXs, t, co) {
    v <- co[["L.y"]] * Ys[t - 1L] + sum(co[paste0("L.", xn)] * Xs[t - 1L, ])
    for (j in seq_len(p)) v <- v + co[[paste0("L", j, ".D.y")]] * dYs[t - j]
    for (m in seq_len(K)) if (best_q[m] >= j0) for (j in j0:best_q[m]) {
      nm <- if (j == 0) paste0("D.", xn[m]) else paste0("L", j, ".D.", xn[m])
      v <- v + co[[nm]] * dXs[t - j, m]
    }
    v
  }
  fitx <- function(Ys, Xs, dYs, dXs, t, e) {
    v <- e$by * Ys[t - 1L] + sum(e$bx * Xs[t - 1L, ])
    for (j in seq_len(p)) v <- v + e$phi[j] * dYs[t - j] + sum(e$th[, j] * dXs[t - j, ])
    v
  }

  ok <- Reduce(`&`, c(lapply(eqY, function(e) !is.na(e$r)),
                      lapply(eqX, function(e) !is.na(e$r))))
  pool <- which(ok)
  t0 <- min(pool); tN <- max(pool)
  detY <- sapply(eqY, function(e) {
    d <- e$xb
    for (t in pool) d[t] <- d[t] - fity(y, X, dy0, dX0, t, e$coef)
    d
  })
  detY <- matrix(detY, T)
  detX <- sapply(eqX, function(e) {
    d <- e$xb
    for (t in pool) d[t] <- d[t] - fitx(y, X, dy0, dX0, t, e)
    d
  })
  detX <- matrix(detX, T)
  poolY <- sapply(eqY, `[[`, "r"); poolY <- matrix(poolY, T)
  poolX <- sapply(eqX, `[[`, "r"); poolX <- matrix(poolX, T)
  npool <- length(pool)

  recurse <- function(eY, eX, h) {
    Ys <- y; Xs <- X; dYs <- dy0; dXs <- dX0
    for (t in t0:tN) {
      for (i in seq_len(K)) {
        dxi <- detX[t, i] + eX[t, i] + fitx(Ys, Xs, dYs, dXs, t, eqX[[i]])
        dXs[t, i] <- dxi
        Xs[t, i] <- Xs[t - 1L, i] + dxi
      }
      dyt <- detY[t, h] + eY[t] + fity(Ys, Xs, dYs, dXs, t, eqY[[h]]$coef)
      dYs[t] <- dyt
      Ys[t] <- Ys[t - 1L] + dyt
    }
    list(y = Ys, X = Xs)
  }

  if (dgpcheck) {
    dev <- vapply(seq_len(nnull), function(h) {
      eY <- rep(0, T); eX <- matrix(0, T, K)
      eY[t0:tN] <- poolY[t0:tN, h]; eX[t0:tN, ] <- poolX[t0:tN, ]
      r <- recurse(eY, eX, h)
      max(abs(r$y - y), abs(r$X - X))
    }, numeric(1))
    return(dev)
  }

  # McNown et al.: residuals recentred once (Stata fbardl recmode 1)
  if (type == "fbardl_mcnown") {
    div <- npool - p - 1
    if (div < 1) div <- npool
    for (h in seq_len(nnull)) poolY[pool, h] <- poolY[pool, h] - sum(poolY[pool, h]) / div
    for (m in seq_len(K)) poolX[pool, m] <- poolX[pool, m] - sum(poolX[pool, m]) / div
  }

  stat_fun <- function(yy, XX) {
    d <- des(yy, XX)
    fit <- stats::lm.fit(d$Xmat, d$Y)
    if (fit$rank < ncol(d$Xmat)) return(c(NA, NA, NA))
    e <- fit$residuals
    s2 <- sum(e^2) / (length(e) - ncol(d$Xmat))
    V <- s2 * chol2inv(qr.R(fit$qr))
    b <- fit$coefficients
    Ftest <- function(nms) {
      i <- match(nms, d$coef_names)
      as.numeric(t(b[i]) %*% solve(V[i, i, drop = FALSE]) %*% b[i]) / length(i)
    }
    iy <- match("L.y", d$coef_names)
    c(Ftest(fov_names), b[iy] / sqrt(V[iy, iy]), Ftest(ind_names))
  }

  stats_mat <- matrix(NA_real_, reps, 3)
  for (bb in seq_len(reps)) {
    last <- 0L; cur <- NULL
    for (s in 1:3) {
      h <- nulluse[s]
      if (h != last) {
        idx <- pool[ceiling(stats::runif(npool) * npool)]
        eY <- rep(0, T); eX <- matrix(0, T, K)
        eY[t0:tN] <- poolY[idx, h]; eX[t0:tN, ] <- poolX[idx, , drop = FALSE]
        if (type != "fbardl_mcnown") {
          eY[t0:tN] <- eY[t0:tN] - mean(eY[t0:tN])
          eX[t0:tN, ] <- sweep(eX[t0:tN, , drop = FALSE], 2,
                               colMeans(eX[t0:tN, , drop = FALSE]))
        }
        r <- recurse(eY, eX, h)
        cur <- tryCatch(stat_fun(r$y, r$X), error = function(e) c(NA, NA, NA))
        last <- h
      }
      stats_mat[bb, s] <- cur[s]
    }
  }

  qhi <- function(v, pr) { v <- sort(v[!is.na(v)]); if (length(v) < 3) NA else v[min(ceiling(pr * length(v)), length(v))] }
  qlo <- function(v, pr) { v <- sort(v[!is.na(v)]); if (length(v) < 3) NA else v[max(floor(pr * length(v)), 1)] }
  Fb <- stats_mat[, 1]; tb <- stats_mat[, 2]; Ib <- stats_mat[, 3]
  lv <- c(0.10, 0.05, 0.025, 0.01)
  F_cv <- sapply(1 - lv, qhi, v = Fb); t_cv <- sapply(lv, qlo, v = tb)
  I_cv <- sapply(1 - lv, qhi, v = Ib)
  names(F_cv) <- names(t_cv) <- names(I_cv) <- c("10%", "5%", "2.5%", "1%")
  Fov_pval <- mean(Fb[!is.na(Fb)] >= Fov_stat)
  t_pval <- mean(tb[!is.na(tb)] <= t_stat)
  Find_pval <- mean(Ib[!is.na(Ib)] >= Find_stat)

  if (Fov_pval < 0.05 && t_pval < 0.05 && Find_pval < 0.05) {
    decision <- "COINTEGRATION detected (all tests significant at 5%)"
  } else if (Fov_pval >= 0.05 && t_pval >= 0.05 && Find_pval >= 0.05) {
    decision <- "NO COINTEGRATION detected at 5% level"
  } else if (Fov_pval < 0.05 && Find_pval < 0.05 && t_pval >= 0.05) {
    decision <- "DEGENERATE CASE #1: Fov & Find significant but t not (y may be I(0))"
  } else if (Fov_pval < 0.05 && t_pval < 0.05 && Find_pval >= 0.05) {
    decision <- "DEGENERATE CASE #2: Fov & t significant but Find not"
  } else {
    decision <- "PARTIAL EVIDENCE: check individual test results"
  }

  list(
    source = "bootstrap",
    Fov.pval = Fov_pval, t.pval = t_pval, Find.pval = Find_pval,
    F.cv = F_cv, t.cv = t_cv, Find.cv = I_cv,
    F.cv05.I0 = NA_real_, F.cv05.I1 = unname(F_cv["5%"]),
    t.cv05.I0 = NA_real_, t.cv05.I1 = unname(t_cv["5%"]),
    Find.cv05 = unname(I_cv["5%"]),
    nvalid = colSums(!is.na(stats_mat)),
    decision = decision,
    Fov.boot = Fb, t.boot = tb, Find.boot = Ib)
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
