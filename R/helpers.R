#' Build ARDL regression data (fixed q for all variables)
#' @noRd
.build_ardl_data <- function(y, X, p, q, fourier_sin = NULL, fourier_cos = NULL) {
  n <- length(y)
  nindep <- ncol(X)

  max_lag <- max(p, q)

  # Start index (need max_lag + 1 observations for differencing and lags)
  start_idx <- max_lag + 2

  if (start_idx > n) return(NULL)

  # Number of usable observations
  nobs <- n - start_idx + 1

  # Dependent variable: first difference of y
  dy <- diff(y)
  Y <- dy[start_idx:n - 1]

  # Build regressor matrix
  regressors <- list()
  coef_names <- c()

  # 1. Lagged levels (ECM terms)
  # L.y
  regressors$Ly <- y[(start_idx - 1):(n - 1)]
  coef_names <- c(coef_names, "L.y")

  # L.x for each independent variable
  for (j in 1:nindep) {
    regressors[[paste0("L.x", j)]] <- X[(start_idx - 1):(n - 1), j]
    coef_names <- c(coef_names, paste0("L.", colnames(X)[j]))
  }

  # 2. Lagged differences of y
  for (lag in 1:p) {
    regressors[[paste0("LD", lag, ".y")]] <- dy[(start_idx - lag - 1):(n - lag - 1)]
    coef_names <- c(coef_names, paste0("L", lag, ".D.y"))
  }

  # 3. Contemporaneous and lagged differences of X
  dX <- diff(X)
  for (j in 1:nindep) {
    for (lag in 0:q) {
      if (lag == 0) {
        regressors[[paste0("D.x", j)]] <- dX[(start_idx - 1):(n - 1), j]
        coef_names <- c(coef_names, paste0("D.", colnames(X)[j]))
      } else {
        regressors[[paste0("LD", lag, ".x", j)]] <- dX[(start_idx - lag - 1):(n - lag - 1), j]
        coef_names <- c(coef_names, paste0("L", lag, ".D.", colnames(X)[j]))
      }
    }
  }

  # 4. Fourier terms (if provided)
  if (!is.null(fourier_sin)) {
    regressors$sin <- fourier_sin[start_idx:n]
    coef_names <- c(coef_names, "sin")
  }
  if (!is.null(fourier_cos)) {
    regressors$cos <- fourier_cos[start_idx:n]
    coef_names <- c(coef_names, "cos")
  }

  # 5. Constant
  regressors$cons <- rep(1, length(Y))
  coef_names <- c(coef_names, "constant")

  # Combine into matrix
  Xmat <- do.call(cbind, regressors)

  return(list(Y = Y, Xmat = Xmat, coef_names = coef_names))
}


#' Build ARDL regression data (flexible q for each variable)
#' @noRd
.build_ardl_data_flex <- function(y, X, p, q_vec, fourier_sin = NULL, fourier_cos = NULL) {
  n <- length(y)
  nindep <- ncol(X)

  if (length(q_vec) != nindep) {
    stop("q_vec must have length equal to number of independent variables")
  }

  max_lag <- max(p, max(q_vec))

  # Start index
  start_idx <- max_lag + 2

  if (start_idx > n) return(NULL)

  # Dependent variable: first difference of y
  dy <- diff(y)
  Y <- dy[(start_idx - 1):(n - 1)]

  # Build regressor matrix
  regressors <- list()
  coef_names <- c()

  # 1. Lagged levels (ECM terms)
  regressors$Ly <- y[(start_idx - 1):(n - 1)]
  coef_names <- c(coef_names, "L.y")

  for (j in 1:nindep) {
    regressors[[paste0("L.x", j)]] <- X[(start_idx - 1):(n - 1), j]
    coef_names <- c(coef_names, paste0("L.", colnames(X)[j]))
  }

  # 2. Lagged differences of y
  for (lag in 1:p) {
    regressors[[paste0("LD", lag, ".y")]] <- dy[(start_idx - lag - 1):(n - lag - 1)]
    coef_names <- c(coef_names, paste0("L", lag, ".D.y"))
  }

  # 3. Contemporaneous and lagged differences of X (variable-specific q)
  dX <- diff(X)
  for (j in 1:nindep) {
    qj <- q_vec[j]
    for (lag in 0:qj) {
      if (lag == 0) {
        regressors[[paste0("D.x", j)]] <- dX[(start_idx - 1):(n - 1), j]
        coef_names <- c(coef_names, paste0("D.", colnames(X)[j]))
      } else {
        regressors[[paste0("LD", lag, ".x", j)]] <- dX[(start_idx - lag - 1):(n - lag - 1), j]
        coef_names <- c(coef_names, paste0("L", lag, ".D.", colnames(X)[j]))
      }
    }
  }

  # 4. Fourier terms
  if (!is.null(fourier_sin)) {
    regressors$sin <- fourier_sin[start_idx:n]
    coef_names <- c(coef_names, "sin")
  }
  if (!is.null(fourier_cos)) {
    regressors$cos <- fourier_cos[start_idx:n]
    coef_names <- c(coef_names, "cos")
  }

  # 5. Constant
  regressors$cons <- rep(1, length(Y))
  coef_names <- c(coef_names, "constant")

  Xmat <- do.call(cbind, regressors)
  colnames(Xmat) <- coef_names

  return(list(Y = Y, Xmat = Xmat, coef_names = coef_names))
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


#' PSS Bounds Test critical values
#' @noRd
.pss_bounds_test <- function(Fov_stat, t_stat, Find_stat, k, case, nobs) {
  # PSS (2001) asymptotic critical values (Case 3: unrestricted intercept, no trend)
  # These are approximate; Kripfganz & Schneider (2020) provide better finite-sample CVs

  # F-test critical values by k (number of independent variables)
  F_cvs <- list(
    "1" = list(I0 = c(4.04, 4.94, 6.84), I1 = c(4.78, 5.73, 7.84)),
    "2" = list(I0 = c(3.17, 3.79, 5.15), I1 = c(4.14, 4.85, 6.36)),
    "3" = list(I0 = c(2.72, 3.23, 4.29), I1 = c(3.77, 4.35, 5.61)),
    "4" = list(I0 = c(2.45, 2.86, 3.74), I1 = c(3.52, 4.01, 5.06)),
    "5" = list(I0 = c(2.26, 2.62, 3.41), I1 = c(3.35, 3.79, 4.68))
  )

  # t-test critical values
  t_cvs <- list(
    "1" = list(I0 = c(-2.57, -2.86, -3.43), I1 = c(-2.91, -3.22, -3.82)),
    "2" = list(I0 = c(-2.57, -2.86, -3.43), I1 = c(-3.21, -3.53, -4.10)),
    "3" = list(I0 = c(-2.57, -2.86, -3.43), I1 = c(-3.46, -3.78, -4.37)),
    "4" = list(I0 = c(-2.57, -2.86, -3.43), I1 = c(-3.66, -3.99, -4.60)),
    "5" = list(I0 = c(-2.57, -2.86, -3.43), I1 = c(-3.82, -4.16, -4.79))
  )

  k_str <- as.character(min(k, 5))

  F_I0 <- F_cvs[[k_str]]$I0  # 10%, 5%, 1%
  F_I1 <- F_cvs[[k_str]]$I1
  t_I0 <- t_cvs[[k_str]]$I0
  t_I1 <- t_cvs[[k_str]]$I1

  # F-test p-value (approximate using bounds)
  if (Fov_stat > F_I1[3]) {
    Fov_pval <- 0.005  # < 1%
  } else if (Fov_stat > F_I1[2]) {
    Fov_pval <- 0.025  # 1-5%
  } else if (Fov_stat > F_I1[1]) {
    Fov_pval <- 0.075  # 5-10%
  } else if (Fov_stat > F_I0[1]) {
    Fov_pval <- 0.15  # Inconclusive
  } else {
    Fov_pval <- 0.5  # > 10%
  }

  # t-test p-value
  if (t_stat < t_I1[3]) {
    t_pval <- 0.005
  } else if (t_stat < t_I1[2]) {
    t_pval <- 0.025
  } else if (t_stat < t_I1[1]) {
    t_pval <- 0.075
  } else if (t_stat < t_I0[1]) {
    t_pval <- 0.15
  } else {
    t_pval <- 0.5
  }

  # Find p-value (use F distribution approximation)
  Find_pval <- pf(Find_stat, k, nobs - k - 1, lower.tail = FALSE)

  # Decision
  Fov_reject <- Fov_stat > F_I1[2]  # Compare to 5% I(1) bound
  t_reject <- t_stat < t_I1[2]
  Find_reject <- Find_pval < 0.05

  if (Fov_reject && t_reject && Find_reject) {
    decision <- "COINTEGRATION detected (all tests significant at 5%)"
  } else if (!Fov_reject && !t_reject && !Find_reject) {
    decision <- "NO COINTEGRATION detected at 5% level"
  } else if (Fov_reject && Find_reject && !t_reject) {
    decision <- "DEGENERATE CASE #1: Fov & Find significant but t not (y may be I(0))"
  } else if (Fov_reject && t_reject && !Find_reject) {
    decision <- "DEGENERATE CASE #2: Fov & t significant but Find not (x not in ECM)"
  } else if (Fov_stat > F_I0[2] && Fov_stat < F_I1[2]) {
    decision <- "INCONCLUSIVE (F-statistic between bounds)"
  } else {
    decision <- "PARTIAL EVIDENCE: check individual test results"
  }

  return(list(
    Fov.pval = Fov_pval,
    t.pval = t_pval,
    Find.pval = Find_pval,
    F.cv05.I0 = F_I0[2],
    F.cv05.I1 = F_I1[2],
    t.cv05.I0 = t_I0[2],
    t.cv05.I1 = t_I1[2],
    Find.cv05 = qf(0.95, k, nobs - k - 1),
    decision = decision
  ))
}


#' Bootstrap ARDL test
#' @noRd
.bootstrap_ardl_test <- function(y, X, best_p, best_q, best_kstar, T,
                                  Fov_stat, t_stat, Find_stat,
                                  type, reps, final_data, coef_names,
                                  level_idx, indep_level_idx, ecm_idx) {

  n <- length(y)
  nindep <- ncol(X)

  # Generate Fourier terms
  ttrend <- 1:T
  if (best_kstar > 0) {
    fourier_sin <- sin(2 * pi * best_kstar * ttrend / T)
    fourier_cos <- cos(2 * pi * best_kstar * ttrend / T)
  } else {
    fourier_sin <- NULL
    fourier_cos <- NULL
  }

  # Original model residuals
  fit_original <- lm(final_data$Y ~ final_data$Xmat - 1)
  resid_orig <- fit_original$residuals
  coefs_orig <- coef(fit_original)

  # Store bootstrap test statistics
  Fov_boot <- numeric(reps)
  t_boot <- numeric(reps)
  Find_boot <- numeric(reps)

  # Bootstrap loop
  for (b in 1:reps) {
    # Resample residuals
    resid_boot <- sample(resid_orig, replace = TRUE)

    # Construct bootstrap dependent variable
    if (type == "fbardl_mcnown") {
      # Unconditional bootstrap: Y* = X*beta + e*
      Y_boot <- final_data$Xmat %*% coefs_orig + resid_boot
    } else {
      # Conditional bootstrap (Bertelli et al.): constrained null
      # Set ECM coefficient to 0 (no cointegration)
      coefs_null <- coefs_orig
      coefs_null[ecm_idx] <- 0
      Y_boot <- final_data$Xmat %*% coefs_null + resid_boot
    }

    # Re-estimate model
    fit_boot <- tryCatch({
      lm(Y_boot ~ final_data$Xmat - 1)
    }, error = function(e) NULL)

    if (is.null(fit_boot)) {
      Fov_boot[b] <- NA
      t_boot[b] <- NA
      Find_boot[b] <- NA
      next
    }

    coefs_b <- coef(fit_boot)
    vcov_b <- tryCatch(vcov(fit_boot), error = function(e) NULL)

    if (is.null(vcov_b) || any(is.na(diag(vcov_b)))) {
      Fov_boot[b] <- NA
      t_boot[b] <- NA
      Find_boot[b] <- NA
      next
    }

    se_b <- sqrt(diag(vcov_b))
    nparams <- length(coefs_b)
    df_resid <- length(Y_boot) - nparams

    # F-overall
    R_overall <- diag(nparams)[level_idx, , drop = FALSE]
    r_overall <- rep(0, length(level_idx))
    Fov_boot[b] <- tryCatch({
      .wald_f_test(coefs_b, vcov_b, R_overall, r_overall, df_resid)
    }, error = function(e) NA)

    # t-statistic
    t_boot[b] <- coefs_b[ecm_idx] / se_b[ecm_idx]

    # F-independent
    R_ind <- diag(nparams)[indep_level_idx, , drop = FALSE]
    r_ind <- rep(0, length(indep_level_idx))
    Find_boot[b] <- tryCatch({
      .wald_f_test(coefs_b, vcov_b, R_ind, r_ind, df_resid)
    }, error = function(e) NA)
  }

  # Remove NAs
  Fov_boot <- Fov_boot[!is.na(Fov_boot)]
  t_boot <- t_boot[!is.na(t_boot)]
  Find_boot <- Find_boot[!is.na(Find_boot)]

  # Compute critical values and p-values
  Fov_cv <- quantile(Fov_boot, c(0.90, 0.95, 0.975, 0.99), na.rm = TRUE)
  t_cv <- quantile(t_boot, c(0.01, 0.025, 0.05, 0.10), na.rm = TRUE)
  Find_cv <- quantile(Find_boot, c(0.90, 0.95, 0.975, 0.99), na.rm = TRUE)

  Fov_pval <- mean(Fov_boot >= Fov_stat, na.rm = TRUE)
  t_pval <- mean(t_boot <= t_stat, na.rm = TRUE)
  Find_pval <- mean(Find_boot >= Find_stat, na.rm = TRUE)

  # Decision
  if (Fov_pval < 0.05 && t_pval < 0.05 && Find_pval < 0.05) {
    decision <- "COINTEGRATION detected (all tests significant at 5%)"
  } else if (Fov_pval >= 0.05 && t_pval >= 0.05 && Find_pval >= 0.05) {
    decision <- "NO COINTEGRATION detected at 5% level"
  } else if (Fov_pval < 0.05 && Find_pval < 0.05 && t_pval >= 0.05) {
    decision <- "DEGENERATE CASE #1: Fov & Find significant but t not (y may be I(0))"
  } else if (Fov_pval < 0.05 && t_pval < 0.05 && Find_pval >= 0.05) {
    decision <- "DEGENERATE CASE #2: Fov & t significant but Find not (x not in ECM)"
  } else {
    decision <- "PARTIAL EVIDENCE: check individual test results"
  }

  return(list(
    Fov.pval = Fov_pval,
    t.pval = t_pval,
    Find.pval = Find_pval,
    F.cv05.I0 = Fov_cv["90%"],
    F.cv05.I1 = Fov_cv["95%"],
    F.cv01 = Fov_cv["99%"],
    t.cv05.I0 = t_cv["10%"],
    t.cv05.I1 = t_cv["5%"],
    t.cv01 = t_cv["1%"],
    Find.cv05 = Find_cv["95%"],
    Find.cv01 = Find_cv["99%"],
    decision = decision,
    Fov.boot = Fov_boot,
    t.boot = t_boot,
    Find.boot = Find_boot
  ))
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
.run_diagnostics <- function(residuals, nobs, nparams) {
  n <- length(residuals)

  # Jarque-Bera normality test
  m3 <- mean((residuals - mean(residuals))^3)
  m2 <- mean((residuals - mean(residuals))^2)
  skew <- m3 / m2^1.5
  m4 <- mean((residuals - mean(residuals))^4)
  kurt <- m4 / m2^2 - 3
  jb_stat <- n * (skew^2 / 6 + kurt^2 / 24)
  jb_pval <- pchisq(jb_stat, 2, lower.tail = FALSE)

  # Breusch-Godfrey serial correlation test (AR1)
  bg1 <- tryCatch({
    resid_lag1 <- c(NA, residuals[-n])
    fit_bg1 <- lm(residuals ~ resid_lag1)
    r2_bg1 <- summary(fit_bg1)$r.squared
    stat <- n * r2_bg1
    pval <- pchisq(stat, 1, lower.tail = FALSE)
    list(stat = stat, pval = pval)
  }, error = function(e) list(stat = NA, pval = NA))

  # Breusch-Godfrey AR(4)
  bg4 <- tryCatch({
    if (n > 8) {
      resid_lags <- embed(c(rep(NA, 4), residuals), 5)
      fit_bg4 <- lm(resid_lags[, 1] ~ resid_lags[, 2:5])
      r2_bg4 <- summary(fit_bg4)$r.squared
      stat <- (n - 4) * r2_bg4
      pval <- pchisq(stat, 4, lower.tail = FALSE)
      list(stat = stat, pval = pval)
    } else {
      list(stat = NA, pval = NA)
    }
  }, error = function(e) list(stat = NA, pval = NA))

  # Breusch-Pagan heteroskedasticity test
  bp <- tryCatch({
    resid_sq <- residuals^2
    fit_bp <- lm(resid_sq ~ seq_along(residuals))
    ess <- sum((fitted(fit_bp) - mean(resid_sq))^2)
    tss <- sum((resid_sq - mean(resid_sq))^2)
    r2_bp <- ess / tss
    stat <- n * r2_bp
    pval <- pchisq(stat, 1, lower.tail = FALSE)
    list(stat = stat, pval = pval)
  }, error = function(e) list(stat = NA, pval = NA))

  # ARCH(1) test
  arch1 <- tryCatch({
    resid_sq <- residuals^2
    resid_sq_lag <- c(NA, resid_sq[-n])
    fit_arch <- lm(resid_sq ~ resid_sq_lag)
    r2_arch <- summary(fit_arch)$r.squared
    stat <- n * r2_arch
    pval <- pchisq(stat, 1, lower.tail = FALSE)
    list(stat = stat, pval = pval)
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
