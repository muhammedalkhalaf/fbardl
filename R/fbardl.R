#' Fourier Bootstrap ARDL Cointegration Test
#'
#' @description
#' Performs the Fourier Bootstrap ARDL (FBARDL) bounds testing approach for
#' cointegration analysis. This function combines the Pesaran, Shin & Smith
#' (2001) ARDL bounds testing framework with Fourier terms to capture structural
#' breaks, and provides bootstrap critical values for robust inference.
#'
#' @param formula A formula of the form \code{y ~ x1 + x2 + ...} specifying the
#'   dependent and independent variables.
#' @param data A data frame containing the time series variables.
#' @param type Character string specifying the test type:
#'   \itemize{
#'     \item \code{"fardl"}: Fourier ARDL with PSS bounds test (default)
#'     \item \code{"fbardl_mcnown"}: Bootstrap ARDL (McNown, Sam & Goh, 2018)
#'     \item \code{"fbardl_bvz"}: Bootstrap ARDL (Bertelli, Vacca & Zoia, 2022)
#'   }
#' @param maxlag Integer. Maximum lag order for grid search (default: 4).
#' @param maxk Numeric. Maximum Fourier frequency (default: 5).
#' @param ic Character string. Information criterion for lag selection:
#'   \code{"aic"} (default) or \code{"bic"}.
#' @param case Integer. PSS case specification (2, 3, 4, or 5). Default is 3
#'   (unrestricted intercept, no trend).
#' @param reps Integer. Number of bootstrap replications (default: 999).
#' @param fourier Logical. Whether to include Fourier terms (default: TRUE).
#' @param level Numeric. Confidence level for intervals (default: 0.95).
#' @param horizon Integer. Horizon for dynamic multipliers (default: 20).
#' @param unconditional Logical. If \code{TRUE}, the contemporaneous
#'   differences of the regressors are excluded (the unconditional form of
#'   Yilanci et al., 2020, and McNown et al., 2018). Default \code{FALSE}.
#'
#' @return An object of class \code{"fbardl"} containing:
#' \describe{
#'   \item{coefficients}{Named vector of estimated coefficients}
#'   \item{std.errors}{Standard errors of coefficients}
#'   \item{t.values}{t-statistics}
#'   \item{p.values}{p-values}
#'   \item{long.run}{Long-run coefficient estimates with standard errors}
#'   \item{short.run}{Short-run coefficient estimates}
#'   \item{ecm.coef}{Error correction coefficient (speed of adjustment)}
#'   \item{best.p}{Selected lag order for dependent variable}
#'   \item{best.q}{Selected lag orders for independent variables}
#'   \item{best.kstar}{Selected Fourier frequency}
#'   \item{F.overall}{F-statistic for overall cointegration test}
#'   \item{t.dependent}{t-statistic on lagged dependent variable}
#'   \item{F.independent}{F-statistic on lagged independent variables}
#'   \item{cointegration}{Cointegration test results with critical values}
#'   \item{diagnostics}{Diagnostic test results}
#'   \item{model.fit}{Model fit statistics (R2, AIC, BIC, etc.)}
#'   \item{residuals}{Model residuals}
#'   \item{fitted.values}{Fitted values}
#'   \item{nobs}{Number of observations}
#'   \item{call}{The matched call}
#' }
#'
#' @details
#' The FBARDL approach extends the standard ARDL bounds testing procedure by:
#' \enumerate{
#'   \item Incorporating Fourier terms to capture smooth structural breaks
#'   \item Using bootstrap methods to generate finite-sample critical values
#'   \item Implementing the McNown et al. (2018) procedure to detect degenerate cases
#' }
#'
#' The PSS case sets the deterministic terms: a constant in all cases and a
#' linear trend in cases 4 and 5. The overall F test restricts the lagged
#' levels, together with the intercept in case 2 and the trend in case 4.
#'
#' With \code{type = "fardl"} the F and t statistics are compared with the
#' finite-sample critical values and approximate p-values of Kripfganz and
#' Schneider (2020), computed from their response surface coefficients
#' (Stata package \code{ardl}) for the sample size, the number of
#' regressors and the number of short-run coefficients (the Fourier terms
#' included). These bounds do not account for the Fourier terms; the
#' F test on the lagged regressors has no tabulated distribution and is
#' reported with the bootstrap types only.
#'
#' The bootstrap types follow the Stata module \code{fbardl} 1.3.0. Data
#' are generated recursively under the null from the restricted
#' equilibrium-correction equation for \eqn{y} and the equations for
#' \eqn{\Delta x}, resampling the residual pairs. \code{"fbardl_mcnown"}
#' uses one null (all lagged levels, plus the restricted deterministic term)
#' for the three statistics and an unrestricted \eqn{\Delta x} equation;
#' \code{"fbardl_bvz"} uses a separate null for each statistic, a marginal
#' \eqn{\Delta x} equation without the lagged level of \eqn{y}, and
#' recentred residuals.
#'
#' The procedure involves three main steps:
#' \enumerate{
#'   \item Selection of optimal Fourier frequency k* by minimum SSR
#'   \item Selection of lag orders (p, q) by AIC or BIC
#'   \item Cointegration testing with bootstrap or PSS critical values
#' }
#'
#' Three test statistics are computed:
#' \itemize{
#'   \item \code{F.overall}: Joint test on all lagged level variables
#'   \item \code{t.dependent}: t-test on lagged dependent variable
#'   \item \code{F.independent}: Joint test on lagged independent variables
#' }
#'
#' @references
#' Kripfganz, S. and Schneider, D. C. (2020). Response surface regressions
#' for critical value bounds and approximate p-values in equilibrium
#' correction models. \emph{Oxford Bulletin of Economics and Statistics},
#' 82(6), 1456-1481. \doi{10.1111/obes.12377}
#'
#' Pesaran, M. H., Shin, Y., & Smith, R. J. (2001). Bounds testing approaches
#' to the analysis of level relationships. \emph{Journal of Applied Econometrics},
#' 16(3), 289-326. \doi{10.1002/jae.616}
#'
#' McNown, R., Sam, C. Y., & Goh, S. K. (2018). Bootstrapping the autoregressive
#' distributed lag test for cointegration. \emph{Applied Economics}, 50(13),
#' 1509-1521. \doi{10.1080/00036846.2017.1366643}
#'
#' Yilanci, V., Bozoklu, S., & Gorus, M. S. (2020). Are BRICS countries
#' pollution havens? Evidence from a bootstrap ARDL bounds testing approach
#' with a Fourier function. \emph{Sustainable Cities and Society}, 55, 102035.
#' \doi{10.1016/j.scs.2020.102035}
#'
#' Kripfganz, S., & Schneider, D. C. (
#' 2020). Response surface regressions for critical value bounds and approximate
#' p-values in equilibrium correction models. \emph{Oxford Bulletin of Economics
#' and Statistics}, 82(6), 1456-1481. \doi{10.1111/obes.12377}
#'
#' @examples
#' \donttest{
#' # Load example data
#' data(fbardl_data)
#'
#' # Basic Fourier ARDL test
#' result <- fbardl(y ~ x1 + x2, data = fbardl_data, type = "fardl")
#' summary(result)
#'
#' # Bootstrap ARDL (McNown approach)
#' result_boot <- fbardl(y ~ x1 + x2, data = fbardl_data,
#'                       type = "fbardl_mcnown", reps = 499)
#' summary(result_boot)
#'
#' # Without Fourier terms
#' result_nofourier <- fbardl(y ~ x1 + x2, data = fbardl_data,
#'                            fourier = FALSE)
#' }
#'
#' @export
fbardl <- function(formula, data, type = c("fardl", "fbardl_mcnown", "fbardl_bvz"),
                   maxlag = 4, maxk = 5, ic = c("aic", "bic"), case = 3,
                   reps = 999, fourier = TRUE, level = 0.95, horizon = 20,
                   unconditional = FALSE) {

  # Match arguments
  type <- match.arg(type)
  ic <- match.arg(ic)

  # Validate inputs
  if (!inherits(formula, "formula")) {
    stop("'formula' must be a formula object")
  }
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }
  if (maxlag < 1 || maxlag > 12) {
    stop("'maxlag' must be between 1 and 12")
  }
  if (!case %in% c(2, 3, 4, 5)) {
    stop("'case' must be 2, 3, 4, or 5")
  }
  hastrend <- case %in% c(4, 5)
  j0 <- if (isTRUE(unconditional)) 1L else 0L

  # Extract variables from formula
  mf <- model.frame(formula, data = data, na.action = na.pass)
  y <- model.response(mf)
  X <- model.matrix(formula, data = mf)[, -1, drop = FALSE]  # Remove intercept

  depvar <- all.vars(formula)[1]
  indepvars <- colnames(X)
  nindep <- ncol(X)

  if (nindep < 1) {
    stop("At least one independent variable required")
  }

  # Convert to matrices, handle NAs
  y <- as.numeric(y)
  X <- as.matrix(X)
  n_original <- length(y)

  if (n_original < 20) {
    stop("Sample size too small: need at least 20 observations")
  }

  # Store call
  cl <- match.call()

  # ============================================================================
  # STEP 1: SELECT FOURIER FREQUENCY k* BY MINIMUM SSR
  # ============================================================================
  message("Step 1: Selecting Fourier frequency k* by minimum SSR...")

  if (fourier) {
    kgrid <- seq(0.1, maxk, by = 0.1)
  } else {
    kgrid <- 0
  }

  best_kstar <- 0
  best_ssr_k <- Inf
  ssr_by_k <- numeric(length(kgrid))

  T <- n_original  # For Fourier terms
  ttrend <- 1:T

  for (i in seq_along(kgrid)) {
    kval <- kgrid[i]

    # Generate Fourier terms
    if (kval > 0) {
      fourier_sin <- sin(2 * pi * kval * ttrend / T)
      fourier_cos <- cos(2 * pi * kval * ttrend / T)
    } else {
      fourier_sin <- NULL
      fourier_cos <- NULL
    }

    # Build regression data for maximum lag
    reg_data <- .build_ardl_data_flex(y, X, maxlag, rep(maxlag, nindep),
                                      fourier_sin, fourier_cos, hastrend, j0)

    if (is.null(reg_data) || length(reg_data$Y) < 10) {
      ssr_by_k[i] <- Inf
      next
    }

    # Estimate OLS
    fit <- tryCatch({
      lm.fit(reg_data$Xmat, reg_data$Y)
    }, error = function(e) NULL)

    if (!is.null(fit)) {
      ssr_by_k[i] <- sum(fit$residuals^2)
      if (ssr_by_k[i] < best_ssr_k) {
        best_ssr_k <- ssr_by_k[i]
        best_kstar <- kval
      }
    } else {
      ssr_by_k[i] <- Inf
    }
  }

  if (!fourier) best_kstar <- 0

  message(sprintf("  Optimal k* = %.1f (min SSR = %.4f)", best_kstar, best_ssr_k))

  # ============================================================================
  # STEP 2: SELECT LAG ORDERS (p, q) BY AIC/BIC
  # ============================================================================
  message(sprintf("Step 2: Selecting lag orders (p, q) by %s...", toupper(ic)))

  # Generate Fourier terms at optimal k*
  if (best_kstar > 0) {
    fourier_sin <- sin(2 * pi * best_kstar * ttrend / T)
    fourier_cos <- cos(2 * pi * best_kstar * ttrend / T)
  } else {
    fourier_sin <- NULL
    fourier_cos <- NULL
  }

  best_ic_val <- Inf
  best_p <- 1
  best_q <- rep(0, nindep)
  names(best_q) <- indepvars
  total_specs <- 0

  # Grid search over p and all combinations of q
  nq <- maxlag + 1  # q can be 0, 1, ..., maxlag

  for (p in 1:maxlag) {
    # Generate all combinations of q for each independent variable
    q_combinations <- expand.grid(replicate(nindep, 0:maxlag, simplify = FALSE))
    colnames(q_combinations) <- indepvars

    for (qidx in 1:nrow(q_combinations)) {
      total_specs <- total_specs + 1
      q_vec <- as.integer(q_combinations[qidx, ])

      # Build regression data
      reg_data <- .build_ardl_data_flex(y, X, p, q_vec, fourier_sin, fourier_cos,
                                        hastrend, j0)

      if (is.null(reg_data) || length(reg_data$Y) < 10) next

      # Estimate OLS
      fit <- tryCatch({
        lm.fit(reg_data$Xmat, reg_data$Y)
      }, error = function(e) NULL)

      if (is.null(fit)) next

      nobs_tmp <- length(fit$residuals)
      k_tmp <- ncol(reg_data$Xmat)
      rss_tmp <- sum(fit$residuals^2)

      # Log-likelihood (assuming Gaussian errors)
      ll_tmp <- -nobs_tmp/2 * (log(2 * pi) + log(rss_tmp/nobs_tmp) + 1)

      # Information criterion
      if (ic == "aic") {
        ic_tmp <- -2 * ll_tmp + 2 * k_tmp
      } else {
        ic_tmp <- -2 * ll_tmp + k_tmp * log(nobs_tmp)
      }

      if (ic_tmp < best_ic_val) {
        best_ic_val <- ic_tmp
        best_p <- p
        best_q <- q_vec
        names(best_q) <- indepvars
      }
    }
  }

  message(sprintf("  Optimal p = %d, q = (%s) (%d models evaluated)",
                  best_p, paste(best_q, collapse = ", "), total_specs))

  # ============================================================================
  # STEP 3: FINAL ESTIMATION
  # ============================================================================
  message(sprintf("Step 3: Estimating final ARDL(%d,%s) model with k* = %.1f...",
                  best_p, paste(best_q, collapse = ","), best_kstar))

  # Build final regression data
  final_data <- .build_ardl_data_flex(y, X, best_p, best_q, fourier_sin, fourier_cos,
                                      hastrend, j0)

  # Final OLS estimation
  final_fit <- lm(final_data$Y ~ final_data$Xmat - 1)

  nobs <- length(final_fit$residuals)
  nparams <- length(coef(final_fit))

  # Extract coefficients and statistics
  coefs <- coef(final_fit)
  vcov_mat <- vcov(final_fit)
  se <- sqrt(diag(vcov_mat))
  tvals <- coefs / se
  pvals <- 2 * pt(abs(tvals), df = nobs - nparams, lower.tail = FALSE)

  # Model fit statistics
  # The design holds its own constant, so lm() is called without one and
  # summary() would report the uncentred R-squared
  tss <- sum((final_data$Y - mean(final_data$Y))^2)
  r2 <- 1 - sum(final_fit$residuals^2) / tss
  r2_adj <- 1 - (1 - r2) * (nobs - 1) / (nobs - nparams)
  rss <- sum(final_fit$residuals^2)
  ll <- -nobs/2 * (log(2 * pi) + log(rss/nobs) + 1)
  aic_val <- -2 * ll + 2 * nparams
  bic_val <- -2 * ll + nparams * log(nobs)
  fstat <- summary(final_fit)$fstatistic
  rmse <- sqrt(rss / (nobs - nparams))

  # ============================================================================
  # EXTRACT ECM AND COMPUTE TEST STATISTICS
  # ============================================================================

  # Get coefficient names
  coef_names <- final_data$coef_names

  # ECM coefficient (coefficient on lagged dependent variable)
  ecm_idx <- which(coef_names == "L.y")
  ecm_coef <- unname(coefs[ecm_idx])
  ecm_se <- unname(se[ecm_idx])
  ecm_t <- unname(tvals[ecm_idx])
  ecm_p <- unname(pvals[ecm_idx])

  # Indices for level variables (for F-tests)
  level_idx <- grep("^L\\.", coef_names)
  indep_level_idx <- setdiff(level_idx, ecm_idx)

  # F-test: overall (all lagged levels)
  # Compute F-statistic manually
  # The F_ov restriction depends on the PSS case: the intercept joins it in
  # case 2 and the trend in case 4
  fov_names <- coef_names[level_idx]
  if (case == 2) fov_names <- c(fov_names, "constant")
  if (case == 4) fov_names <- c(fov_names, "trend")
  fov_idx <- match(fov_names, coef_names)
  R_overall <- diag(nparams)[fov_idx, , drop = FALSE]
  r_overall <- rep(0, length(fov_idx))
  Fov_stat <- .wald_f_test(coefs, vcov_mat, R_overall, r_overall, nobs - nparams)

  # F-test: independent variables only
  R_ind <- diag(nparams)[indep_level_idx, , drop = FALSE]
  r_ind <- rep(0, length(indep_level_idx))
  Find_stat <- .wald_f_test(coefs, vcov_mat, R_ind, r_ind, nobs - nparams)

  # t-statistic on dependent variable
  t_stat <- ecm_t

  # ============================================================================
  # COINTEGRATION TESTING
  # ============================================================================

  if (type == "fardl") {
    # Short-run coefficients: every regressor other than the constant, the
    # trend and the k + 1 lagged levels (Fourier terms included)
    sr <- nparams - 1L - (nindep + 1L) - as.integer(hastrend)
    coint_result <- .pss_bounds_test(Fov_stat, t_stat, nindep, case, nobs, sr)
  } else {
    coint_result <- .bootstrap_ardl_test(
      y, X, best_p, best_q, fourier_sin, fourier_cos, hastrend, j0,
      case, type, reps, Fov_stat, t_stat, Find_stat, fov_names)
  }

  # ============================================================================
  # LONG-RUN COEFFICIENTS
  # ============================================================================

  long_run <- .compute_long_run(coefs, vcov_mat, coef_names, indepvars,
                                 ecm_idx, indep_level_idx, nobs - nparams, level)

  # ============================================================================
  # SHORT-RUN COEFFICIENTS
  # ============================================================================

  short_run <- .extract_short_run(coefs, se, tvals, pvals, coef_names,
                                   depvar, indepvars, best_p, best_q)

  # ============================================================================
  # DIAGNOSTIC TESTS
  # ============================================================================

  diagnostics <- .run_diagnostics(final_fit$residuals, nobs, nparams,
                                 final_data$Xmat, final_fit$fitted.values)

  # ============================================================================
  # ASSEMBLE RESULTS
  # ============================================================================

  result <- list(
    coefficients = coefs,
    std.errors = se,
    t.values = tvals,
    p.values = pvals,
    coef.names = coef_names,
    vcov = vcov_mat,
    long.run = long_run,
    short.run = short_run,
    ecm.coef = ecm_coef,
    ecm.se = ecm_se,
    ecm.t = ecm_t,
    ecm.p = ecm_p,
    best.p = best_p,
    best.q = best_q,
    best.kstar = best_kstar,
    F.overall = Fov_stat,
    t.dependent = t_stat,
    F.independent = Find_stat,
    cointegration = coint_result,
    diagnostics = diagnostics,
    model.fit = list(
      r.squared = r2,
      adj.r.squared = r2_adj,
      log.likelihood = ll,
      aic = aic_val,
      bic = bic_val,
      F.statistic = if (!is.null(fstat)) fstat[1] else NA,
      rmse = rmse,
      ic.val = best_ic_val
    ),
    residuals = final_fit$residuals,
    fitted.values = final_fit$fitted.values,
    nobs = nobs,
    nparams = nparams,
    depvar = depvar,
    indepvars = indepvars,
    type = type,
    ic = ic,
    case = case,
    unconditional = isTRUE(unconditional),
    fourier = fourier,
    level = level,
    reps = if (type != "fardl") reps else NULL,
    total.specs = total_specs,
    call = cl
  )

  class(result) <- "fbardl"
  return(result)
}


#' @export
print.fbardl <- function(x, ...) {
  cat("\nFourier Bootstrap ARDL Cointegration Test\n")
  cat(rep("=", 60), "\n", sep = "")

  cat("\nCall:\n")
  print(x$call)

  cat("\nModel: ARDL(", x$best.p, ",",
      paste(x$best.q, collapse = ","), ")\n", sep = "")
  cat("Fourier frequency (k*): ", x$best.kstar, "\n", sep = "")
  cat("Test type: ", x$type, "\n", sep = "")
  cat("Sample size: ", x$nobs, "\n", sep = "")

  cat("\nCointegration Test Results:\n")
  cat(rep("-", 50), "\n", sep = "")
  co <- x$cointegration
  pv <- function(v) if (length(v) == 2) sprintf("p I(0)/I(1): %6.4f / %6.4f", v[1], v[2]) else
    if (is.na(v)) "p-value: n/a" else sprintf("p-value: %6.4f", v)
  cat(sprintf("F-overall:     %8.4f  (%s)\n", x$F.overall, pv(co$Fov.pval)))
  cat(sprintf("t-dependent:   %8.4f  (%s)\n", x$t.dependent, pv(co$t.pval)))
  cat(sprintf("F-independent: %8.4f  (%s)\n", x$F.independent, pv(co$Find.pval)))
  cat("Critical values: ", co$source, "\n", sep = "")
  cat("\nDecision: ", x$cointegration$decision, "\n", sep = "")

  cat("\nError Correction Coefficient: ", sprintf("%.6f", x$ecm.coef), "\n", sep = "")

  cat("\nModel Fit:\n")
  cat(sprintf("R-squared: %.4f, Adj. R-squared: %.4f\n",
              x$model.fit$r.squared, x$model.fit$adj.r.squared))
  cat(sprintf("AIC: %.4f, BIC: %.4f\n",
              x$model.fit$aic, x$model.fit$bic))

  invisible(x)
}


#' @export
summary.fbardl <- function(object, ...) {
  cat("\n")
  cat(rep("=", 78), "\n", sep = "")

  # Header based on type
  if (object$type == "fardl") {
    cat("  Fourier ARDL (FARDL) Cointegration Analysis\n")
  } else if (object$type == "fbardl_mcnown") {
    cat("  Fourier Bootstrap ARDL - McNown, Sam & Goh (2018)\n")
  } else {
    cat("  Fourier Bootstrap ARDL - Bertelli, Vacca & Zoia (2022)\n")
  }

  cat(rep("=", 78), "\n", sep = "")
  cat("  Dependent variable  : ", object$depvar, "\n", sep = "")
  cat("  Independent var(s)  : ", paste(object$indepvars, collapse = ", "), "\n", sep = "")
  cat("  Sample size (T)     : ", object$nobs, "\n", sep = "")
  cat("  Max lag order       : ", object$best.p, "\n", sep = "")
  if (object$fourier) {
    cat("  Fourier frequency   : ", object$best.kstar, "\n", sep = "")
  } else {
    cat("  Fourier terms       : excluded\n")
  }
  cat("  Information crit.   : ", toupper(object$ic), "\n", sep = "")
  cat("  PSS Case            : Case ", object$case, "\n", sep = "")
  if (object$type != "fardl") {
    cat("  Bootstrap reps      : ", object$reps, "\n", sep = "")
  }
  cat(rep("=", 78), "\n", sep = "")

  # Table 1: Model Selection Summary
  cat("\nTable 1: Model Selection Summary\n")
  cat(rep("-", 68), "\n", sep = "")
  cat(sprintf("  ARDL Specification            : ARDL(%d,%s)\n",
              object$best.p, paste(object$best.q, collapse = ",")))
  cat(sprintf("  Fourier Frequency (k*)        : %.1f\n", object$best.kstar))
  cat(sprintf("  Observations                  : %d\n", object$nobs))
  cat(sprintf("  R-squared                     : %.6f\n", object$model.fit$r.squared))
  cat(sprintf("  Adjusted R-squared            : %.6f\n", object$model.fit$adj.r.squared))
  cat(sprintf("  Log-Likelihood                : %.4f\n", object$model.fit$log.likelihood))
  cat(sprintf("  AIC                           : %.4f\n", object$model.fit$aic))
  cat(sprintf("  BIC                           : %.4f\n", object$model.fit$bic))
  cat(sprintf("  RMSE                          : %.6f\n", object$model.fit$rmse))
  cat(sprintf("  Models evaluated              : %d\n", object$total.specs))
  cat(rep("-", 68), "\n", sep = "")

  # Table 2: Coefficient Estimates
  cat("\nTable 2: ARDL Regression Results (EC Representation)\n")
  cat(rep("-", 68), "\n", sep = "")

  # Speed of adjustment
  cat("\n  ADJ - Speed of Adjustment:\n")
  cat(sprintf("    L.%s: %.6f (SE: %.6f, t: %.4f, p: %.4f)%s\n",
              object$depvar, object$ecm.coef, object$ecm.se,
              object$ecm.t, object$ecm.p, .stars(object$ecm.p)))

  # Long-run coefficients
  cat("\n  LR - Long-Run Coefficients:\n")
  for (i in seq_len(nrow(object$long.run))) {
    cat(sprintf("    %s: %.6f (SE: %.6f, z: %.4f, p: %.4f)%s\n",
                rownames(object$long.run)[i],
                object$long.run[i, "estimate"],
                object$long.run[i, "std.error"],
                object$long.run[i, "z.value"],
                object$long.run[i, "p.value"],
                .stars(object$long.run[i, "p.value"])))
  }

  # Short-run coefficients
  cat("\n  SR - Short-Run Coefficients:\n")
  for (i in seq_len(nrow(object$short.run))) {
    cat(sprintf("    %s: %.6f (SE: %.6f, t: %.4f, p: %.4f)%s\n",
                rownames(object$short.run)[i],
                object$short.run[i, "estimate"],
                object$short.run[i, "std.error"],
                object$short.run[i, "t.value"],
                object$short.run[i, "p.value"],
                .stars(object$short.run[i, "p.value"])))
  }

  cat("\n  Stars: *** p<0.01, ** p<0.05, * p<0.10\n")
  cat(rep("-", 68), "\n", sep = "")

  # Table 3: Cointegration Test Results
  cat("\nTable 3: Cointegration Test Results\n")
  cat(rep("-", 68), "\n", sep = "")

  coint <- object$cointegration
  cat("  Critical values: ", coint$source, "\n", sep = "")
  if (coint$source == "bootstrap") {
    row <- function(lab, stat, cv, pval) {
      cat(sprintf("  %-14s %9.4f  CV 10%%/5%%/1%%: %8.3f %8.3f %8.3f  p = %.4f%s\n",
                  lab, stat, cv["10%"], cv["5%"], cv["1%"], pval, .stars(pval)))
    }
    row("F-overall:", object$F.overall, coint$F.cv, coint$Fov.pval)
    row("t-dependent:", object$t.dependent, coint$t.cv, coint$t.pval)
    row("F-independent:", object$F.independent, coint$Find.cv, coint$Find.pval)
  } else {
    row <- function(lab, stat, cv, pval) {
      cat(sprintf("  %-14s %9.4f  5%% bounds [%.3f, %.3f]  1%% bounds [%.3f, %.3f]  p I(0)/I(1) = %.4f / %.4f\n",
                  lab, stat, cv["I0", "5%"], cv["I1", "5%"], cv["I0", "1%"],
                  cv["I1", "1%"], pval[1], pval[2]))
    }
    row("F-overall:", object$F.overall, coint$F.cv, coint$Fov.pval)
    row("t-dependent:", object$t.dependent, coint$t.cv, coint$t.pval)
    cat(sprintf("  %-14s %9.4f  (no tabulated critical values; use a bootstrap type)\n",
                "F-independent:", object$F.independent))
  }
  cat("\n  Decision: ", coint$decision, "\n", sep = "")
  cat(rep("-", 68), "\n", sep = "")

  # Table 4: Diagnostic Tests
  cat("\nTable 4: Diagnostic Tests\n")
  cat(rep("-", 68), "\n", sep = "")
  diag <- object$diagnostics
  cat(sprintf("  Jarque-Bera (Normality)      : %.4f  (p = %.4f)%s\n",
              diag$jb.stat, diag$jb.pval, .stars(1 - diag$jb.pval)))
  cat(sprintf("  Breusch-Godfrey (AR 1)       : %.4f  (p = %.4f)%s\n",
              diag$bg1.stat, diag$bg1.pval, .stars(1 - diag$bg1.pval)))
  cat(sprintf("  Breusch-Godfrey (AR 4)       : %.4f  (p = %.4f)%s\n",
              diag$bg4.stat, diag$bg4.pval, .stars(1 - diag$bg4.pval)))
  cat(sprintf("  Breusch-Pagan (Heterosked.)  : %.4f  (p = %.4f)%s\n",
              diag$bp.stat, diag$bp.pval, .stars(1 - diag$bp.pval)))
  cat(sprintf("  ARCH(1) Test                 : %.4f  (p = %.4f)%s\n",
              diag$arch1.stat, diag$arch1.pval, .stars(1 - diag$arch1.pval)))
  cat(rep("-", 68), "\n", sep = "")

  cat("\n")
  cat(rep("=", 78), "\n", sep = "")
  cat("  fbardl v1.0.0\n")
  cat(rep("=", 78), "\n", sep = "")

  invisible(object)
}


#' Significance stars helper
#' @noRd
.stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.01) return(" ***")
  if (p < 0.05) return(" **")
  if (p < 0.10) return(" *")
  return("")
}
