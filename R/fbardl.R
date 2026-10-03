#' Fourier Bootstrap ARDL Cointegration Test
#'
#' @description
#' ARDL bounds tests for cointegration (Pesaran, Shin and Smith, 2001) with
#' optional Fourier terms for smooth breaks and bootstrap critical values
#' (McNown, Sam and Goh, 2018; Bertelli, Vacca and Zoia, 2022).
#'
#' @param formula A formula of the form \code{y ~ x1 + x2 + ...} specifying the
#'   dependent and independent variables.
#' @param data A data frame containing the time series variables.
#' @param type Character string specifying the test type:
#'   \itemize{
#'     \item \code{"fbardl_bvz"} (default): bootstrap with separate nulls
#'       (Bertelli, Vacca and Zoia, 2022)
#'     \item \code{"fbardl_mcnown"}: bootstrap with the null of the overall
#'       F test for all statistics (McNown, Sam and Goh, 2018)
#'     \item \code{"fardl"}: Kripfganz and Schneider (2020) bounds; valid only
#'       without Fourier terms (see Details)
#'   }
#' @param maxlag Integer. Maximum lag order for grid search (default: 4).
#' @param maxk Numeric. Maximum Fourier frequency (default: 5).
#' @param ic Character string. Information criterion for lag selection:
#'   \code{"aic"} (default) or \code{"bic"}.
#' @param case Integer. PSS case specification (2, 3, 4, or 5). Default is 3
#'   (unrestricted intercept, no trend).
#' @param reps Integer. Number of bootstrap replications (default: 999).
#' @param fourier Logical. Whether to include Fourier terms (default: TRUE).
#' @param level Numeric. Confidence level of the long-run coefficient
#'   intervals (default: 0.95). Test decisions are made at the 5\% level.
#' @param horizon Not used; kept for compatibility with earlier versions
#'   (no dynamic multipliers are computed).
#' @param unconditional Logical. If \code{TRUE}, the contemporaneous
#'   differences of the regressors are excluded from the estimated model, as
#'   in the simulation model of McNown, Sam and Goh (2018, eq. 12). Default
#'   \code{FALSE}.
#' @param kgrid Grid of candidate Fourier frequencies when k* is selected:
#'   \code{"integer"} (default), \eqn{k = 1, \dots,} \code{maxk}, or
#'   \code{"fractional"}, \eqn{k = 0.1, 0.2, \dots,} \code{maxk} (the only
#'   grid of fbardl 1.1.0 and earlier). See Details.
#' @param kstar Optional Fourier frequency fixed by the user in advance
#'   (requires \code{fourier = TRUE}). It is then not selected, neither on
#'   the data nor in the bootstrap.
#' @param lags Optional list with elements \code{p} (number of lagged
#'   differences of y) and \code{q} (highest lag of the differences of each
#'   regressor, one value per regressor or a single value for all) fixing the
#'   lags in advance. They are then not selected, neither on the data nor in
#'   the bootstrap.
#' @param seed Optional integer seed for the bootstrap; the random number
#'   generator state of the session is restored afterwards. With
#'   \code{NULL} (default) the session's generator is used.
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
#'   \item{best.p}{Selected (or fixed) lag order for the dependent variable}
#'   \item{best.q}{Selected (or fixed) lag orders for the independent variables}
#'   \item{best.kstar}{Selected (or fixed) Fourier frequency}
#'   \item{kgrid, kstar.fixed, lags.fixed, ssr.by.k}{Frequency grid used,
#'     whether k* and the lags were fixed by the user, and the SSR of each
#'     candidate frequency}
#'   \item{F.overall}{F-statistic for overall cointegration test}
#'   \item{t.dependent}{t-statistic on lagged dependent variable}
#'   \item{F.independent}{F-statistic on lagged independent variables}
#'   \item{cointegration}{Cointegration test results: p-values, critical
#'     values, \code{decision} (text) and \code{decision.code}; for the
#'     bootstrap types also the bootstrap distributions (\code{Fov.boot},
#'     \code{t.boot}, \code{Find.boot}), the numbers of valid and failed
#'     replications, \code{reselect} (what was selected again in each
#'     replication) and \code{dgpcheck}}
#'   \item{diagnostics}{Diagnostic test results}
#'   \item{model.fit}{Model fit statistics (R2, AIC, BIC, etc.)}
#'   \item{residuals}{Model residuals}
#'   \item{fitted.values}{Fitted values}
#'   \item{nobs}{Number of observations}
#'   \item{call}{The matched call}
#' }
#'
#' @details
#' \strong{Model.} The equilibrium-correction model regresses
#' \eqn{\Delta y_t} on \eqn{y_{t-1}}, \eqn{x_{t-1}}, \eqn{p} lagged
#' differences of \eqn{y}, the differences of each regressor at lags 0 (1
#' if \code{unconditional = TRUE}) to \eqn{q_j}, the Fourier terms
#' \eqn{\sin(2\pi k^* t/T)} and \eqn{\cos(2\pi k^* t/T)}, a constant and, in
#' PSS cases 4 and 5, a linear trend. The overall F test (Fov) restricts the
#' lagged levels, together with the intercept in case 2 and the trend in case
#' 4; t is the t statistic on \eqn{y_{t-1}} and Find the F test on the
#' lagged regressors.
#'
#' \strong{Selection.} Step 1: unless \code{kstar} is given, \eqn{k^*}
#' minimises the sum of squared residuals of the model with every lag at
#' \code{maxlag}. With \code{kgrid = "integer"} the candidates are
#' \eqn{1, \dots,} \code{maxk}, the integer frequencies chosen by minimum SSR
#' in Enders and Lee (2012); with \code{kgrid = "fractional"} they are
#' \eqn{0.1, 0.2, \dots,} \code{maxk}. The integer default, and whether the
#' fractional grid is appropriate, are package choices pending verification
#' against Yilanci, Bozoklu and Gorus (2020) and Omay (2015). Step 2: unless
#' \code{lags} is given, \eqn{p} (1 to \code{maxlag}) and each \eqn{q_j} (0
#' to \code{maxlag}) minimise the AIC or BIC with \eqn{k^*} fixed.
#'
#' \strong{Bounds (\code{type = "fardl"}).} Without Fourier terms
#' (\code{fourier = FALSE}) Fov and t are compared with the finite-sample
#' critical values and approximate p-values of Kripfganz and Schneider
#' (2020), for the sample size, the number of regressors and the number of
#' short-run coefficients, and the 5\% bounds decision is reported. With
#' Fourier terms no valid bounds exist: no p-values and no decision are
#' given, and the bounds printed are the bounds for the model without
#' Fourier terms (short-run coefficients counted without the Fourier terms);
#' they are not valid with Fourier terms. Use a bootstrap type for
#' inference. Find has no tabulated distribution.
#'
#' \strong{Bootstrap (\code{"fbardl_bvz"}, \code{"fbardl_mcnown"}).} The
#' bootstrap uses the recursive engine shared with the package ardlverse
#' (file \code{R/ardl_boot_engine.R}, engine version 1.1.0). Once, on the
#' data: the restricted equation for \eqn{\Delta y} (under the null of each
#' statistic) and the equation for \eqn{\Delta x} are estimated with the
#' selected (or fixed) \eqn{k^*}, \eqn{p} and \eqn{q}, and their residuals
#' are saved. \code{"fbardl_bvz"} follows Bertelli, Vacca and Zoia (2022,
#' Section 3): separate nulls for Fov, t and Find, the marginal VECM for
#' \eqn{\Delta x} (\eqn{x_{t-1}} and \eqn{p} lags of \eqn{\Delta y} and
#' \eqn{\Delta x}), and initial values drawn as a random block of the data.
#' \code{"fbardl_mcnown"} uses the null of Fov for all statistics (McNown,
#' Sam and Goh, 2018, Steps 1 to 8); MSG write the equation for \eqn{y}
#' without the contemporaneous differences (their eq. 12), which is the form
#' used with \code{unconditional = TRUE}; with the default conditional ECM,
#' applying their joint null to it is a package choice. The
#' \eqn{\Delta x} equation is unrestricted (it also contains \eqn{y_{t-1}}),
#' and the initial values are the first observations of the data (with
#' \code{"fbardl_bvz"} the levels in the random block start the recursion).
#' In each replication: (1) the residual pairs are resampled with
#' replacement and recentred (after each draw for \code{"fbardl_bvz"}; once,
#' before resampling, for \code{"fbardl_mcnown"}); (2) \eqn{x^*} and
#' \eqn{y^*} are generated recursively from the estimated equations, with
#' the Fourier terms (at the data's \eqn{k^*}) and the trend held fixed; (3)
#' \eqn{k^*} is selected again on the bootstrap sample by Step 1 (unless
#' \code{kstar} is given or \code{fourier = FALSE}) and the lags by Step 2
#' (unless \code{lags} is given); (4) the model is estimated with the
#' selected specification and Fov, t and Find are computed exactly as on the
#' data. Repeating the selection in each replication is a package choice
#' (BVZ step 5(c) re-estimates the unrestricted model but does not discuss
#' re-selection; McNown, Sam and Goh do not discuss it either); it makes the bootstrap distribution account for the data-based selection.
#' In the package's Monte Carlo experiments (n = 100, independent random
#' walks, \code{maxlag = 1}, 200 samples, 199 replications, 5\% level), the
#' default call gave sizes 0.075 (Fov), 0.135 (t), 0.110 (Find) and 0.040
#' for the combined decision; holding \eqn{k^*} and the lags fixed at their
#' data-selected values (fbardl 1.1.0) gave 0.375, 0.420 and 0.270. With
#' \code{unconditional = TRUE} the generating equation for \eqn{\Delta y}
#' also omits the contemporaneous differences of the regressors, so the
#' bootstrap model is the estimated model in both forms. P-values are the shares of bootstrap statistics at
#' least as extreme as the observed one; critical values are order
#' statistics (MSG eqs. 15-16, BVZ eqs. 24-25). The decision (5\% level)
#' requires Fov, t and Find to reject. Fov rejecting and t not is the
#' degenerate case of the first type (whatever Find gives); Fov and t
#' rejecting and Find not is the degenerate case of the second type
#' (terminology of Bertelli, Vacca and Zoia 2022, eqs. 8-9; McNown, Sam and
#' Goh 2018 number the two degenerate cases the other way). The same labels
#' are used for both bootstrap types. Failed replications are set to \code{NA} and counted.
#' The run time grows with \code{reps}, \code{maxlag}, the number of
#' regressors and the frequency grid, because the selection is repeated in
#' every replication; fixing \code{kstar} and \code{lags} removes it.
#'
#' @references
#' Bertelli, S., Vacca, G. and Zoia, M. (2022). Bootstrap cointegration
#' tests in ARDL models. \emph{Economic Modelling}, 116, 105987.
#' \doi{10.1016/j.econmod.2022.105987}
#'
#' Enders, W. and Lee, J. (2012). A unit root test using a Fourier series to
#' approximate smooth breaks. \emph{Oxford Bulletin of Economics and
#' Statistics}, 74(4), 574-599. \doi{10.1111/j.1468-0084.2011.00662.x}
#'
#' Kripfganz, S. and Schneider, D. C. (2020). Response surface regressions
#' for critical value bounds and approximate p-values in equilibrium
#' correction models. \emph{Oxford Bulletin of Economics and Statistics},
#' 82(6), 1456-1481. \doi{10.1111/obes.12377}
#'
#' McNown, R., Sam, C. Y. and Goh, S. K. (2018). Bootstrapping the
#' autoregressive distributed lag test for cointegration. \emph{Applied
#' Economics}, 50(13), 1509-1521. \doi{10.1080/00036846.2017.1366643}
#'
#' Omay, T. (2015). Fractional frequency flexible Fourier form to
#' approximate smooth breaks in unit root testing. \emph{Economics Letters},
#' 134, 123-126. \doi{10.1016/j.econlet.2015.07.010}
#'
#' Pesaran, M. H., Shin, Y. and Smith, R. J. (2001). Bounds testing
#' approaches to the analysis of level relationships. \emph{Journal of
#' Applied Econometrics}, 16(3), 289-326. \doi{10.1002/jae.616}
#'
#' Yilanci, V., Bozoklu, S. and Gorus, M. S. (2020). Are BRICS countries
#' pollution havens? Evidence from a bootstrap ARDL bounds testing approach
#' with a Fourier function. \emph{Sustainable Cities and Society}, 55, 102035.
#' \doi{10.1016/j.scs.2020.102035}
#'
#' @examples
#' \donttest{
#' data(fbardl_data)
#'
#' # Bootstrap (Bertelli, Vacca and Zoia), k* and lags selected again in
#' # each replication; few replications to keep the example short
#' result <- fbardl(y ~ x1 + x2, data = fbardl_data, maxlag = 2, maxk = 3,
#'                  reps = 49, seed = 1)
#' result
#'
#' # McNown, Sam and Goh bootstrap with k* and lags fixed in advance
#' result_msg <- fbardl(y ~ x1 + x2, data = fbardl_data,
#'                      type = "fbardl_mcnown", kstar = 1,
#'                      lags = list(p = 1, q = c(1, 1)), reps = 99, seed = 1)
#' summary(result_msg)
#'
#' # Bounds test without Fourier terms
#' result_bounds <- fbardl(y ~ x1 + x2, data = fbardl_data, type = "fardl",
#'                         fourier = FALSE)
#' }
#'
#' @export
fbardl <- function(formula, data, type = c("fbardl_bvz", "fbardl_mcnown", "fardl"),
                   maxlag = 4, maxk = 5, ic = c("aic", "bic"), case = 3,
                   reps = 999, fourier = TRUE, level = 0.95, horizon = 20,
                   unconditional = FALSE, kgrid = c("integer", "fractional"),
                   kstar = NULL, lags = NULL, seed = NULL) {

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
  if (type != "fardl" && (!is.numeric(reps) || length(reps) != 1 || reps < 1 ||
                          reps != round(reps)))
    stop("'reps' must be a positive integer")
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

  kgrid <- match.arg(kgrid)
  if (!is.null(kstar)) {
    if (!isTRUE(fourier) || !is.numeric(kstar) || length(kstar) != 1 ||
        !is.finite(kstar) || kstar <= 0)
      stop("'kstar' must be a single positive number and requires fourier = TRUE")
  }
  if (isTRUE(fourier) && is.null(kstar)) {
    if (kgrid == "integer" && (maxk < 1 || maxk != round(maxk)))
      stop("'maxk' must be a positive integer with kgrid = \"integer\"")
    if (kgrid == "fractional" && maxk < 0.1)
      stop("'maxk' must be at least 0.1 with kgrid = \"fractional\"")
  }
  if (!is.null(lags)) {
    if (!is.list(lags) || !all(c("p", "q") %in% names(lags)))
      stop("'lags' must be a list with elements 'p' and 'q'")
    lp <- lags$p
    lq <- lags$q
    if (length(lq) == 1) lq <- rep(lq, nindep)
    if (length(lp) != 1 || !is.numeric(lp) || lp < 0 || lp != round(lp))
      stop("'lags$p' must be a non-negative integer")
    if (length(lq) != nindep || !is.numeric(lq) || any(lq < 0) || any(lq != round(lq)))
      stop("'lags$q' must hold one non-negative integer per regressor")
    lags <- list(p = as.integer(lp), q = as.integer(lq))
  }

  setup <- list(
    maxlag = maxlag, hastrend = hastrend, j0 = j0, case = case, ic = ic,
    fourier = isTRUE(fourier), kstar = kstar, lags = lags,
    kvalues = if (kgrid == "integer") seq_len(max(0, floor(maxk))) else
      seq(0.1, max(0.1, maxk), by = 0.1))

  T <- n_original  # For Fourier terms

  # ============================================================================
  # STEPS 1 AND 2: SELECT k* BY MINIMUM SSR, THEN (p, q) BY AIC/BIC
  # ============================================================================
  if (!fourier) {
    message("Step 1: Fourier terms excluded (k* = 0)")
  } else if (!is.null(kstar)) {
    message(sprintf("Step 1: Fourier frequency fixed by the user: k* = %g", kstar))
  } else {
    message(sprintf("Step 1: Selecting Fourier frequency k* by minimum SSR (%s grid)...",
                    kgrid))
  }
  sel <- .fbardl_select(y, X, setup)
  best_kstar <- sel$kstar
  best_ssr_k <- sel$best_ssr_k
  if (fourier && is.null(kstar))
    message(sprintf("  Optimal k* = %.1f (min SSR = %.4f)", best_kstar, best_ssr_k))

  if (is.null(lags)) {
    message(sprintf("Step 2: Selecting lag orders (p, q) by %s...", toupper(ic)))
  } else {
    message("Step 2: Lag orders fixed by the user")
  }
  best_p <- sel$p
  best_q <- sel$q
  names(best_q) <- indepvars
  best_ic_val <- sel$ic
  total_specs <- sel$total_specs
  message(sprintf("  p = %d, q = (%s) (%d models evaluated)",
                  best_p, paste(best_q, collapse = ", "), total_specs))

  fs <- .fbardl_fourier(T, best_kstar)
  fourier_sin <- fs$sin
  fourier_cos <- fs$cos

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

  ind_names <- coef_names[indep_level_idx]
  if (type == "fardl") {
    # Short-run coefficients: every regressor other than the constant, the
    # trend, the k + 1 lagged levels and the Fourier terms. With Fourier terms
    # no valid bounds exist: the bounds of the model without them are shown
    # for reference only, with no p-values and no decision.
    has_fourier <- best_kstar > 0
    sr <- nparams - 1L - (nindep + 1L) - as.integer(hastrend) - 2L * has_fourier
    coint_result <- .pss_bounds_test(Fov_stat, t_stat, nindep, case, nobs, sr,
                                     valid = !has_fourier)
  } else {
    message(sprintf("Step 4: Bootstrap (%d replications)...", reps))
    coint_result <- .fbardl_bootstrap(y, X, sel, setup, type, reps,
                                      fov_names, ind_names, seed = seed)
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
    kgrid = if (fourier && is.null(kstar)) kgrid else NA_character_,
    kstar.fixed = !is.null(kstar),
    lags.fixed = !is.null(lags),
    ssr.by.k = sel$ssr_by_k,
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
      paste(x$best.q, collapse = ","), ")",
      if (isTRUE(x$lags.fixed)) " (lags fixed by the user)" else "", "\n", sep = "")
  cat("Fourier frequency (k*): ", x$best.kstar,
      if (isTRUE(x$kstar.fixed)) " (fixed by the user)" else
        if (isTRUE(x$fourier)) paste0(" (", x$kgrid, " grid)") else " (no Fourier terms)",
      "\n", sep = "")
  cat("Test type: ", x$type, "\n", sep = "")
  cat("Sample size: ", x$nobs, "\n", sep = "")

  cat("\nCointegration Test Results:\n")
  cat(rep("-", 50), "\n", sep = "")
  co <- x$cointegration
  pv <- function(v) if (length(v) == 2) {
    if (anyNA(v)) "p-value: n/a" else sprintf("p I(0)/I(1): %6.4f / %6.4f", v[1], v[2])
  } else if (is.na(v)) "p-value: n/a" else sprintf("p-value: %6.4f", v)
  cat(sprintf("F-overall:     %8.4f  (%s)\n", x$F.overall, pv(co$Fov.pval)))
  cat(sprintf("t-dependent:   %8.4f  (%s)\n", x$t.dependent, pv(co$t.pval)))
  cat(sprintf("F-independent: %8.4f  (%s)\n", x$F.independent, pv(co$Find.pval)))
  cat("Critical values: ", co$source, "\n", sep = "")
  if (identical(co$source, "bootstrap")) {
    cat("Scheme: ", co$scheme, "; ", x$reps, " replications", "\n", sep = "")
    redo <- c(if (co$reselect["kstar"]) "k*", if (co$reselect["lags"]) "lags")
    cat("Selected again in each replication: ",
        if (length(redo)) paste(redo, collapse = " and ") else "nothing (k* and lags fixed)",
        "\n", sep = "")
  } else if (isFALSE(co$valid)) {
    cat("Bounds for the model without Fourier terms; not valid with Fourier terms.\n")
    cat("Use type = \"fbardl_bvz\" or \"fbardl_mcnown\" for inference.\n")
  }
  cat("\nDecision: ", x$cointegration$decision, "\n", sep = "")
  if (identical(co$source, "bootstrap") && grepl("^DEGENERATE", co$decision.code))
    cat("(Degenerate cases numbered as in Bertelli, Vacca and Zoia, 2022;",
        "McNown, Sam and Goh, 2018, number them the other way.)\n")

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
    cat("  Fourier Bootstrap ARDL - McNown, Sam and Goh (2018)\n")
  } else {
    cat("  Fourier Bootstrap ARDL - Bertelli, Vacca and Zoia (2022)\n")
  }

  cat(rep("=", 78), "\n", sep = "")
  cat("  Dependent variable  : ", object$depvar, "\n", sep = "")
  cat("  Independent var(s)  : ", paste(object$indepvars, collapse = ", "), "\n", sep = "")
  cat("  Sample size (T)     : ", object$nobs, "\n", sep = "")
  cat("  Lags of D.y (p)     : ", object$best.p, "\n", sep = "")
  if (object$fourier) {
    cat("  Fourier frequency   : ", object$best.kstar, "\n", sep = "")
  } else {
    cat("  Fourier terms       : excluded\n")
  }
  cat("  Information crit.   : ", toupper(object$ic), "\n", sep = "")
  cat("  PSS Case            : Case ", object$case, "\n", sep = "")
  if (object$type != "fardl") {
    cat("  Bootstrap reps      : ", object$reps, "\n", sep = "")
    redo <- c(if (object$cointegration$reselect["kstar"]) "k*",
              if (object$cointegration$reselect["lags"]) "lags")
    cat("  Re-selected per rep : ",
        if (length(redo)) paste(redo, collapse = ", ") else "none (fixed)", "\n", sep = "")
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
    if (isFALSE(coint$valid)) {
      cat("  These are the bounds for the model without Fourier terms; not valid\n")
      cat("  with Fourier terms (no p-values, no decision). Use type = \"fbardl_bvz\"\n")
      cat("  or \"fbardl_mcnown\" for inference.\n")
    }
  }
  cat("\n  Decision: ", coint$decision, "\n", sep = "")
  if (identical(coint$source, "bootstrap") && grepl("^DEGENERATE", coint$decision.code))
    cat("  (Degenerate cases numbered as in Bertelli, Vacca and Zoia, 2022;\n",
        "  McNown, Sam and Goh, 2018, number them the other way.)\n", sep = "")
  cat(rep("-", 68), "\n", sep = "")

  # Table 4: Diagnostic Tests
  cat("\nTable 4: Diagnostic Tests\n")
  cat(rep("-", 68), "\n", sep = "")
  diag <- object$diagnostics
  cat(sprintf("  Jarque-Bera (Normality)      : %.4f  (p = %.4f)%s\n",
              diag$jb.stat, diag$jb.pval, .stars(diag$jb.pval)))
  cat(sprintf("  Breusch-Godfrey (AR 1)       : %.4f  (p = %.4f)%s\n",
              diag$bg1.stat, diag$bg1.pval, .stars(diag$bg1.pval)))
  cat(sprintf("  Breusch-Godfrey (AR 4)       : %.4f  (p = %.4f)%s\n",
              diag$bg4.stat, diag$bg4.pval, .stars(diag$bg4.pval)))
  cat(sprintf("  Breusch-Pagan (Heterosked.)  : %.4f  (p = %.4f)%s\n",
              diag$bp.stat, diag$bp.pval, .stars(diag$bp.pval)))
  cat(sprintf("  ARCH(1) Test                 : %.4f  (p = %.4f)%s\n",
              diag$arch1.stat, diag$arch1.pval, .stars(diag$arch1.pval)))
  cat(rep("-", 68), "\n", sep = "")

  cat("\n")
  cat(rep("=", 78), "\n", sep = "")
  cat("  fbardl ", getNamespaceVersion("fbardl"), "\n", sep = "")
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
