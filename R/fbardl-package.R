#' @keywords internal
#' @importFrom stats coef embed fitted lm lm.fit model.frame model.matrix
#'   model.response na.pass pchisq pf pnorm pt qf qnorm quantile vcov
"_PACKAGE"

#' fbardl: Fourier Bootstrap ARDL Cointegration Test
#'
#' @description
#' The fbardl package implements the Fourier Bootstrap ARDL bounds testing
#' approach for cointegration analysis. It combines the Pesaran, Shin & Smith
#' (2001) ARDL framework with Fourier terms to capture structural breaks,
#' and provides bootstrap critical values for robust inference.
#'
#' @section Main Function:
#' \itemize{
#'   \item \code{\link{fbardl}}: Perform Fourier Bootstrap ARDL cointegration test
#' }
#'
#' @section Test Types:
#' \itemize{
#'   \item \code{"fardl"}: Standard Fourier ARDL with PSS bounds test
#'   \item \code{"fbardl_mcnown"}: Bootstrap ARDL (McNown, Sam & Goh, 2018)
#'   \item \code{"fbardl_bvz"}: Bootstrap ARDL (Bertelli, Vacca & Zoia, 2022)
#' }
#'
#' @section Data:
#' \itemize{
#'   \item \code{\link{fbardl_data}}: Example dataset for demonstration
#' }
#'
#' @references
#' Pesaran, M. H., Shin, Y., & Smith, R. J. (2001). Bounds testing approaches
#' to the analysis of level relationships. \emph{Journal of Applied Econometrics},
#' 16(3), 289-326. \doi{10.1002/jae.616}
#'
#' McNown, R., Sam, C. Y., & Goh, S. K. (2018). Bootstrapping the autoregressive
#' distributed lag test for cointegration. \emph{Applied Economics}, 50(13),
#' 1509-1521. \doi{10.1080/00036846.2017.1366643}
#'
#' @docType package
#' @name fbardl-package
NULL
