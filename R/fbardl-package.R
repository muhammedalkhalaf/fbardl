#' @keywords internal
#' @importFrom stats coef embed fitted lm lm.fit model.frame model.matrix
#'   model.response na.pass pchisq pf pnorm pt qf qnorm quantile vcov
"_PACKAGE"

#' fbardl: Fourier Bootstrap ARDL Cointegration Test
#'
#' @description
#' The fbardl package implements ARDL bounds tests for cointegration
#' (Pesaran, Shin and Smith, 2001) with optional Fourier terms for smooth
#' breaks and recursive bootstrap critical values (McNown, Sam and Goh, 2018;
#' Bertelli, Vacca and Zoia, 2022). With Fourier terms only the bootstrap
#' types give valid inference.
#'
#' @section Main Function:
#' \itemize{
#'   \item \code{\link{fbardl}}: Perform Fourier Bootstrap ARDL cointegration test
#' }
#'
#' @section Test Types:
#' \itemize{
#'   \item \code{"fbardl_bvz"} (default): Bootstrap ARDL (Bertelli, Vacca and Zoia, 2022)
#'   \item \code{"fbardl_mcnown"}: Bootstrap ARDL (McNown, Sam and Goh, 2018)
#'   \item \code{"fardl"}: Kripfganz and Schneider (2020) bounds, valid
#'     without Fourier terms only
#' }
#'
#' @section Data:
#' \itemize{
#'   \item \code{\link{fbardl_data}}: Example dataset for demonstration
#' }
#'
#' @references
#' Bertelli, S., Vacca, G. and Zoia, M. (2022). Bootstrap cointegration
#' tests in ARDL models. \emph{Economic Modelling}, 116, 105987.
#' \doi{10.1016/j.econmod.2022.105987}
#'
#' McNown, R., Sam, C. Y. and Goh, S. K. (2018). Bootstrapping the
#' autoregressive distributed lag test for cointegration. \emph{Applied
#' Economics}, 50(13), 1509-1521. \doi{10.1080/00036846.2017.1366643}
#'
#' Pesaran, M. H., Shin, Y. and Smith, R. J. (2001). Bounds testing
#' approaches to the analysis of level relationships. \emph{Journal of
#' Applied Econometrics}, 16(3), 289-326. \doi{10.1002/jae.616}
#'
#' @docType package
#' @name fbardl-package
NULL
