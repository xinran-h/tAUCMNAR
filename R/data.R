#' Demo data
#'
#' @description
#' A data used for demo. 
#'
#' @format ## `demo`
#' A dataframe of 1000 rows and 7 columns:
#' \describe{
#'   \item{ID}{Patient ID}
#'   \item{Y}{Observed event time}
#'   \item{delta}{Censoring indicator. 0 = censored, 1 = event 1}
#'   \item{z}{A numeric vector indicating a continuous non-instrument covariate}
#'   \item{u}{A numeric vector indicating a discrete instrument covariate}
#'   \item{M}{Baseline biomarker values}
#'   \item{V}{A missing indicator for biomarker values, V = 1 if observed, V = 0 if missing}
#' }
"demo"


