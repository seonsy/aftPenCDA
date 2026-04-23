#' Simulated right-censored data
#'
#' Example dataset for penalized AFT model fitting with right-censored survival data.
#'
#' @format A data frame with 100 rows and 12 variables:
#' \describe{
#'   \item{y}{Observed survival or censoring time}
#'   \item{d}{Censoring indicator (1 = event observed, 0 = censored)}
#'   \item{X1}{First covariate}
#'   \item{X2}{Second covariate}
#'   \item{X3}{Third covariate}
#'   \item{X4}{Fourth covariate}
#'   \item{X5}{Fifth covariate}
#'   \item{X6}{Sixth covariate}
#'   \item{X7}{Seventh covariate}
#'   \item{X8}{Eighth covariate}
#'   \item{X9}{Ninth covariate}
#'   \item{X10}{Tenth covariate}
#' }
#' @source Simulated data
"simdat_rc"


#' Simulated clustered partly interval-censored data
#'
#' Example dataset for penalized AFT model fitting with clustered partly interval-censored survival data.
#'
#' @format A data frame with 6 variables:
#' \describe{
#'   \item{L}{Left endpoint of observation interval}
#'   \item{R}{Right endpoint of observation interval}
#'   \item{delta}{Exact observation indicator (1 = exact, 0 = interval)}
#'   \item{id}{Cluster identifier}
#'   \item{x1}{First covariate}
#'   \item{x2}{Second covariate}
#' }
#' @source Simulated data
"simdat_pic"
