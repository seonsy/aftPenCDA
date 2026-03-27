#' @importFrom stats coef lm
NULL

#' Internal BAR thresholding operator
#'
#' Computes the thresholding update used in the BAR penalty.
#'
#' @param resid A scalar update quantity.
#' @param lambda A tuning parameter.
#' @param n Sample size used in the thresholding rule.
#'
#' @return A scalar thresholded value.
#' @keywords internal
#'
#'
BAR_threshold <- function(resid, lambda, n) {
  if (resid <= 2 * sqrt(lambda / n) && resid >= -2 * sqrt(lambda / n)) {
    0
  } else {
    resid / 2 + sqrt(resid^2 / 4 - lambda / n)
  }
}

#' Internal soft-thresholding operator
#'
#' Computes the soft-thresholding update used in SCAD, ALASSO, and LASSO penalties.
#'
#' @param z A numeric scalar or vector.
#' @param lambda A nonnegative thresholding parameter.
#'
#' @return A numeric value (or vector) after soft-thresholding.
#' @keywords internal

Soft <- function(z, lambda) {
  ifelse(z > lambda, z - lambda,
         ifelse(z < -lambda, z + lambda, 0))
}

#' Penalized AFT estimation for right-censored data
#'
#' Fits a penalized accelerated failure time (AFT) model for right-censored
#' survival data using induced smoothing and a penalized coordinate descent
#' algorithm. Supported penalties include BAR, LASSO, adaptive LASSO, and SCAD.
#'
#' @param dt A data frame whose first two columns are \code{y} and \code{d},
#'   where \code{y} is the observed survival or censoring time and \code{d}
#'   is the event indicator (\code{1} for observed event and \code{0} for
#'   right-censoring). The remaining columns are treated as covariates.
#' @param lambda A nonnegative tuning parameter controlling the amount of
#'   penalization.
#' @param se A character string specifying the variance estimation method.
#'   \describe{
#'     \item{"CF"}{Closed-form (plug-in) variance estimator based on the analytic expression of the estimating function (see Equation (6)).}
#'     \item{"ZL"}{Perturbation-based variance estimator using the resampling approach of Zeng and Lin (2008).}
#'   }
#' @param type Penalty type. One of \code{"BAR"}, \code{"LASSO"},
#'   \code{"ALASSO"}, or \code{"SCAD"}.
#' @param r A positive tuning constant used in the SCAD penalty. Ignored unless
#'   \code{type = "SCAD"}. The default is \code{3.7}.
#' @param eps Convergence tolerance for the outer penalized coordinate descent
#'   iterations. The default is \code{1e-8}.
#' @param max.iter Maximum number of iterations for the outer penalized
#'   coordinate descent algorithm. The default is \code{100}.
#'
#' @return A list containing the following components:
#' \itemize{
#'   \item \code{beta}: final coefficient estimate on the original scale.
#' }
#'
#' @details
#' The function first calls the Rcpp backend \code{is_aft_cpp()} to obtain
#' an initial estimator together with gradient and Hessian information.
#' A Cholesky-based transformation is then applied, followed by coordinate-wise
#' penalized updates.
#'
#' For \code{type = "BAR"}, the update uses the internal \code{BAR_threshold()} operator.
#' For \code{type = "LASSO"}, \code{"ALASSO"}, and \code{"SCAD"}, soft-thresholding-based updates are used.
#' @examples
#' \dontrun{
#' n = 100
#' p = 5
#' beta0 = rep(1,p)
#' x = matrix(rnorm(n * p), n, p)
#' T = exp(x%*%beta0 + rnorm(n))
#' C = rexp(n, rate = exp(-2))
#' d = 1*(T<C)
#' y = pmin(T,C)
#' dt = data.frame(y,d,x)
#' fit <- aftpen(dt, lambda = 0.1, se = "CF", type = "BAR")
#' fit$beta
#' }
#'
#' @export
aftpen <- function(dt, lambda, se, type = c("BAR", "LASSO", "ALASSO", "SCAD"), r = 3.7, eps = 1e-8, max.iter = 100) {
  type <- match.arg(type)
  fit  <- is_aft_cpp(dt$y, dt$d, as.matrix(dt[, -c(1, 2)]), se)
  init <- fit$beta
  h    <- fit$hess
  g    <- fit$grad

  x_chol <- chol(h)
  y_tilde <- forwardsolve(t(x_chol), h %*% init - g)

  n <- length(y_tilde)
  p <- ncol(x_chol)

  u <- y_tilde - mean(y_tilde)

  z <- t(t(x_chol) - apply(x_chol, 2, mean))
  norm.z <- apply(z^2, 2, mean)
  z <- t(t(z) / sqrt(norm.z))

  if (type == "BAR") {
    beta <- init
    beta_prev <- init
    resid <- u - z %*% beta

    for (t in seq_len(max.iter)) {
      new.beta <- beta

      for (j in seq_len(p)) {
        zj <- (crossprod(z[, j], resid) / n) + beta[j]
        new.beta[j] <- BAR_threshold(zj, lambda, n)
        resid <- resid - z[, j] * (new.beta[j] - beta[j])
      }

      if (max(abs(beta - new.beta)) < eps) break
      beta_prev <- beta
      beta <- new.beta
    }

    return(list(
      beta  = beta / sqrt(norm.z)
    ))
  }

  if (type == "ALASSO") beta_weight <- init
  if (type == "LASSO")  beta_weight <- rep(1, length(init))

  init2 <- coef(lm(u ~ -1 + z))
  init2 <- ifelse(is.na(init2), eps, init2)

  beta <- init2
  beta_prev <- init2
  resid <- u - z %*% beta

  for (t in seq_len(max.iter)) {
    new.beta <- beta

    for (j in seq_len(p)) {
      zj <- crossprod(z[, j], resid) / n + beta[j]

      new.beta[j] <- if (type == "SCAD") {
        Soft(zj, lambda) * (abs(zj) <= 2 * lambda) +
          Soft(zj, r * lambda / (r - 1)) / (1 - 1 / (r - 1)) *
          (2 * lambda < abs(zj) & abs(zj) <= r * lambda) +
          zj * (abs(zj) > r * lambda)
      } else {
        Soft(zj, lambda / abs(beta_weight[j]))
      }

      resid <- resid - z[, j] * (new.beta[j] - beta[j])
    }

    if (max(abs(beta - new.beta)) < eps) break
    beta_prev <- beta
    beta <- new.beta
  }

  list(
    beta  = beta / sqrt(norm.z)
  )
}

#' Penalized AFT estimation for clustered-partly interval censored data
#'
#' Fits a penalized accelerated failure time (AFT) model for partly interval censored
#' survival data using induced smoothing and a penalized coordinate descent
#' algorithm. Supported penalties include BAR, LASSO, adaptive LASSO, and SCAD.
#'
#' @param dt A data frame containing PIC survival data. It must include
#'   \code{L}, \code{R}, \code{delta}, and \code{id}, where \code{L} and
#'   \code{R} define the observation interval, \code{delta} indicates
#'   whether the failure time is exactly observed (\code{1}) or censored
#'   (\code{0}), and \code{id} is the cluster identifier. The remaining
#'   columns are treated as covariates.
#' @param lambda A nonnegative tuning parameter controlling the amount of
#'   penalization.
#' @param se A character string specifying the variance estimation method.
#'   \describe{
#'     \item{"CF"}{Closed-form (analytic plug-in) variance estimator based on the estimating function.}
#'     \item{"ZL"}{Perturbation-resampling variance estimator following Zeng and Lin (2008).}
#'   }
#' @param type Penalty type. One of \code{"BAR"}, \code{"LASSO"}, \code{"ALASSO"}, or \code{"SCAD"}.
#' @param r A positive tuning constant used in the SCAD penalty. Ignored unless
#'   \code{type = "SCAD"}. The default is \code{3.7}.
#' @param eps Convergence tolerance for the outer penalized coordinate descent
#'   iterations. The default is \code{1e-8}.
#' @param max.iter Maximum number of iterations for the outer penalized
#'   coordinate descent algorithm. The default is \code{100}.
#'
#' @return A list containing the following components:
#' \itemize{
#'   \item \code{beta}: final coefficient estimate on the original scale.
#' }
#'
#' @details
#' The input data \code{dt} are assumed to arise from clustered partly
#' interval-censored survival data with informative cluster sizes.
#'
#' Specifically, observations are grouped into clusters, where each cluster
#' shares a latent frailty variable that affects both the failure times and the cluster size.
#' As a result, the number of observations within each
#' cluster is not fixed but depends on the underlying frailty, leading to an informative cluster size structure.
#'
#' For each subject, the failure time follows an accelerated failure time (AFT)
#' model, and the observed data consist of an interval \eqn{(L, R)} together
#' with an indicator \code{delta}. When \eqn{L = R} (i.e., \code{delta = 1}),
#' the observation is exact; otherwise (\code{delta = 0}), the observation is
#' censored and may correspond to left-censoring, right-censoring, or
#' interval-censoring depending on the relationship between the true failure
#' time and the inspection times.
#'
#' The function first calls the Rcpp backend \code{is_aftp_pic_cpp()} to obtain
#' an initial estimator together with gradient and Hessian information.
#' A Cholesky-based transformation is then applied, followed by coordinate-wise
#' penalized updates.
#'
#' For \code{type = "BAR"}, the update uses the internal
#' \code{BAR_threshold()} operator. For \code{"LASSO"},
#' \code{"ALASSO"}, and \code{"SCAD"}, soft-thresholding-based updates are used.
#' @examples
#' \dontrun{
#' set.seed(1)
#'
#' ## simplified generator for clustered partly interval-censored data
#' n <- 100
#' p <- 2
#' beta0 <- c(1,1)
#' clu_rate <- 0.5
#' exactrates <- 0.8
#' left <- 0.001
#' right <- 0.01
#'
#' ## cluster-level frailty and informative cluster sizes
#' eta <- 1 / clu_rate
#' v <- rgamma(n, shape = eta, rate = eta)
#' m <- ifelse(v > median(v), 5, 3)
#' id <- rep(seq_len(n), m)
#' vi <- rep(v, m)
#'
#' ## subject-level covariates and failure times
#' N <- sum(m)
#' x <- matrix(rnorm(N * p), ncol = p)
#' T <- as.vector(exp(x %*% beta0 + vi * log(rexp(N))))
#'
#' ## build (L, R, delta)
#' L <- R <- delta <- numeric(N)
#' index <- rbinom(N, 1, exactrates)
#'
#' for (i in seq_len(N)) {
#'   if (index[i] == 1) {
#'     L[i] <- T[i]
#'     R[i] <- T[i]
#'     delta[i] <- 1
#'   } else {
#'     U <- cumsum(c(1e-8, runif(10, left, right)))
#'     LL <- U[-length(U)]
#'     RR <- U[-1]
#'
#'     if (T[i] < min(LL)) {
#'       L[i] <- 1e-8
#'       R[i] <- min(LL)
#'       delta[i] <- 0
#'     } else if (T[i] > max(RR)) {
#'       L[i] <- max(RR)
#'       R[i] <- 1e8
#'       delta[i] <- 0
#'     } else {
#'       idd <- which(T[i] > LL & T[i] < RR)
#'       if (length(idd) == 1) {
#'         L[i] <- LL[idd]
#'         R[i] <- RR[idd]
#'         delta[i] <- 0
#'       } else {
#'         L[i] <- T[i]
#'         R[i] <- T[i]
#'         delta[i] <- 1
#'       }
#'     }
#'   }
#' }
#'
#' dt <- data.frame(
#'   L = L, R = R, delta = delta, id = id,
#'   x1 = x[, 1], x2 = x[, 2]
#' )
#'
#' fit <- aftpen_pic(dt, lambda = 0.1, se = "CF", type = "BAR")
#' fit$beta
#' }
#' @export
aftpen_pic <- function(dt, lambda, se,
                         type = c("BAR", "LASSO", "ALASSO", "SCAD"), r = 3.7, eps = 1e-8, max.iter = 100) {
  type <- match.arg(type)
  fit  <- is_aftp_pic_cpp(dt, se)
  init <- fit$beta
  h    <- fit$hess
  g    <- fit$grad

  x_chol <- chol(h)
  y_tilde <- forwardsolve(t(x_chol), h %*% init - g)

  n <- length(y_tilde)
  p <- ncol(x_chol)

  u <- y_tilde - mean(y_tilde)

  z <- t(t(x_chol) - apply(x_chol, 2, mean))
  norm.z <- apply(z^2, 2, mean)
  z <- t(t(z) / sqrt(norm.z))

  if (type == "BAR") {
    beta <- init
    beta_prev <- init
    resid <- u - z %*% beta

    for (t in seq_len(max.iter)) {
      new.beta <- beta

      for (j in seq_len(p)) {
        zj <- (crossprod(z[, j], resid) / n) + beta[j]
        new.beta[j] <- BAR_threshold(zj, lambda, n)
        resid <- resid - z[, j] * (new.beta[j] - beta[j])
      }

      if (max(abs(beta - new.beta)) < eps) break
      beta_prev <- beta
      beta <- new.beta
    }

    return(list(
      beta  = beta / sqrt(norm.z)
    ))
  }

  if (type == "ALASSO") beta_weight = init
  if (type == "LASSO")  beta_weight = rep(1, length(init))

  init2 <- coef(lm(u ~ -1 + z))
  init2 <- ifelse(is.na(init2), eps, init2)

  beta <- init2
  beta_prev <- init2
  resid <- u - z %*% beta

  for (t in seq_len(max.iter)) {
    new.beta <- beta

    for (j in seq_len(p)) {
      zj <- crossprod(z[, j], resid) / n + beta[j]

      new.beta[j] <- if (type == "SCAD") {
        Soft(zj, lambda) * (abs(zj) <= 2 * lambda) +
          Soft(zj, r * lambda / (r - 1)) / (1 - 1 / (r - 1)) *
          (2 * lambda < abs(zj) & abs(zj) <= r * lambda) +
          zj * (abs(zj) > r * lambda)
      } else {
        Soft(zj, lambda / abs(beta_weight[j]))
      }

      resid <- resid - z[, j] * (new.beta[j] - beta[j])
    }

    if (max(abs(beta - new.beta)) < eps) break
    beta_prev <- beta
    beta <- new.beta
  }

  list(
    beta  = beta / sqrt(norm.z))
}
