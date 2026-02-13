#' @export
BAR_threshold <- function(resid, lambda, n) {
  if (resid <= 2 * sqrt(lambda / n) && resid >= -2 * sqrt(lambda / n)) {
    0
  } else {
    resid / 2 + sqrt(resid^2 / 4 - lambda / n)
  }
}
#' @export
Soft <- function(z, lambda) {
  dplyr::case_when(
    z >  lambda ~ z - lambda,
    z < -lambda ~ z + lambda,
    TRUE        ~ 0
  )
}

#' @export
aftPenCDA <- function(dt, lambda, se, type, r = 3.7, eps = 1e-8, max.iter = 100) {
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
    tilde_beta <- init
    resid <- u - z %*% beta

    for (t in seq_len(max.iter)) {
      new.beta <- beta

      for (j in seq_len(p)) {
        zj <- (crossprod(z[, j], resid) / n) + beta[j]
        new.beta[j] <- BAR_threshold(zj, lambda, n)
        resid <- resid - z[, j] * (new.beta[j] - beta[j])
      }

      if (max(abs(beta - new.beta)) < eps) break
      tilde_beta <- beta
      beta <- new.beta
    }

    return(list(
      beta  = beta / sqrt(norm.z),
      tbeta = tilde_beta / sqrt(norm.z)
    ))
  }

  if (type == "ALASSO") tilde.beta <- init
  if (type == "LASSO")  tilde.beta <- rep(1, length(init))

  init2 <- coef(lm(u ~ -1 + z))
  init2 <- ifelse(is.na(init2), eps, init2)

  beta <- init2
  tilde_beta <- init2
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
        Soft(zj, lambda / abs(tilde.beta[j]))
      }

      resid <- resid - z[, j] * (new.beta[j] - beta[j])
    }

    if (max(abs(beta - new.beta)) < eps) break
    tilde_beta <- beta
    beta <- new.beta
  }

  list(
    beta  = beta / sqrt(norm.z),
    tbeta = tilde_beta / sqrt(norm.z)
  )
}
