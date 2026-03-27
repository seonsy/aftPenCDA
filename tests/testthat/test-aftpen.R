test_that("aftpen runs and returns beta", {
  set.seed(1)

  n <- 100
  p <- 10
  beta0 <- rep(1, p)
  x <- matrix(rnorm(n * p), n, p)
  T <- exp(x %*% beta0 + rnorm(n))
  C <- rexp(n, rate = exp(-2))
  d <- 1 * (T < C)
  y <- pmin(T, C)

  dt <- data.frame(y = y, d = d, x)

  fit <- aftpen(dt, lambda = 0.1, se = "CF", type = "BAR")

  expect_true(is.list(fit))
  expect_true("beta" %in% names(fit))
  expect_equal(length(fit$beta), 10)
  expect_true(all(is.finite(fit$beta)))
})
