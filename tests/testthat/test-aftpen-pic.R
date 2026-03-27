test_that("aftpen_pic runs and returns valid output", {
  set.seed(1)

  make_pic_data <- function() {
    n <- 100
    p <- 2
    beta0 <- c(1, 1)
    clu_rate <- 0.5
    exactrates <- 0.8
    left <- 0.001
    right <- 0.01

    eta <- 1 / clu_rate
    v <- rgamma(n, shape = eta, rate = eta)
    m <- ifelse(v > median(v), 5, 3)
    id <- rep(seq_len(n), m)
    vi <- rep(v, m)

    N <- sum(m)
    x <- matrix(rnorm(N * p), ncol = p)
    T <- as.vector(exp(x %*% beta0 + vi * log(rexp(N))))

    L <- R <- delta <- numeric(N)
    index <- rbinom(N, 1, exactrates)

    for (i in seq_len(N)) {
      if (index[i] == 1) {
        L[i] <- T[i]
        R[i] <- T[i]
        delta[i] <- 1
      } else {
        U <- cumsum(c(1e-8, runif(10, left, right)))
        LL <- U[-length(U)]
        RR <- U[-1]

        if (T[i] < min(LL)) {
          L[i] <- 1e-8
          R[i] <- min(LL)
          delta[i] <- 0
        } else if (T[i] > max(RR)) {
          L[i] <- max(RR)
          R[i] <- 1e8
          delta[i] <- 0
        } else {
          idd <- which(T[i] > LL & T[i] < RR)
          if (length(idd) == 1) {
            L[i] <- LL[idd]
            R[i] <- RR[idd]
            delta[i] <- 0
          } else {
            L[i] <- T[i]
            R[i] <- T[i]
            delta[i] <- 1
          }
        }
      }
    }

    data.frame(
      L = L,
      R = R,
      delta = delta,
      id = id,
      x1 = x[, 1],
      x2 = x[, 2]
    )
  }

  dt <- make_pic_data()
  fit <- aftpen_pic(dt, lambda = 0.1, se = "CF", type = "BAR")

  expect_true(is.list(fit))
  expect_true("beta" %in% names(fit))
  expect_equal(length(fit$beta), 2)
  expect_true(all(is.finite(fit$beta)))
})
