## =========================================
## Simulated data generation for aftPenCDA
## =========================================

set.seed(429)

## =========================
## 1. Right-censored data
## =========================

n <- 100
p <- 10

beta0 <- c(1, 1, 1, rep(0, p - 3))
x <- matrix(rnorm(n * p), n, p)

T <- exp(x %*% beta0 + rnorm(n))
C <- rexp(n, rate = 0.5)

y <- pmin(T, C)
d <- as.numeric(T <= C)

simdat_rc <- data.frame(y = y, d = d, x)


## =========================
## 2. PIC data (clustered)
## =========================

n <- 100
p <- 2
beta0 <- c(1, 1)

clu_rate <- 0.5
exactrates <- 0.8
left <- 0.001
right <- 0.01

## cluster-level frailty
eta <- 1 / clu_rate
v <- rgamma(n, shape = eta, rate = eta)
m <- ifelse(v > median(v), 5, 3)
id <- rep(seq_len(n), m)
vi <- rep(v, m)

## subject-level
N <- sum(m)
x <- matrix(rnorm(N * p), ncol = p)
colnames(x) <- paste0("x", seq_len(p))

T <- as.vector(exp(x %*% beta0 + vi * log(rexp(N))))

## interval construction
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

simdat_pic <- data.frame(
  L = L,
  R = R,
  delta = delta,
  id = id,
  x1 = x[, 1],
  x2 = x[, 2]
)


## =========================
## 3. Save to package data
## =========================

usethis::use_data(simdat_rc, simdat_pic, overwrite = TRUE)