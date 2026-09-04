# Shared prototype library: pluggable-quadrature JT recursion + exact reference (sourced by analysis scripts)
suppressPackageStartupMessages({library(gsDesign); library(mvtnorm)})
# ---- grids -------------------------------------------------------------------
jt_grid <- function(r, mu, a, b) {   # replicate JT / gsDesign gridpts (Simpson on JT points)
  x <- c(mu - 3 - 4 * log(r / (1:(r - 1))), mu - 3 + 3 * (0:(4 * r)) / 2 / r, mu + 3 + 4 * log(r / (r - 1):1))
  if (min(x) < a) x <- c(a, x[x > a]); if (max(x) > b) x <- c(x[x < b], b)
  m <- length(x)
  if (m < 2) return(list(z = x, w = 0))
  y <- (x[2:m] + x[1:(m - 1)]) / 2
  wodd <- if (m == 2) c(x[2] - x[1], x[2] - x[1]) / 6 else { i <- 2:(m - 1); c(x[2] - x[1], x[i + 1] - x[i - 1], x[m] - x[m - 1]) / 6 }
  weven <- 4 * (x[2:m] - x[1:(m - 1)]) / 6
  z <- w <- numeric(2 * m - 1); z[2 * (1:m) - 1] <- x; z[2 * (1:(m - 1))] <- y; w[2 * (1:m) - 1] <- wodd; w[2 * (1:(m - 1))] <- weven
  list(z = z, w = w)
}
jt_method <- function(r) { force(r); function(mu, a, b, sig) jt_grid(r, mu, a, b) }
uniform_grid <- function(n, rule, L = 8.5) { force(n); force(rule); force(L)
  p <- switch(rule, simpson = 2, boole = 4, weddle = 6)
  wp <- switch(rule, simpson = c(1, 4, 1) / 3, boole = c(7, 32, 12, 32, 7) * 2 / 45, weddle = c(41, 216, 27, 272, 27, 216, 41) / 140)
  function(mu, a, b, sig) {
    lo <- max(a, mu - L); hi <- min(b, mu + L); if (hi <= lo) return(list(z = lo, w = 0))
    npan <- max(1, round((n - 1) / p)); m <- npan * p + 1; h <- (hi - lo) / (m - 1)
    z <- lo + h * (0:(m - 1)); w <- numeric(m)
    for (j in 0:(npan - 1)) { idx <- j * p + 1:(p + 1); w[idx] <- w[idx] + wp }
    list(z = z, w = w * h) } }
gl_nodes <- function(n) { if (n == 1) return(list(x = 0, w = 2)); i <- 1:(n - 1); beta <- i / sqrt(4 * i^2 - 1); J <- matrix(0, n, n); J[cbind(i, i + 1)] <- beta; J[cbind(i + 1, i)] <- beta
  e <- eigen(J, symmetric = TRUE); list(x = e$values, w = 2 * e$vectors[1, ]^2) }
GL <- lapply(1:200, gl_nodes)
gauss_grid <- function(n, L = 8.5) { force(n); force(L); g <- GL[[n]]; function(mu, a, b, sig) { lo <- max(a, mu - L); hi <- min(b, mu + L)
  if (hi <= lo) return(list(z = lo, w = 0)); list(z = (hi + lo) / 2 + (hi - lo) / 2 * g$x, w = (hi - lo) / 2 * g$w) } }
# adaptive GL: n = clamp(ceil(c0 + c1 * width / sig), nmin, nmax), sig = kernel width sqrt(Delta_{k+1}/I_k)
gauss_adaptive <- function(c0, c1, nmin = 8, nmax = 200, L = 8.5) { force(c0); force(c1); force(nmin); force(nmax); force(L)
  function(mu, a, b, sig) { lo <- max(a, mu - L); hi <- min(b, mu + L); if (hi <= lo) return(list(z = lo, w = 0))
    n <- min(nmax, max(nmin, ceiling(c0 + c1 * (hi - lo) / sig))); g <- GL[[n]]
    list(z = (hi + lo) / 2 + (hi - lo) / 2 * g$x, w = (hi - lo) / 2 * g$w) } }
# ---- generic recursion --------------------------------------------------------
recurse <- function(I, a, b, theta, grid_fn) {
  K <- length(I); mu <- theta * sqrt(I); up <- lo <- numeric(K); npts <- integer(0)
  up[1] <- pnorm(b[1] - mu[1], lower.tail = FALSE); lo[1] <- pnorm(a[1] - mu[1])
  if (K == 1) return(list(upper = up, lower = lo, npts = npts))
  kw <- function(k) { # kernel width: min of incoming and outgoing sqrt(dI / I)
    inn <- if (k > 1) sqrt((I[k] - I[k - 1]) / I[k - 1]) else Inf; out <- if (k < K) sqrt((I[k + 1] - I[k]) / I[k]) else Inf; min(inn, out) }
  g <- grid_fn(mu[1], a[1], b[1], kw(1)); z <- g$z; h <- g$w * dnorm(z - mu[1]); npts <- length(z)
  for (k in 2:K) {
    D <- I[k] - I[k - 1]; s <- sqrt(I[k]); sp <- sqrt(I[k - 1]); sd <- sqrt(D)
    up[k] <- sum(h * pnorm((z * sp + theta * D - b[k] * s) / sd))
    lo[k] <- sum(h * pnorm((a[k] * s - z * sp - theta * D) / sd))
    if (k < K) { g2 <- grid_fn(mu[k], a[k], b[k], kw(k)); z2 <- g2$z; npts <- c(npts, length(z2))
      X <- outer(z2 * s, z * sp + theta * D, "-") / sd
      h <- g2$w * as.vector((dnorm(X) * s / sd) %*% h); z <- z2 }
  }
  list(upper = up, lower = lo, npts = npts)
}
ref_miwa <- function(I, a, b, theta) {
  K <- length(I); corr <- outer(seq_len(K), seq_len(K), function(i, j) sqrt(pmin(I[i], I[j]) / pmax(I[i], I[j]))); mu <- theta * sqrt(I)
  up <- lo <- numeric(K); up[1] <- pnorm(b[1] - mu[1], lower.tail = FALSE); lo[1] <- pnorm(a[1] - mu[1])
  for (k in 2:K) { idx <- 1:k
    up[k] <- pmvnorm(lower = c(a[1:(k-1)], b[k]), upper = c(b[1:(k-1)], Inf), mean = mu[idx], corr = corr[idx, idx], algorithm = Miwa(steps = 4096))
    lo[k] <- pmvnorm(lower = c(a[1:(k-1)], -Inf), upper = c(b[1:(k-1)], a[k]), mean = mu[idx], corr = corr[idx, idx], algorithm = Miwa(steps = 4096)) }
  list(upper = as.numeric(up), lower = as.numeric(lo))
}
