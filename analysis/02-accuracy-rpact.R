suppressPackageStartupMessages({library(gsDesign); library(rpact); library(mvtnorm)})
options(digits = 10)
# Reference: exact multivariate normal probabilities via mvtnorm (Miwa, deterministic) for K<=5
ref_probs <- function(I, a, b, theta) {
  K <- length(I); corr <- outer(seq_len(K), seq_len(K), function(i, j) sqrt(pmin(I[i], I[j]) / pmax(I[i], I[j])))
  mu <- theta * sqrt(I)
  up <- lo <- numeric(K)
  for (k in seq_len(K)) {
    if (k == 1) { up[1] <- pnorm(b[1] - mu[1], lower.tail = FALSE); lo[1] <- pnorm(a[1] - mu[1]); next }
    idx <- 1:k
    up[k] <- pmvnorm(lower = c(a[1:(k-1)], b[k]), upper = c(b[1:(k-1)], Inf), mean = mu[idx], corr = corr[idx, idx], algorithm = Miwa(steps = 4096))
    lo[k] <- pmvnorm(lower = c(a[1:(k-1)], -Inf), upper = c(b[1:(k-1)], a[k]), mean = mu[idx], corr = corr[idx, idx], algorithm = Miwa(steps = 4096))
  }
  list(upper = as.numeric(up), lower = as.numeric(lo))
}
# rpact one-sided: decision matrix rows = (futility, efficacy) shifted by drift; rpact clamps at -6 / 8 internally
rpact_probs <- function(I, a, b, theta) {
  mu <- theta * sqrt(I); K <- length(I)
  fut <- a - mu; fut[a <= -6] <- -6            # rpact convention: floor stays unshifted when "no futility"
  dm <- matrix(c(fut, b - mu), nrow = 2, byrow = TRUE)
  p <- rpact:::.getGroupSequentialProbabilitiesCpp(dm, I / I[K])
  list(upper = p[3, ] - p[2, ], lower = p[1, ])
}
gs_probs <- function(I, a, b, theta, r) { z <- gsProbability(k = length(I), theta = theta, n.I = I, a = a, b = b, r = r); list(upper = as.numeric(z$upper$prob), lower = as.numeric(z$lower$prob)) }
designs <- list(
  "OF-like one-sided K=3 (no futility)" = { x <- gsDesign(k = 3, test.type = 1, sfu = sfLDOF); list(I = x$n.I, a = rep(-20, 3), b = x$upper$bound, delta = x$delta) },
  "Pocock-like one-sided K=5 (no futility)" = { x <- gsDesign(k = 5, test.type = 1, sfu = sfLDPocock); list(I = x$n.I, a = rep(-20, 5), b = x$upper$bound, delta = x$delta) },
  "default asymmetric K=3 binding futility" = { x <- gsDesign(k = 3, test.type = 3); list(I = x$n.I, a = x$lower$bound, b = x$upper$bound, delta = x$delta) },
  "asymmetric K=5 non-binding tt4" = { x <- gsDesign(k = 5, test.type = 4); list(I = x$n.I, a = x$lower$bound, b = x$upper$bound, delta = x$delta) },
  "symmetric two-sided K=4 OF" = { x <- gsDesign(k = 4, test.type = 2, sfu = "OF"); list(I = x$n.I, a = x$lower$bound, b = x$upper$bound, delta = x$delta) },
  "unequal timing K=4 tt4" = { x <- gsDesign(k = 4, test.type = 4, timing = c(0.2, 0.5, 0.9)); list(I = x$n.I, a = x$lower$bound, b = x$upper$bound, delta = x$delta) }
)
res <- list()
for (nm in names(designs)) {
  d <- designs[[nm]]
  for (th in c(0, d$delta, 2 * d$delta)) {
    ref <- ref_probs(d$I, d$a, d$b, th)
    err <- function(p) c(upper = max(abs(p$upper - ref$upper)), lower = max(abs(p$lower - ref$lower)))
    row <- c(gs_r6 = err(gs_probs(d$I, d$a, d$b, th, 6)), gs_r12 = err(gs_probs(d$I, d$a, d$b, th, 12)), gs_r18 = err(gs_probs(d$I, d$a, d$b, th, 18)),
             gs_r36 = err(gs_probs(d$I, d$a, d$b, th, 36)), gs_r80 = err(gs_probs(d$I, d$a, d$b, th, 80)), rpact = err(rpact_probs(d$I, d$a, d$b, th)))
    res[[length(res) + 1]] <- data.frame(design = nm, theta = round(th / max(d$delta, 1e-9), 1), t(row))
  }
}
res <- do.call(rbind, res)
print(format(res, digits = 3, scientific = TRUE), row.names = FALSE)
saveRDS(res, "analysis/results-02-accuracy.rds")
