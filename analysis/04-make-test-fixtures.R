# Generate fixtures capturing gsDesign 3.11.0 behavior (installed package) plus exact mvtnorm references.
suppressPackageStartupMessages({library(mvtnorm)})
# Set GSDESIGN_SRC to a package source directory to generate fixtures from a development version.
src_dir <- Sys.getenv("GSDESIGN_SRC", unset = "")
if (nzchar(src_dir)) suppressPackageStartupMessages(devtools::load_all(src_dir, quiet = TRUE)) else suppressPackageStartupMessages(library(gsDesign))
out_file <- Sys.getenv("GSDESIGN_FIXTURE_OUT", unset = "analysis/numint-fixtures.R")
use_inf <- nzchar(Sys.getenv("GSDESIGN_FIXTURE_INF", unset = ""))  # use -Inf/Inf instead of -20/20 in inputs
NB <- if (use_inf) Inf else 20
ref_miwa <- function(I, a, b, theta) {
  K <- length(I); corr <- outer(seq_len(K), seq_len(K), function(i, j) sqrt(pmin(I[i], I[j]) / pmax(I[i], I[j]))); mu <- theta * sqrt(I)
  up <- lo <- numeric(K); up[1] <- pnorm(b[1] - mu[1], lower.tail = FALSE); lo[1] <- pnorm(a[1] - mu[1])
  for (k in seq_len(K)[-1]) { idx <- 1:k
    up[k] <- pmvnorm(lower = c(a[1:(k-1)], b[k]), upper = c(b[1:(k-1)], Inf), mean = mu[idx], corr = corr[idx, idx], algorithm = Miwa(steps = 4096))
    lo[k] <- pmvnorm(lower = c(a[1:(k-1)], -Inf), upper = c(b[1:(k-1)], a[k]), mean = mu[idx], corr = corr[idx, idx], algorithm = Miwa(steps = 4096)) }
  list(upper = as.numeric(up), lower = as.numeric(lo))
}
sig <- function(x) signif(x, 13)
fx <- list()
# --- A. gsDesign() outputs across test types and boundary families ---------------------------
calls <- list(
  tt1_default = quote(gsDesign(k = 3, test.type = 1)),
  tt2_default = quote(gsDesign(k = 3, test.type = 2)),
  tt3_default = quote(gsDesign(k = 3, test.type = 3)),
  tt4_default = quote(gsDesign(k = 3, test.type = 4)),
  tt5_default = quote(gsDesign(k = 3, test.type = 5)),
  tt6_default = quote(gsDesign(k = 3, test.type = 6)),
  tt7_default = quote(gsDesign(k = 3, test.type = 7)),
  tt8_default = quote(gsDesign(k = 3, test.type = 8)),
  tt1_OF = quote(gsDesign(k = 4, test.type = 1, sfu = "OF")),
  tt2_Pocock = quote(gsDesign(k = 4, test.type = 2, sfu = "Pocock")),
  tt2_WT = quote(gsDesign(k = 4, test.type = 2, sfu = "WT", sfupar = 0.25)),
  tt4_k5_LDOF = quote(gsDesign(k = 5, test.type = 4, sfu = sfLDOF, sfl = sfLDOF)),
  tt4_k10 = quote(gsDesign(k = 10, test.type = 4)),
  tt3_timing = quote(gsDesign(k = 4, test.type = 3, timing = c(0.25, 0.5, 0.75), sfl = sfHSD, sflpar = -1)),
  tt4_nfix = quote(gsDesign(k = 3, test.type = 4, n.fix = 800, delta = 0)),
  tt4_delta = quote(gsDesign(k = 3, test.type = 4, delta = 0.25)),
  tt4_nI = quote(gsDesign(k = 3, test.type = 4, n.I = c(300, 600, 860), maxn.IPlan = 900)),
  tt2_nI = quote(gsDesign(k = 3, test.type = 2, n.fix = 800, n.I = c(300, 600, 860), maxn.IPlan = 900)),
  tt4_usTime = quote(gsDesign(k = 3, test.type = 4, usTime = c(0.3, 0.7), lsTime = c(0.4, 0.8))),
  tt4_sfPoints = quote(gsDesign(k = 4, test.type = 4, timing = c(.1, .4, .7), sfu = sfPoints, sfupar = c(.033333, .063367, .1), sfl = sfPoints, sflpar = c(.25, .5, .75))),
  tt1_sfTruncated = quote(gsDesign(k = 4, test.type = 1, sfu = sfTruncated, sfupar = list(sf = sfHSD, param = -4, trange = c(0.4, 1)))),
  tt4_sfGapped = quote(gsDesign(k = 4, test.type = 4, sfu = sfGapped, sfupar = list(sf = sfHSD, param = -4, trange = c(0.3, 0.8)))),
  tt4_skipUpper = quote(gsDesign(k = 3, test.type = 4, testUpper = c(FALSE, TRUE, TRUE))),
  tt4_skipLower = quote(gsDesign(k = 3, test.type = 4, testLower = c(TRUE, FALSE, FALSE))),
  tt3_skipLower = quote(gsDesign(k = 3, test.type = 3, testLower = c(FALSE, TRUE, FALSE))),
  tt8_skipHarm = quote(gsDesign(k = 3, test.type = 8, astar = 0.05, testHarm = c(TRUE, TRUE, FALSE))),
  tt7_skipMixed = quote(gsDesign(k = 3, test.type = 7, astar = 0.05, testLower = c(TRUE, FALSE, TRUE), testHarm = c(TRUE, TRUE, FALSE))),
  tt4_r6 = quote(gsDesign(k = 3, test.type = 4, r = 6)),
  tt4_r80 = quote(gsDesign(k = 3, test.type = 4, r = 80)),
  tt4_alpha05 = quote(gsDesign(k = 4, test.type = 4, alpha = 0.05, beta = 0.2, sfu = sfPower, sfupar = 2, sfl = sfHSD, sflpar = 1)),
  tt6_astar = quote(gsDesign(k = 3, test.type = 6, astar = 0.1)),
  tt5_astar = quote(gsDesign(k = 3, test.type = 5, astar = 0.1))
)
fx$gsDesign <- lapply(calls, function(cl) {
  x <- eval(cl)
  out <- list(call = deparse(cl), k = x$k, n.I = sig(x$n.I), upper = sig(x$upper$bound), upper_prob = sig(x$upper$prob), en = sig(x$en), theta = sig(x$theta), delta = sig(x$delta))
  if (!is.null(x$lower)) { out$lower <- sig(x$lower$bound); out$lower_prob <- sig(x$lower$prob) }
  if (!is.null(x$harm)) { out$harm <- sig(x$harm$bound); out$harm_prob <- sig(x$harm$prob) }
  if (!is.null(x$falseposnb)) out$falseposnb <- sig(x$falseposnb)
  out
})
# --- B. gsProbability at diverse theta incl. large drift & infinite-ish bounds ---------------
x4 <- gsDesign(k = 4, test.type = 4); x1 <- gsDesign(k = 5, test.type = 1, sfu = sfLDOF)
probs <- list(
  tt4_theta = list(k = 4, theta = x4$delta * c(-1, 0, 0.5, 1, 2, 3), n.I = x4$n.I, a = x4$lower$bound, b = x4$upper$bound),
  tt1_theta = list(k = 5, theta = x1$delta * c(0, 1, 2, 4), n.I = x1$n.I, a = rep(-NB, 5), b = x1$upper$bound),
  skip = list(k = 4, theta = x4$delta * c(0, 1, 3), n.I = x4$n.I, a = c(x4$lower$bound[1], -NB, x4$lower$bound[3:4]), b = c(NB, x4$upper$bound[2:4])),
  r6 = list(k = 4, theta = x4$delta * c(0, 1), n.I = x4$n.I, a = x4$lower$bound, b = x4$upper$bound, r = 6),
  r80 = list(k = 4, theta = x4$delta * c(0, 1), n.I = x4$n.I, a = x4$lower$bound, b = x4$upper$bound, r = 80),
  overrun = list(k = 4, theta = x4$delta * c(0, 1), n.I = x4$n.I, a = x4$lower$bound, b = x4$upper$bound, overrun = c(10, 20, 30)),
  k1 = list(k = 1, theta = c(0, 1, 2), n.I = 1, a = -NB, b = 1.96),
  crossed = list(k = 3, theta = c(0, 1), n.I = 1:3, a = c(0, 2.5, 1.9), b = c(2.5, 2.2, 1.9)),
  close = list(k = 3, theta = c(0, 3.24), n.I = c(0.5, 0.98, 1), a = c(-0.5, 1.0, 2.0), b = c(3, 2.5, 2.0))
)
fx$gsProbability <- lapply(probs, function(p) { z <- do.call(gsProbability, p); list(args = p, upper = sig(z$upper$prob), lower = sig(z$lower$prob), en = sig(z$en)) })
# --- C. gsBound / gsBound1 direct calls (incl. zero spend -> +/-20, error flags) --------------
bcalls <- list(
  gsBound_basic = quote(gsBound(I = c(1, 2, 3) / 3, trueneg = rep(.02, 3), falsepos = rep(.01, 3))),
  gsBound_zero = quote(gsBound(I = c(1, 2, 3), trueneg = c(0.01, 0, 0.02), falsepos = c(0.001, 0, 0.024))),
  gsBound_zero_close = quote(gsBound(I = c(1, 1.01, 3), trueneg = c(0.01, 0, 0.02), falsepos = c(0.001, 0, 0.024))),
  gsBound_r = quote(gsBound(I = c(1, 2, 3) / 3, trueneg = rep(.02, 3), falsepos = rep(.01, 3), r = 6, tol = 1e-8)),
  gsBound1_basic = quote(gsBound1(theta = 0, I = c(1, 2, 3) / 3, a = rep(-20, 3), probhi = c(.001, .009, .015))),
  gsBound1_theta = quote(gsBound1(theta = -2.2, I = c(1, 2, 3) / 3, a = -c(3.0902, 2.5105, 2.0011), probhi = rep(.05, 3))),
  gsBound1_zero = quote(gsBound1(theta = 0, I = c(1, 2, 3), a = rep(-20, 3), probhi = c(0.001, 0, 0.024))),
  gsBound1_zero_close = quote(gsBound1(theta = 0, I = c(1, 1.01, 3), a = rep(-20, 3), probhi = c(0.001, 0, 0.024))),
  gsBound1_zero_first = quote(gsBound1(theta = 0, I = c(1, 2, 3), a = rep(-20, 3), probhi = c(0, 0.01, 0.015))),
  gsBound1_bigdrift = quote(gsBound1(theta = 5, I = c(1, 2, 3), a = rep(-20, 3), probhi = c(0, 0.5, 0.4))),
  gsBound1_crossing = quote(gsBound1(theta = 0, I = c(1, 2, 3), a = c(-1, 1, 2.5), probhi = c(0.01, 0.02, 0.1))),
  gsBound1_k1 = quote(gsBound1(theta = 0.5, I = 2, a = -20, probhi = 0.025))
)
if (use_inf) bcalls <- lapply(bcalls, function(cl) parse(text = gsub("-20", "-Inf", paste(deparse(cl), collapse = "")))[[1]])
fx$gsBound <- lapply(bcalls, function(cl) { y <- eval(cl); out <- list(call = deparse(cl), a = sig(y$a), b = sig(y$b), error = y$error); if (!is.null(y$problo)) out$problo <- sig(y$problo); out })
# --- D. gsDensity, normalGrid, gsCP/gsPP/gsPI/gsPOS/gsCPOS/gsBoundCP, sequentialPValue -------
xd <- gsDesign(k = 4, test.type = 4)
fx$gsDensity <- list(
  i2 = { d <- gsDensity(xd, theta = c(0, xd$delta), i = 2, zi = seq(-3, 4, 0.5)); list(theta = sig(d$theta), zi = d$zi, density = sig(d$density)) },
  i4_r10 = { d <- gsDensity(xd, theta = c(0, xd$delta, 2 * xd$delta), i = 4, zi = c(-2, 0, 1.5, 2.5, 3.5), r = 10); list(density = sig(d$density)) },
  i1 = { d <- gsDensity(xd, theta = 0.1, i = 1, zi = c(-1, 0, 1)); list(density = sig(d$density)) }
)
fx$normalGrid <- list(
  default = { g <- normalGrid(); list(z = sig(g$z), gridwgts = sig(g$gridwgts), wgts_sum = sig(sum(g$wgts))) },
  r3 = { g <- normalGrid(r = 3, mu = 2, sigma = 3); list(z = sig(g$z), gridwgts = sig(g$gridwgts)) },
  bounds = { g <- normalGrid(bounds = c(-1, 1)); list(n = length(g$z), wgts_sum = sig(sum(g$wgts))) }
)
fx$gsCP <- list(
  gsCP = { y <- gsCP(xd, i = 2, zi = 1.2, theta = c(0, xd$delta)); list(upper = sig(y$upper$prob), lower = sig(y$lower$prob), n.I = sig(y$n.I)) },
  gsPP = sig(gsPP(xd, i = 2, zi = 1.2)),
  gsPI = sig(gsPI(xd, i = 2, zi = 1.2, j = 3)),
  gsPOS = sig(gsPOS(xd, theta = c(0, xd$delta), wgts = c(.5, .5))),
  gsCPOS = sig(gsCPOS(i = 2, xd, theta = c(0, xd$delta), wgts = c(.5, .5))),
  gsBoundCP = sig(as.matrix(gsBoundCP(xd))),
  gsBoundCP1 = sig(gsBoundCP(gsDesign(k = 3, test.type = 1)))
)
fx$sequentialPValue <- sig(sequentialPValue(gsD = xd, n.I = c(100, 200, 300, 400), Z = c(1, 2, 2.5, 3), usTime = c(.25, .5, .75, 1)))
# --- E. exact references (mvtnorm Miwa) for gsProbability correctness tests -----------------
exact_cases <- list(
  tt1_OF = { x <- gsDesign(k = 3, test.type = 1, sfu = sfLDOF); list(k = 3, n.I = x$n.I, a = rep(-NB, 3), b = x$upper$bound, theta = c(0, x$delta, 2 * x$delta)) },
  tt3_k3 = { x <- gsDesign(k = 3, test.type = 3); list(k = 3, n.I = x$n.I, a = x$lower$bound, b = x$upper$bound, theta = c(0, x$delta, 2 * x$delta)) },
  tt4_k5 = { x <- gsDesign(k = 5, test.type = 4); list(k = 5, n.I = x$n.I, a = x$lower$bound, b = x$upper$bound, theta = c(0, x$delta, 2 * x$delta)) },
  tt2_k4_OF = { x <- gsDesign(k = 4, test.type = 2, sfu = "OF"); list(k = 4, n.I = x$n.I, a = x$lower$bound, b = x$upper$bound, theta = c(0, x$delta)) },
  unequal_k4 = { x <- gsDesign(k = 4, test.type = 4, timing = c(0.2, 0.5, 0.9)); list(k = 4, n.I = x$n.I, a = x$lower$bound, b = x$upper$bound, theta = c(0, x$delta, 2 * x$delta)) },
  close_k3 = { x <- gsDesign(k = 3, test.type = 4, timing = c(0.5, 0.98)); list(k = 3, n.I = x$n.I, a = x$lower$bound, b = x$upper$bound, theta = c(0, x$delta)) }
)
fx$exact <- lapply(exact_cases, function(p) { rf <- lapply(p$theta, function(th) ref_miwa(p$n.I, p$a, p$b, th))
  list(args = p, upper = signif(do.call(cbind, lapply(rf, `[[`, "upper")), 12), lower = signif(do.call(cbind, lapply(rf, `[[`, "lower")), 12)) })
fx$exact_note <- "Reference probabilities computed with mvtnorm 1.4.2 pmvnorm(algorithm = Miwa(steps = 4096)); accuracy ~1e-10."
con <- file(out_file, "w")
writeLines(c(sprintf("# Generated by analysis/04-make-test-fixtures.R from gsDesign %s and mvtnorm %s.", as.character(packageVersion("gsDesign")), as.character(packageVersion("mvtnorm"))), "# Do not edit by hand.", "numint_fixtures <- "), con)
dput(fx, con, control = c("keepNA", "keepInteger", "niceNames", "showAttributes", "digits17"))
close(con)
cat("fixture file:", out_file, file.size(out_file), "bytes\n")
