suppressPackageStartupMessages(library(gsDesign))
options(digits = 17)
# 1. Does the current C code accept +/-Inf bounds in probrej?
x <- gsDesign(k = 3, test.type = 1)
p20  <- gsProbability(k = 3, theta = c(0, x$delta, 3 * x$delta), n.I = x$n.I, a = rep(-20, 3), b = x$upper$bound)
pInf <- gsProbability(k = 3, theta = c(0, x$delta, 3 * x$delta), n.I = x$n.I, a = rep(-Inf, 3), b = x$upper$bound)
cat("probrej: max |p20 - pInf| upper:", max(abs(p20$upper$prob - pInf$upper$prob)), " lower:", max(abs(p20$lower$prob - pInf$lower$prob)), "\n")
# with an efficacy analysis skipped (b = 20 vs Inf) and big drift
b <- c(20, x$upper$bound[2:3]); bI <- c(Inf, x$upper$bound[2:3])
for (th in c(0, x$delta, 2 * x$delta, 3 * x$delta, 4 * x$delta)) {
  q20 <- gsProbability(k = 3, theta = th, n.I = x$n.I, a = rep(-20, 3), b = b)
  qI  <- gsProbability(k = 3, theta = th, n.I = x$n.I, a = rep(-Inf, 3), b = bI)
  cat(sprintf("theta=%.3f mu_k=%.2f  upper diff: %.3e  lower diff: %.3e  P(b1=20 crossed)=%.3e\n", th, th * sqrt(x$n.I[3]),
      max(abs(q20$upper$prob - qI$upper$prob)), max(abs(q20$lower$prob - qI$lower$prob)), q20$upper$prob[1]))
}
# 2. gsBound1 with probhi = 0 at interim and close timings -> NaN?
for (I in list(c(1, 2, 3), c(1, 1.01, 3), c(1, 1.0001, 3), c(1, 1.000001, 3))) {
  y <- gsBound1(theta = 0, I = I, a = rep(-20, 3), probhi = c(0.001, 0, 0.024))
  cat("gsBound1 I=", format(I), " b=", format(y$b), " error=", y$error, "\n")
}
# 3. gsBound with trueneg = 0 / falsepos = 0 at interim & close timing
for (I in list(c(1, 2, 3), c(1, 1.01, 3), c(1, 1.0001, 3), c(1, 1.000001, 3))) {
  y <- gsBound(I = I, trueneg = c(0.01, 0, 0.02), falsepos = c(0.001, 0, 0.024))
  cat("gsBound  I=", format(I), " a=", format(y$a), " b=", format(y$b), " error=", y$error, "\n")
}
# 4. gsBound1 with a = -Inf passed in (currently allowed?)
y1 <- gsBound1(theta = 0, I = c(1, 2, 3), a = rep(-Inf, 3), probhi = c(0.001, 0.009, 0.015))
y2 <- gsBound1(theta = 0, I = c(1, 2, 3), a = rep(-20, 3), probhi = c(0.001, 0.009, 0.015))
cat("gsBound1 a=-Inf vs -20: b diff", max(abs(y1$b - y2$b)), " problo(Inf)=", y1$problo, "\n")
# 5. large drift with gsBound1: mu = theta*sqrt(I) large, probhi = 0 at analysis 1 -> b[0]=20 truncates grid?
for (theta in c(1, 3, 5, 8, 12)) {
  y <- gsBound1(theta = theta, I = c(1, 2, 3), a = rep(-Inf, 3), probhi = c(0, 0.5, 0.4))
  z <- gsProbability(k = 3, theta = theta, n.I = c(1, 2, 3), a = rep(-Inf, 3), b = y$b)
  cat(sprintf("theta=%5.1f  b=%s  achieved probhi=%s  err=%d\n", theta, paste(format(y$b, digits = 8), collapse = ","), paste(format(z$upper$prob, digits = 8), collapse = ","), y$error))
}
# 6. gsDesign with testUpper skipping and a large theta in gsProbability
d <- gsDesign(k = 4, test.type = 4, testUpper = c(FALSE, FALSE, TRUE, TRUE), testLower = c(TRUE, FALSE, TRUE, TRUE))
cat("bounds upper:", d$upper$bound, " lower:", d$lower$bound, "\n")
