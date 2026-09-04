# Where does replacing the finite sentinel 20 by Inf change results? Compare gsProbability with +/-20 vs +/-Inf
# for absent bounds, as a function of drift mu_K = theta * sqrt(I_K) and of the grid parameter r.
suppressPackageStartupMessages(devtools::load_all("work/gsDesign", quiet = TRUE))
x <- gsDesign(k = 4, test.type = 4)
res <- data.frame()
for (r in c(18, 40, 80)) for (muK in c(0, 2, 4, 5, 6, 8, 10, 12, 15)) {
  th <- muK / sqrt(x$n.I[4])
  a20 <- c(x$lower$bound[1], -20, x$lower$bound[3:4]); b20 <- c(20, x$upper$bound[2:4])
  aI <- a20; aI[2] <- -Inf; bI <- b20; bI[1] <- Inf
  p20 <- gsProbability(k = 4, theta = th, n.I = x$n.I, a = a20, b = b20, r = r)
  pI  <- gsProbability(k = 4, theta = th, n.I = x$n.I, a = aI, b = bI, r = r)
  one20 <- gsProbability(k = 4, theta = th, n.I = x$n.I, a = rep(-20, 4), b = x$upper$bound, r = r)
  oneI  <- gsProbability(k = 4, theta = th, n.I = x$n.I, a = rep(-Inf, 4), b = x$upper$bound, r = r)
  res <- rbind(res, data.frame(r = r, muK = muK, skip_upper_diff = max(abs(p20$upper$prob - pI$upper$prob)), skip_lower_diff = max(abs(p20$lower$prob - pI$lower$prob)),
                               onesided_diff = max(abs(one20$upper$prob - oneI$upper$prob)), grid_extent = 3 + 4 * log(r)))
}
print(format(res, digits = 3), row.names = FALSE)
# bound derivation: gsBound1 with zero spend at analysis 1 under drift (theta), 20 vs Inf sentinel effect on later bounds
cat("\ngsBound1 with probhi[1] = 0 under drift: bounds b[2:3] with old (20) vs new (Inf) behaviour, via gsProbability check\n")
for (theta in c(0, 2, 4, 6, 8)) {
  y <- gsBound1(theta = theta, I = c(1, 2, 3), a = rep(-Inf, 3), probhi = c(0, 0.3, 0.3))
  # emulate old behaviour: truncate grid at 20 by passing b[1] = 20 to gsProbability evaluation of the new bounds
  pI <- gsProbability(k = 3, theta = theta, n.I = 1:3, a = rep(-Inf, 3), b = y$b)
  p20 <- gsProbability(k = 3, theta = theta, n.I = 1:3, a = rep(-20, 3), b = c(20, y$b[2:3]))
  cat(sprintf("theta=%.0f (mu_1=%.1f)  b=%s  achieved(Inf)=%s  achieved(20)=%s\n", theta, theta, paste(format(y$b, digits = 7), collapse = ","), paste(format(pI$upper$prob, digits = 7), collapse = ","), paste(format(p20$upper$prob, digits = 7), collapse = ",")))
}
saveRDS(res, "analysis/results-06-sentinel.rds")
