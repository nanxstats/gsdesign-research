suppressPackageStartupMessages({library(gsDesign); library(rpact); library(microbenchmark)})
rows <- list()
tm <- function(expr, times = 200) { e <- substitute(expr); mb <- microbenchmark(list = list(e), times = times); median(mb$time) / 1e3 }
for (K in c(3, 5, 10)) {
  x <- gsDesign(k = K, test.type = 4); x1 <- gsDesign(k = K, test.type = 1, sfu = sfLDOF)
  I <- x$n.I / x$n.I[K]
  dm4 <- matrix(c(x$lower$bound, x$upper$bound), nrow = 2, byrow = TRUE)
  dm1 <- matrix(c(rep(-6, K), x1$upper$bound), nrow = 2, byrow = TRUE)
  t_rp4 <- tm(rpact:::.getGroupSequentialProbabilitiesCpp(dm4, I)); t_rp1 <- tm(rpact:::.getGroupSequentialProbabilitiesCpp(dm1, I))
  t_gs4 <- tm(.C("probrej", as.integer(K), 1L, 0, I, x$lower$bound, x$upper$bound, double(K), double(K), 18L, PACKAGE = "gsDesign"))
  t_gs1 <- tm(.C("probrej", as.integer(K), 1L, 0, I, rep(-20, K), x1$upper$bound, double(K), double(K), 18L, PACKAGE = "gsDesign"))
  cat(sprintf("K=%2d  crossing probs (one theta): rpact two-sided %7.1f us | gsDesign probrej %7.1f us || rpact one-sided %7.1f us | gsDesign probrej one-sided %7.1f us\n", K, t_rp4, t_gs4, t_rp1, t_gs1))
  rows[[length(rows) + 1]] <- data.frame(task = "crossing probabilities", K = K, rpact_two_sided_us = t_rp4, gsDesign_two_sided_us = t_gs4, rpact_one_sided_us = t_rp1, gsDesign_one_sided_us = t_gs1)
}
# design derivation: alpha-spending bounds (one-sided), OF-like spending, K = 3, 5, 10
for (K in c(3, 5, 10)) {
  t_rp <- tm(getDesignGroupSequential(kMax = K, alpha = 0.025, sided = 1, typeOfDesign = "asOF"), times = 20)
  t_rpc <- tm(rpact:::.getDesignGroupSequentialAlphaSpendingCpp(K, 0.025, NA_real_, "asOF", 1, (1:K) / K, FALSE, rep(-6, K), 1e-8), times = 50)
  t_gs <- tm(gsDesign(k = K, test.type = 1, sfu = sfLDOF), times = 20)
  fp <- gsDesign(k = K, test.type = 1, sfu = sfLDOF)$upper$spend
  t_gsc <- tm(.C("gsbound1", as.integer(K), 0, (1:K) / K, rep(-20, K), double(K), double(K), fp, 1e-6, 18L, 0L, 0L, PACKAGE = "gsDesign"), times = 200)
  cat(sprintf("K=%2d  OF-spending bounds: rpact getDesignGroupSequential %8.1f us (C++ core %7.1f us) | gsDesign() %8.1f us (C core gsbound1 %7.1f us)\n", K, t_rp, t_rpc, t_gs, t_gsc))
  rows[[length(rows) + 1]] <- data.frame(task = "OF-spending bounds", K = K, rpact_R_us = t_rp, rpact_cpp_us = t_rpc, gsDesign_R_us = t_gs, gsDesign_c_us = t_gsc)
}
# bound agreement rpact vs gsDesign, OF spending K=5
d <- getDesignGroupSequential(kMax = 5, alpha = 0.025, sided = 1, typeOfDesign = "asOF"); g <- gsDesign(k = 5, test.type = 1, sfu = sfLDOF)
cat("bounds rpact - gsDesign (asOF K=5):", format(d$criticalValues - g$upper$bound, digits = 3), "\n")
saveRDS(list(timing = rows, bound_diff = d$criticalValues - g$upper$bound), "analysis/results-05-rpact-timing.rds")
