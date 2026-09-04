# Prototype of the JT recursion with pluggable quadrature rules; accuracy vs exact (Miwa) / high-resolution references.
options(digits = 10)
source("analysis/quad-lib.R")
mk <- function(x, a = NULL) list(I = x$n.I, a = if (is.null(a)) x$lower$bound else a, b = x$upper$bound, delta = x$delta)
designs <- list(
  "OF one-sided K=3"        = mk(gsDesign(k = 3, test.type = 1, sfu = sfLDOF), a = rep(-Inf, 3)),
  "Pocock one-sided K=5"    = mk(gsDesign(k = 5, test.type = 1, sfu = sfLDPocock), a = rep(-Inf, 5)),
  "asym K=3 tt3"            = mk(gsDesign(k = 3, test.type = 3)),
  "asym K=5 tt4"            = mk(gsDesign(k = 5, test.type = 4)),
  "sym K=4 OF"              = mk(gsDesign(k = 4, test.type = 2, sfu = "OF")),
  "unequal K=4 tt4"         = mk(gsDesign(k = 4, test.type = 4, timing = c(0.2, 0.5, 0.9))),
  "close final K=3 (.5,.98)"= mk(gsDesign(k = 3, test.type = 4, timing = c(0.5, 0.98))),
  "skip eff IA1/fut IA2 K=4"= { x <- gsDesign(k = 4, test.type = 4, testUpper = c(FALSE, TRUE, TRUE, TRUE), testLower = c(TRUE, FALSE, TRUE, TRUE)); a <- x$lower$bound; b <- x$upper$bound; a[a <= -20] <- -Inf; b[b >= 20] <- Inf; list(I = x$n.I, a = a, b = b, delta = x$delta) },
  "one-sided K=10 HSD"      = mk(gsDesign(k = 10, test.type = 1), a = rep(-Inf, 10)),
  "asym K=10 tt4"           = mk(gsDesign(k = 10, test.type = 4)),
  "asym K=20 tt4"           = mk(gsDesign(k = 20, test.type = 4)),
  "one-sided K=20 Pocock"   = mk(gsDesign(k = 20, test.type = 1, sfu = sfLDPocock), a = rep(-Inf, 20))
)
thetas <- function(d) c(0, d$delta, 2 * d$delta, 3 * d$delta, -d$delta)
# references: Miwa for K<=5, else agreement of JT r=80 and GL n=160 (report their discrepancy)
refs <- list(); cat("Reference cross-checks (K>5): |JT r=80 - GL n=160|\n")
for (nm in names(designs)) { d <- designs[[nm]]; for (th in thetas(d)) {
  if (length(d$I) <= 5) refs[[paste(nm, th)]] <- ref_miwa(d$I, d$a, d$b, th) else {
    p1 <- recurse(d$I, d$a, d$b, th, jt_method(80)); p2 <- recurse(d$I, d$a, d$b, th, gauss_grid(160))
    cat(sprintf("  %-26s theta=%5.2f  %.2e\n", nm, th, max(abs(p1$upper - p2$upper), abs(p1$lower - p2$lower))))
    refs[[paste(nm, th)]] <- list(upper = (p1$upper + p2$upper) / 2, lower = (p1$lower + p2$lower) / 2) } } }
methods <- list()
for (r in c(6, 9, 12, 18, 24, 36)) methods[[sprintf("JT-Simpson r=%d (%d pts)", r, 12 * r - 3)]] <- jt_method(r)
for (n in c(25, 37, 49, 61, 91, 121, 181)) methods[[sprintf("Uniform-Weddle n=%d", n)]] <- uniform_grid(n, "weddle")
for (n in c(25, 49, 91, 181)) methods[[sprintf("Uniform-Boole n=%d", n)]] <- uniform_grid(n, "boole")
for (n in c(49, 91, 181)) methods[[sprintf("Uniform-Simpson n=%d", n)]] <- uniform_grid(n, "simpson")
for (n in c(12, 16, 20, 24, 32, 40, 48, 64, 96)) methods[[sprintf("Gauss-Legendre n=%d", n)]] <- gauss_grid(n)
for (cc in list(c(8, 3), c(8, 4), c(12, 4), c(8, 6), c(12, 6))) methods[[sprintf("GL adaptive n=%d+%d*W/sig", cc[1], cc[2])]] <- gauss_adaptive(cc[1], cc[2])
out <- data.frame(); detail <- list()
for (mn in names(methods)) {
  errs <- c(); npts <- c()
  for (nm in names(designs)) { d <- designs[[nm]]; for (th in thetas(d)) {
    p <- recurse(d$I, d$a, d$b, th, methods[[mn]]); rf <- refs[[paste(nm, th)]]
    errs[paste(nm, th)] <- max(abs(p$upper - rf$upper), abs(p$lower - rf$lower)); npts[paste(nm, th)] <- mean(p$npts) } }
  detail[[mn]] <- errs
  out <- rbind(out, data.frame(method = mn, mean_pts = round(mean(npts)), max_err = max(errs), median_err = median(errs), worst_case = names(errs)[which.max(errs)]))
}
print(format(out, digits = 3), row.names = FALSE)
cat("\nPer-design errors for selected methods:\n")
sel <- c("JT-Simpson r=18 (213 pts)", "Uniform-Weddle n=91", "Gauss-Legendre n=24", "Gauss-Legendre n=32", "Gauss-Legendre n=48", "GL adaptive n=8+4*W/sig", "GL adaptive n=12+6*W/sig")
tab <- do.call(cbind, lapply(sel, function(m) detail[[m]])); colnames(tab) <- c("JT18", "Wed91", "GL24", "GL32", "GL48", "Ad8+4", "Ad12+6")
print(format(tab, digits = 2, scientific = TRUE), quote = FALSE)
saveRDS(list(summary = out, detail = detail), "analysis/results-03-quadrature.rds")
