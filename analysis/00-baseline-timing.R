sink("analysis/results-00-baseline-timing.txt", split = TRUE)
suppressPackageStartupMessages({library(gsDesign); library(microbenchmark)})
cat("gsDesign", as.character(packageVersion("gsDesign")), "\n")
tm <- function(expr, times = 20) {
  e <- substitute(expr)
  mb <- microbenchmark(list = list(e), times = times)
  sprintf("%9.3f ms (median)", median(mb$time) / 1e6)
}
cat("gsDesign() default (k=3, tt=4):     ", tm(gsDesign()), "\n")
for (tt in 1:6) cat(sprintf("gsDesign(k=3, test.type=%d):         ", tt), tm(gsDesign(k = 3, test.type = tt)), "\n")
for (tt in c(1, 2, 3, 4)) cat(sprintf("gsDesign(k=10, test.type=%d):        ", tt), tm(gsDesign(k = 10, test.type = tt), times = 5), "\n")
cat("gsDesign(k=3, tt=1, sfu='OF'):       ", tm(gsDesign(k = 3, test.type = 1, sfu = "OF")), "\n")
cat("gsSurv() default:                    ", tm(gsSurv(), times = 5), "\n")
x <- gsDesign()
cat("gsProbability(d=x, 1 theta):         ", tm(gsProbability(d = x, theta = 0.1)), "\n")
cat("gsProbability(d=x, 50 theta):        ", tm(gsProbability(d = x, theta = seq(0, 0.5, length.out = 50))), "\n")
x10 <- gsDesign(k = 10)
cat("gsProbability(d=x10, 50 theta):      ", tm(gsProbability(d = x10, theta = seq(0, 0.5, length.out = 50))), "\n")
# raw .C timings
k <- 3L; r <- 18L
I <- x$n.I; a <- x$lower$bound; b <- x$upper$bound
plo <- double(k); phi <- double(k)
cat("raw .C probrej (k=3, 1 theta):       ", tm(.C("probrej", k, 1L, 0.2, I, a, b, plo, phi, r, PACKAGE = "gsDesign"), times = 200), "\n")
k10 <- 10L; I10 <- x10$n.I; a10 <- x10$lower$bound; b10 <- x10$upper$bound
plo10 <- double(k10); phi10 <- double(k10)
cat("raw .C probrej (k=10, 1 theta):      ", tm(.C("probrej", k10, 1L, 0.2, I10, a10, b10, plo10, phi10, r, PACKAGE = "gsDesign"), times = 100), "\n")
th50 <- seq(0, 0.5, length.out = 50)
cat("raw .C probrej (k=10, 50 theta):     ", tm(.C("probrej", k10, 50L, th50, I10, a10, b10, double(500), double(500), r, PACKAGE = "gsDesign"), times = 20), "\n")
fp <- x$upper$spend
cat("raw .C gsbound1 (k=3):               ", tm(.C("gsbound1", k, 0, x$timing, rep(-20, k), rep(0, k), rep(0, k), fp, 1e-6, r, 0L, 0L, PACKAGE = "gsDesign"), times = 200), "\n")
fp10 <- x10$upper$spend
cat("raw .C gsbound1 (k=10):              ", tm(.C("gsbound1", k10, 0, x10$timing, rep(-20, k10), rep(0, k10), rep(0, k10), fp10, 1e-6, r, 0L, 0L, PACKAGE = "gsDesign"), times = 100), "\n")
cat("raw .C gsbound (k=10):               ", tm(.C("gsbound", k10, x10$timing, rep(0, k10), rep(0, k10), fp10, fp10, 1e-6, r, 0L, 0L, PACKAGE = "gsDesign"), times = 100), "\n")
cat("raw .C probrej r=80 (k=10, 1 theta): ", tm(.C("probrej", k10, 1L, 0.2, I10, a10, b10, plo10, phi10, 80L, PACKAGE = "gsDesign"), times = 20), "\n")
# count internal calls
ns <- asNamespace("gsDesign")
cnt <- c(gsprob = 0L, gsBound1 = 0L, gsBound = 0L)
for (f in names(cnt)) trace(f, quote(cnt[FN] <<- cnt[FN] + 1L), where = ns, print = FALSE)
# replace FN placeholder trick: define separate tracers
untrace("gsprob", where = ns); untrace("gsBound1", where = ns); untrace("gsBound", where = ns)
trace("gsprob", quote(cnt["gsprob"] <<- cnt["gsprob"] + 1L), where = ns, print = FALSE)
trace("gsBound1", quote(cnt["gsBound1"] <<- cnt["gsBound1"] + 1L), where = ns, print = FALSE)
trace("gsBound", quote(cnt["gsBound"] <<- cnt["gsBound"] + 1L), where = ns, print = FALSE)
for (spec in list(list(k=3,tt=4), list(k=10,tt=4), list(k=3,tt=3), list(k=3,tt=1), list(k=3,tt=2), list(k=3,tt=6))) {
  cnt[] <- 0L; invisible(gsDesign(k = spec$k, test.type = spec$tt))
  cat(sprintf("calls in gsDesign(k=%d, tt=%d): gsprob=%d gsBound1=%d gsBound=%d\n", spec$k, spec$tt, cnt["gsprob"], cnt["gsBound1"], cnt["gsBound"]))
}
untrace("gsprob", where = ns); untrace("gsBound1", where = ns); untrace("gsBound", where = ns)
# profile share of .C
for (spec in list(list(k=3,tt=4), list(k=10,tt=3))) {
  Rprof(tmp <- tempfile(), interval = 0.0005)
  for (i in 1:20) gsDesign(k = spec$k, test.type = spec$tt)
  Rprof(NULL)
  p <- summaryRprof(tmp)$by.self
  cat(sprintf("\n--- self-time profile gsDesign(k=%d, tt=%d) ---\n", spec$k, spec$tt))
  print(head(p[order(-p$self.pct), c("self.time", "self.pct")], 8))
}

sink()
