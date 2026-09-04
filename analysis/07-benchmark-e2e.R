# End-to-end timings of user-facing functions. Run in a fresh R process per configuration:
#   CONFIG=baseline  -> installed gsDesign 3.11.0
#   CONFIG=jt        -> work/gsDesign (patched), default quadrature (Jennison & Turnbull grid)
#   CONFIG=gl        -> work/gsDesign (patched), options(gsDesign.quadrature = "gl")
suppressPackageStartupMessages(library(microbenchmark))
cfg <- Sys.getenv("CONFIG", "baseline")
if (cfg == "baseline") suppressPackageStartupMessages(library(gsDesign)) else suppressPackageStartupMessages(devtools::load_all("work/gsDesign", quiet = TRUE))
if (cfg == "gl") options(gsDesign.quadrature = "gl")
tm <- function(expr, times = 20) { e <- substitute(expr); mb <- microbenchmark(list = list(e), times = times); median(mb$time) / 1e6 }
x3 <- gsDesign(k = 3); x10 <- gsDesign(k = 10)
th50 <- seq(0, 0.5, length.out = 50)
res <- data.frame(config = cfg, task = c(
  "gsDesign(k=3, test.type=1)", "gsDesign(k=3, test.type=2)", "gsDesign(k=3, test.type=3)", "gsDesign(k=3, test.type=4)", "gsDesign(k=3, test.type=5)", "gsDesign(k=3, test.type=6)",
  "gsDesign(k=10, test.type=1)", "gsDesign(k=10, test.type=3)", "gsDesign(k=10, test.type=4)", "gsDesign(k=5, test.type=1, sfu='OF')",
  "gsSurv() default", "gsSurv(k=5, test.type=3)",
  "gsProbability(k=3 design, 50 theta)", "gsProbability(k=10 design, 50 theta)", "gsBoundSummary(k=3 design)", "gsCP + gsPP (k=10 design)", "sequentialPValue (k=4)"),
  ms = c(
  tm(gsDesign(k = 3, test.type = 1)), tm(gsDesign(k = 3, test.type = 2)), tm(gsDesign(k = 3, test.type = 3)), tm(gsDesign(k = 3, test.type = 4)), tm(gsDesign(k = 3, test.type = 5)), tm(gsDesign(k = 3, test.type = 6)),
  tm(gsDesign(k = 10, test.type = 1), 10), tm(gsDesign(k = 10, test.type = 3), 5), tm(gsDesign(k = 10, test.type = 4), 10), tm(gsDesign(k = 5, test.type = 1, sfu = "OF")),
  tm(gsSurv(), 5), tm(gsSurv(k = 5, test.type = 3), 5),
  tm(gsProbability(d = x3, theta = th50)), tm(gsProbability(d = x10, theta = th50)), tm(gsBoundSummary(x3)), tm({gsCP(x10, i = 3, zi = 1); gsPP(x10, i = 3, zi = 1)}),
  tm(sequentialPValue(gsD = gsDesign(k = 4), n.I = c(100, 200, 300, 400), Z = c(1, 2, 2.5, 3), usTime = c(.25, .5, .75, 1)), 10)))
out <- sprintf("analysis/results-07-e2e-%s.rds", cfg); saveRDS(res, out); print(res, digits = 4); cat("saved", out, "\n")
