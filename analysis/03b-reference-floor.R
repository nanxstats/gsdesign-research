suppressPackageStartupMessages({library(gsDesign)})
source("analysis/quad-lib.R")
# tighter references and stress designs
mk <- function(x, a = NULL) list(I = x$n.I, a = if (is.null(a)) x$lower$bound else a, b = x$upper$bound, delta = x$delta)
stress <- list(
  "one-sided K=10 HSD"      = mk(gsDesign(k = 10, test.type = 1), a = rep(-Inf, 10)),
  "asym K=20 tt4"           = mk(gsDesign(k = 20, test.type = 4)),
  "one-sided K=20 Pocock"   = mk(gsDesign(k = 20, test.type = 1, sfu = sfLDPocock), a = rep(-Inf, 20)),
  "K=50 constant bounds"    = list(I = (1:50) / 50, a = c(rep(-1, 49), 2.5), b = c(rep(3.2, 49), 2.5), delta = 3.24),
  "K=3 (.5,.995,1) dI/I=.005" = { x <- gsDesign(k = 3, test.type = 4, timing = c(0.5, 0.9)); list(I = x$n.I[3] * c(0.5, 0.995, 1), a = x$lower$bound, b = x$upper$bound, delta = x$delta) },
  "K=4 (.1,.2,.99,1)"       = { x <- gsDesign(k = 4, test.type = 4, timing = c(0.1, 0.2, 0.9)); list(I = x$n.I[4] * c(0.1, 0.2, 0.99, 1), a = x$lower$bound, b = x$upper$bound, delta = x$delta) },
  "no-stop IA2 (a=-Inf,b=Inf) K=4" = { x <- gsDesign(k = 4, test.type = 4); a <- x$lower$bound; b <- x$upper$bound; a[2] <- -Inf; b[2] <- Inf; list(I = x$n.I, a = a, b = b, delta = x$delta) },
  "asym K=5 tt4"            = mk(gsDesign(k = 5, test.type = 4)),
  "narrow-then-wide (1,1.01,3)" = { x <- gsDesign(k = 3, test.type = 4, timing = c(0.5, 0.9)); list(I = x$n.I[3] * c(1, 1.01, 3) / 3, a = c(x$lower$bound[1], -Inf, x$lower$bound[3]), b = c(x$upper$bound[1], Inf, x$upper$bound[3]), delta = x$delta) },
  "narrow-then-wide K=4 (.3,.303,.6,1)" = { x <- gsDesign(k = 4, test.type = 4, timing = c(0.3, 0.6, 0.8)); list(I = x$n.I[4] * c(0.3, 0.303, 0.6, 1), a = x$lower$bound, b = x$upper$bound, delta = x$delta) }
)
floor_tab <- data.frame(); err_tab <- data.frame()
cat("Reference floor: pairwise max |diff| over theta in {0, d, 2d, 5d, -d}\n")
cat(sprintf("%-32s %12s %12s %12s %12s\n", "design", "GL200-GL160", "GL200-JT160", "GL200-Wed721", "n_adapt(8+3)"))
refs <- list()
for (nm in names(stress)) { d <- stress[[nm]]; e1 <- e2 <- e3 <- 0; nn <- c()
  for (th in c(0, d$delta, 2 * d$delta, 5 * d$delta, -d$delta)) {
    p200 <- recurse(d$I, d$a, d$b, th, gauss_grid(200)); p160 <- recurse(d$I, d$a, d$b, th, gauss_grid(160))
    pjt <- recurse(d$I, d$a, d$b, th, jt_method(160)); pw <- recurse(d$I, d$a, d$b, th, uniform_grid(721, "weddle"))
    pa <- recurse(d$I, d$a, d$b, th, gauss_adaptive(8, 3)); nn <- c(nn, pa$npts)
    dd <- function(p, q) max(abs(p$upper - q$upper), abs(p$lower - q$lower))
    e1 <- max(e1, dd(p200, p160)); e2 <- max(e2, dd(p200, pjt)); e3 <- max(e3, dd(p200, pw))
    refs[[paste(nm, th)]] <- p200 }
  cat(sprintf("%-32s %12.2e %12.2e %12.2e %12s\n", nm, e1, e2, e3, paste(range(nn), collapse = "-")))
  floor_tab <- rbind(floor_tab, data.frame(design = nm, gl200_gl160 = e1, gl200_jt160 = e2, gl200_wed721 = e3, n_adapt_min = min(nn), n_adapt_max = max(nn))) }
cat("\nErrors vs GL200 reference:\n")
methods <- list("JT r=18" = jt_method(18), "JT r=36" = jt_method(36), "JT r=80" = jt_method(80), "Weddle n=91" = uniform_grid(91, "weddle"), "Weddle n=181" = uniform_grid(181, "weddle"),
  "GL n=32" = gauss_grid(32), "GL n=48" = gauss_grid(48), "GL n=64" = gauss_grid(64), "GLad 8+3W/s" = gauss_adaptive(8, 3), "GLad 6+2.5W/s" = gauss_adaptive(6, 2.5), "GLad 8+2W/s" = gauss_adaptive(8, 2), "GLad 10+4W/s" = gauss_adaptive(10, 4))
cat(sprintf("%-32s", "design")); for (m in names(methods)) cat(sprintf(" %13s", m)); cat("\n")
for (nm in names(stress)) { d <- stress[[nm]]; cat(sprintf("%-32s", nm))
  for (m in names(methods)) { e <- 0; np <- c()
    for (th in c(0, d$delta, 2 * d$delta, 5 * d$delta, -d$delta)) { p <- recurse(d$I, d$a, d$b, th, methods[[m]]); rf <- refs[[paste(nm, th)]]; e <- max(e, max(abs(p$upper - rf$upper), abs(p$lower - rf$lower))); np <- c(np, p$npts) }
    cat(sprintf(" %8.1e(%3d)", e, round(mean(np))))
    err_tab <- rbind(err_tab, data.frame(design = nm, method = m, max_err = e, mean_pts = mean(np))) }
  cat("\n") }
saveRDS(list(floor = floor_tab, errors = err_tab), "analysis/results-03b-stress.rds")
