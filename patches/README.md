# Patch series against gsDesign 3.11.0

Generated with `git format-patch` from a working copy of
[`keaven/gsDesign@v3.11.0`](https://github.com/keaven/gsDesign/tree/v3.11.0).
Apply in order from the root of a gsDesign checkout at v3.11.0:

```sh
git am patches/0001-*.patch patches/0002-*.patch patches/0003-*.patch patches/0004-*.patch
```

Each patch leaves the package building and its full test suite passing
(`devtools::test()`), so they can also be submitted as separate pull
requests in this order.

| # | Scope | Behavior change |
|---|-------|-----------------|
| 0001 | Tests only: pins the 3.11.0 outputs of the integration layer (fixtures at 13 significant digits), adds exact-reference tests against `mvtnorm::pmvnorm(Miwa)`, and snapshot tests for analyses without a bound. | None |
| 0002 | `EXTREMEZ` sentinel: absent bounds are `-Inf`/`Inf` end to end (C and R), Newton iteration cap and clamp become named constants, `0/0` guard, `NAOK = TRUE` in `.C()` calls, NaN-safe convergence checks, `gsBoundCP()`/`gsBoundSummary()` report `NA` at absent bounds. | Breaking for code comparing bounds with 20; numerically negligible (see report) |
| 0003 | Performance, same algorithm: restructured density update, lower-tail probability hoisted out of the Newton loop, constants hoisted, prototypes in `gsDesign.h`, optional `-DGS_USE_ERFC`. | None (agreement with 3.11.0 to about 1e-16) |
| 0004 | Gauss-Legendre quadrature on the continuation region, opt-in via `options(gsDesign.quadrature = "gl")`. | None by default; opt-in changes results within the accuracy of the default grid (about 1e-7) while being about 100x more accurate and several times faster |

See `docs/index.html` (the research report) for the evidence behind each patch.
