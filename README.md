# gsDesign numerical integration research

Research on the core numerical integration algorithm of the
[gsDesign](https://github.com/keaven/gsDesign) R package (v3.11.0):

1. replacing the finite `EXTREMEZ = 20` sentinel for absent bounds by `Inf`;
2. a review of rpact's integration routine against Jennison & Turnbull (2000, Chapter 19);
3. faster implementations of the recursion and a Gauss-Legendre quadrature that is
   both faster and more accurate.

- **Report**: [`docs/index.html`](docs/index.html) (rendered from `report/index.qmd`).
- **Patches**: `patches/` (a four-patch series against gsDesign v3.11.0, see `patches/README.md`).
- **Analysis scripts and results**: `analysis/` (R scripts `00`-`07`, a standalone C harness under
  `analysis/harness/`, and `results-*` files consumed by the report).
- **Formal checks**: `lean/gsquad/` (Lean 4 + Mathlib proofs about the guarded Newton update).

`deps-src/` (vendored sources, managed by `okr`) and `work/` (the git working copy of gsDesign
in which the patches were developed) are not tracked.
