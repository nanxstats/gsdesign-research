# gsquad: formal checks for gsDesign's boundary search

Lean 4 + Mathlib proofs about the guarded Newton-Raphson update used in
gsDesign's `gsbound.c` / `gsbound1.c` after the `EXTREMEZ` change
(see `Gsquad/NewtonStep.lean`):

- the division branch of the update is only reached when the derivative is
  non-zero (no `0/0` step for either the lower or the upper bound);
- every branch moves the iterate by at most one unit;
- after `n` iterations the iterate is within `n` of its starting value, so the
  iteration cap `GS_MAXITER = 20` bounds the excursion of the search.

Build with `lake build` (Mathlib cache: `lake exe cache get`).
