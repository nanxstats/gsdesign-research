/* Harness: pristine (baseline) vs v2 routines (JT & GL quadrature). Baseline symbols are renamed via macros. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "R.h"
#include "Rmath.h"
#include "gsDesign.h"   /* v2 header */
static double now(void) { struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts); return ts.tv_sec + ts.tv_nsec * 1e-9; }
/* baseline entry points (compiled from the pristine sources with renamed symbols) */
void base_probrej(int *, int *, double *, double *, double *, double *, double *, double *, int *);
void base_gsbound1(int *, double *, double *, double *, double *, double *, double *, double *, int *, int *, int *);
void base_gsbound(int *, double *, double *, double *, double *, double *, double *, int *, int *, int *);

static double I3[3]  = {0.3333333333333333, 0.6666666666666666, 1.0};
static double a3[3]  = {-0.2387178, 0.9403609, 2.0031906};
static double b3[3]  = {3.0106843, 2.5464966, 1.9992261};
static double I10[10] = {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};
static double a10[10] = {-0.9718, -0.3560, 0.0938, 0.4562, 0.7596, 1.0233, 1.2604, 1.4802, 1.6911, 2.0177};
static double b10[10] = {3.7495, 3.4103, 3.1719, 2.9905, 2.8434, 2.7189, 2.6101, 2.5127, 2.4238, 2.3418};
static double fp3[3] = {0.001, 0.0056, 0.0184}, fp10[10] = {0.00009,0.00028,0.00061,0.00110,0.00178,0.00266,0.00378,0.00516,0.00683,0.00871};
static double ninf[10], m20[10];

#define TIME(label, reps, stmt) do { double t0 = now(); for (int _i = 0; _i < (reps); _i++) { stmt; } double t1 = now(); printf("%-46s %9.2f us\n", label, (t1 - t0) / (reps) * 1e6); } while (0)
static double maxdiff(double *x, double *y, int n) { double d = 0; for (int i = 0; i < n; i++) if (fabs(x[i] - y[i]) > d) d = fabs(x[i] - y[i]); return d; }

int main(void) {
  int r = 18, k, nt = 1, jt = GS_QUAD_JT, gl = GS_QUAD_GL, ret, pe = 0; double tol = 1e-6, th, plo[10], phi[10], plo2[10], phi2[10], b[10], b2[10], pl[10], pl2[10], a[10];
  for (int i = 0; i < 10; i++) { ninf[i] = -INFINITY; m20[i] = -20; }
  /* correctness: v2-JT vs baseline */
  th = 0.2; k = 10; base_probrej(&k, &nt, &th, I10, a10, b10, plo, phi, &r); probrej(&k, &nt, &th, I10, a10, b10, plo2, phi2, &r, &jt);
  printf("probrej k=10 two-sided  v2JT vs base: max|dp| = %.2e (upper) %.2e (lower)\n", maxdiff(phi, phi2, 10), maxdiff(plo, plo2, 10));
  base_probrej(&k, &nt, &th, I10, m20, b10, plo, phi, &r); probrej(&k, &nt, &th, I10, ninf, b10, plo2, phi2, &r, &jt);
  printf("probrej k=10 one-sided  v2JT vs base: max|dp| = %.2e (upper) %.2e (lower)\n", maxdiff(phi, phi2, 10), maxdiff(plo, plo2, 10));
  probrej(&k, &nt, &th, I10, a10, b10, plo2, phi2, &r, &gl); base_probrej(&k, &nt, &th, I10, a10, b10, plo, phi, &r);
  printf("probrej k=10 two-sided  GL   vs base: max|dp| = %.2e (upper) %.2e (lower)  [JT r=18 error ~1e-7]\n", maxdiff(phi, phi2, 10), maxdiff(plo, plo2, 10));
  { int r80 = 80; base_probrej(&k, &nt, &th, I10, a10, b10, plo, phi, &r80); }
  printf("probrej k=10 two-sided  GL   vs base r=80: max|dp| = %.2e (upper) %.2e (lower)\n", maxdiff(phi, phi2, 10), maxdiff(plo, plo2, 10));
  th = 0.; base_gsbound1(&k, &th, I10, m20, b, pl, fp10, &tol, &r, &ret, &pe); gsbound1(&k, &th, I10, ninf, b2, pl2, fp10, &tol, &r, &ret, &pe, &jt);
  printf("gsbound1 k=10 one-sided v2JT vs base: max|db| = %.2e  max|dplo| = %.2e\n", maxdiff(b, b2, 10), maxdiff(pl, pl2, 10));
  gsbound1(&k, &th, I10, ninf, b2, pl2, fp10, &tol, &r, &ret, &pe, &gl);
  printf("gsbound1 k=10 one-sided GL   vs base: max|db| = %.2e  (b_GL[9]=%.8f b_base[9]=%.8f)\n", maxdiff(b, b2, 10), b2[9], b[9]);
  base_gsbound(&k, I10, a, b, fp10, fp10, &tol, &r, &ret, &pe); { double a2[10]; gsbound(&k, I10, a2, b2, fp10, fp10, &tol, &r, &ret, &pe, &jt);
  printf("gsbound  k=10 symmetric v2JT vs base: max|da| = %.2e max|db| = %.2e\n", maxdiff(a, a2, 10), maxdiff(b, b2, 10)); }
  printf("\n--- timings (median-ish over repeats) ---\n");
  th = 0.2;
  k = 3;  TIME("probrej k=3 two-sided   base", 5000, base_probrej(&k, &nt, &th, I3, a3, b3, plo, phi, &r));
          TIME("probrej k=3 two-sided   v2 JT", 5000, probrej(&k, &nt, &th, I3, a3, b3, plo, phi, &r, &jt));
          TIME("probrej k=3 two-sided   v2 GL", 5000, probrej(&k, &nt, &th, I3, a3, b3, plo, phi, &r, &gl));
          TIME("probrej k=3 one-sided   base", 2000, base_probrej(&k, &nt, &th, I3, m20, b3, plo, phi, &r));
          TIME("probrej k=3 one-sided   v2 JT", 2000, probrej(&k, &nt, &th, I3, ninf, b3, plo, phi, &r, &jt));
          TIME("probrej k=3 one-sided   v2 GL", 2000, probrej(&k, &nt, &th, I3, ninf, b3, plo, phi, &r, &gl));
  k = 10; TIME("probrej k=10 two-sided  base", 1000, base_probrej(&k, &nt, &th, I10, a10, b10, plo, phi, &r));
          TIME("probrej k=10 two-sided  v2 JT", 1000, probrej(&k, &nt, &th, I10, a10, b10, plo, phi, &r, &jt));
          TIME("probrej k=10 two-sided  v2 GL", 1000, probrej(&k, &nt, &th, I10, a10, b10, plo, phi, &r, &gl));
          TIME("probrej k=10 one-sided  base", 300, base_probrej(&k, &nt, &th, I10, m20, b10, plo, phi, &r));
          TIME("probrej k=10 one-sided  v2 JT", 300, probrej(&k, &nt, &th, I10, ninf, b10, plo, phi, &r, &jt));
          TIME("probrej k=10 one-sided  v2 GL", 300, probrej(&k, &nt, &th, I10, ninf, b10, plo, phi, &r, &gl));
  th = 0.;
  k = 3;  TIME("gsbound1 k=3 one-sided  base", 2000, base_gsbound1(&k, &th, I3, m20, b, pl, fp3, &tol, &r, &ret, &pe));
          TIME("gsbound1 k=3 one-sided  v2 JT", 2000, gsbound1(&k, &th, I3, ninf, b, pl, fp3, &tol, &r, &ret, &pe, &jt));
          TIME("gsbound1 k=3 one-sided  v2 GL", 2000, gsbound1(&k, &th, I3, ninf, b, pl, fp3, &tol, &r, &ret, &pe, &gl));
          TIME("gsbound1 k=3 with lower base", 5000, base_gsbound1(&k, &th, I3, a3, b, pl, fp3, &tol, &r, &ret, &pe));
          TIME("gsbound1 k=3 with lower v2 JT", 5000, gsbound1(&k, &th, I3, a3, b, pl, fp3, &tol, &r, &ret, &pe, &jt));
          TIME("gsbound1 k=3 with lower v2 GL", 5000, gsbound1(&k, &th, I3, a3, b, pl, fp3, &tol, &r, &ret, &pe, &gl));
  k = 10; TIME("gsbound1 k=10 one-sided base", 300, base_gsbound1(&k, &th, I10, m20, b, pl, fp10, &tol, &r, &ret, &pe));
          TIME("gsbound1 k=10 one-sided v2 JT", 300, gsbound1(&k, &th, I10, ninf, b, pl, fp10, &tol, &r, &ret, &pe, &jt));
          TIME("gsbound1 k=10 one-sided v2 GL", 300, gsbound1(&k, &th, I10, ninf, b, pl, fp10, &tol, &r, &ret, &pe, &gl));
          TIME("gsbound  k=10 symmetric base", 300, base_gsbound(&k, I10, a, b, fp10, fp10, &tol, &r, &ret, &pe));
          TIME("gsbound  k=10 symmetric v2 JT", 300, gsbound(&k, I10, a, b, fp10, fp10, &tol, &r, &ret, &pe, &jt));
          TIME("gsbound  k=10 symmetric v2 GL", 300, gsbound(&k, I10, a, b, fp10, fp10, &tol, &r, &ret, &pe, &gl));
  return 0;
}
