#include <stdio.h>
#include <math.h>
#include <time.h>
#include "R.h"
#include "Rmath.h"
#include "gsDesign.h"
static double now(void) { struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts); return ts.tv_sec + ts.tv_nsec * 1e-9; }
void base_probrej(int *, int *, double *, double *, double *, double *, double *, double *, int *);
void base_gsbound1(int *, double *, double *, double *, double *, double *, double *, double *, int *, int *, int *);
void base_gsbound(int *, double *, double *, double *, double *, double *, double *, int *, int *, int *);
static double I10[10] = {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};
static double a10[10] = {-0.9718, -0.3560, 0.0938, 0.4562, 0.7596, 1.0233, 1.2604, 1.4802, 1.6911, 2.0177};
static double b10[10] = {3.7495, 3.4103, 3.1719, 2.9905, 2.8434, 2.7189, 2.6101, 2.5127, 2.4238, 2.3418};
static double fp10[10] = {0.00009,0.00028,0.00061,0.00110,0.00178,0.00266,0.00378,0.00516,0.00683,0.00871};
static double maxdiff(double *x, double *y, int n) { double d = 0; for (int i = 0; i < n; i++) if (fabs(x[i] - y[i]) > d) d = fabs(x[i] - y[i]); return d; }
#define TIME(label, reps, stmt) do { double t0 = now(); for (int _i = 0; _i < (reps); _i++) { stmt; } double t1 = now(); printf("%-40s %9.2f us\n", label, (t1 - t0) / (reps) * 1e6); } while (0)
int main(void) { int r = 18, k = 10, nt = 1, ret, pe = 0; double th = 0.2, tol = 1e-6, plo[10], phi[10], plo2[10], phi2[10], b[10], b2[10], a[10], a2[10], pl[10], pl2[10], m20[10]; for (int i = 0; i < 10; i++) m20[i] = -20;
  base_probrej(&k, &nt, &th, I10, a10, b10, plo, phi, &r); probrej(&k, &nt, &th, I10, a10, b10, plo2, phi2, &r);
  printf("probrej  patch3 vs base: %.2e %.2e\n", maxdiff(phi, phi2, 10), maxdiff(plo, plo2, 10));
  th = 0; base_gsbound1(&k, &th, I10, m20, b, pl, fp10, &tol, &r, &ret, &pe); gsbound1(&k, &th, I10, m20, b2, pl2, fp10, &tol, &r, &ret, &pe);
  printf("gsbound1 patch3 vs base: %.2e %.2e\n", maxdiff(b, b2, 10), maxdiff(pl, pl2, 10));
  base_gsbound(&k, I10, a, b, fp10, fp10, &tol, &r, &ret, &pe); gsbound(&k, I10, a2, b2, fp10, fp10, &tol, &r, &ret, &pe);
  printf("gsbound  patch3 vs base: %.2e %.2e\n", maxdiff(a, a2, 10), maxdiff(b, b2, 10));
  th = 0.2; TIME("probrej k=10 two-sided base", 1000, base_probrej(&k, &nt, &th, I10, a10, b10, plo, phi, &r)); TIME("probrej k=10 two-sided patch3", 1000, probrej(&k, &nt, &th, I10, a10, b10, plo, phi, &r));
  TIME("probrej k=10 one-sided base", 300, base_probrej(&k, &nt, &th, I10, m20, b10, plo, phi, &r)); TIME("probrej k=10 one-sided patch3", 300, probrej(&k, &nt, &th, I10, m20, b10, plo, phi, &r));
  th = 0; TIME("gsbound1 k=10 one-sided base", 300, base_gsbound1(&k, &th, I10, m20, b, pl, fp10, &tol, &r, &ret, &pe)); TIME("gsbound1 k=10 one-sided patch3", 300, gsbound1(&k, &th, I10, m20, b, pl, fp10, &tol, &r, &ret, &pe));
  TIME("gsbound k=10 symmetric base", 300, base_gsbound(&k, I10, a, b, fp10, fp10, &tol, &r, &ret, &pe)); TIME("gsbound k=10 symmetric patch3", 300, gsbound(&k, I10, a, b, fp10, fp10, &tol, &r, &ret, &pe));
  return 0; }
