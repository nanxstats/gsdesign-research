/* Standalone harness: times the pristine gsDesign 3.11.0 C routines. Links against libR for pnorm/qnorm/Rprintf. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "R.h"
#include "Rmath.h"
#include "gsDesign.h"

static double now(void) { struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts); return ts.tv_sec + ts.tv_nsec * 1e-9; }
int gridpts(int, double, double, double, double *, double *);
void h1(double, int, double *, double, double *, double *);
void hupdate(double, double *, int, double, double *, double *, int, double, double *, double *);

/* designs exported from R: gsDesign(k=3) tt=4 default and k=10 tt=4 */
static double I3[3]  = {0.3333333333333333, 0.6666666666666666, 1.0};
static double a3[3]  = {-0.2387178, 0.9403609, 2.0031906};
static double b3[3]  = {3.0106843, 2.5464966, 1.9992261};
static double I10[10] = {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};
static double a10[10] = {-0.9718, -0.3560, 0.0938, 0.4562, 0.7596, 1.0233, 1.2604, 1.4802, 1.6911, 2.0177};
static double b10[10] = {3.7495, 3.4103, 3.1719, 2.9905, 2.8434, 2.7189, 2.6101, 2.5127, 2.4238, 2.3418};

int main(void) {
  int r = 18, k, n, reps; double t0, t1;
  /* primitives */
  { volatile double acc = 0; double x = 0.1; n = 2000000;
    t0 = now(); for (int i = 0; i < n; i++) { acc += exp(-x * x / 2); x += 1e-6; } t1 = now();
    printf("exp():        %6.2f ns/call\n", (t1 - t0) / n * 1e9);
    x = -3; t0 = now(); for (int i = 0; i < n; i++) { acc += pnorm(x, 0., 1., 1, 0); x += 3e-6; } t1 = now();
    printf("Rf_pnorm5():  %6.2f ns/call (x in [-3,3])\n", (t1 - t0) / n * 1e9);
    x = -3; t0 = now(); for (int i = 0; i < n; i++) { acc += 0.5 * erfc(-x * M_SQRT1_2); x += 3e-6; } t1 = now();
    printf("erfc()-based: %6.2f ns/call (x in [-3,3])\n", (t1 - t0) / n * 1e9);
    x = -10; t0 = now(); for (int i = 0; i < n; i++) { acc += pnorm(x, 0., 1., 1, 0); x += 1e-5; } t1 = now();
    printf("Rf_pnorm5():  %6.2f ns/call (x in [-10,10])\n", (t1 - t0) / n * 1e9);
    x = -10; t0 = now(); for (int i = 0; i < n; i++) { acc += 0.5 * erfc(-x * M_SQRT1_2); x += 1e-5; } t1 = now();
    printf("erfc()-based: %6.2f ns/call (x in [-10,10])\n", (t1 - t0) / n * 1e9);
    /* accuracy erfc vs pnorm */
    double maxrel = 0; for (double xx = -38; xx < 8; xx += 0.001) { double p1 = pnorm(xx, 0., 1., 1, 0), p2 = 0.5 * erfc(-xx * M_SQRT1_2); if (p1 > 0) { double rel = fabs(p1 - p2) / p1; if (rel > maxrel) maxrel = rel; } }
    printf("max rel diff pnorm vs erfc on [-38,8]: %.3e\n", maxrel);
  }
  /* grid sizes */
  { double z[1000], w[1000];
    printf("grid size m+1: full=%d  [a3,b3] at analysis1 theta=0: %d  one-sided [-20,b]: %d  [-Inf,b]: %d\n",
      gridpts(r, 0., -20, 20, z, w) + 1, gridpts(r, 0., a3[0], b3[0], z, w) + 1, gridpts(r, 0., -20, b3[0], z, w) + 1, gridpts(r, 0., -INFINITY, b3[0], z, w) + 1);
  }
  /* hupdate alone: full grids (one-sided-like) and continuation-region grids */
  { double z1[1000], w1[1000], h[1000], z2[1000], w2[1000], h2[1000]; int m1, m2;
    m1 = gridpts(r, 0., -20, b3[0], z1, w1); h1(0., m1, w1, I3[0], z1, h); m2 = gridpts(r, 0., -20, b3[1], z2, w2);
    reps = 2000; t0 = now(); for (int i = 0; i < reps; i++) hupdate(0., w2, m1, I3[0], z1, h, m2, I3[1], z2, h2); t1 = now();
    printf("hupdate one-sided grids m1=%d m2=%d: %8.2f us  (%.2f ns per pair)\n", m1 + 1, m2 + 1, (t1 - t0) / reps * 1e6, (t1 - t0) / reps / ((m1 + 1.0) * (m2 + 1.0)) * 1e9);
    m1 = gridpts(r, 0., a3[0], b3[0], z1, w1); h1(0., m1, w1, I3[0], z1, h); m2 = gridpts(r, 0., a3[1], b3[1], z2, w2);
    t0 = now(); for (int i = 0; i < reps; i++) hupdate(0., w2, m1, I3[0], z1, h, m2, I3[1], z2, h2); t1 = now();
    printf("hupdate two-sided grids m1=%d m2=%d: %8.2f us  (%.2f ns per pair)\n", m1 + 1, m2 + 1, (t1 - t0) / reps * 1e6, (t1 - t0) / reps / ((m1 + 1.0) * (m2 + 1.0)) * 1e9);
    m1 = gridpts(r, 0., -20, 20, z1, w1); h1(0., m1, w1, I3[0], z1, h); m2 = gridpts(r, 0., -20, 20, z2, w2);
    t0 = now(); for (int i = 0; i < reps; i++) hupdate(0., w2, m1, I3[0], z1, h, m2, I3[1], z2, h2); t1 = now();
    printf("hupdate full grids m1=%d m2=%d:      %8.2f us  (%.2f ns per pair)\n", m1 + 1, m2 + 1, (t1 - t0) / reps * 1e6, (t1 - t0) / reps / ((m1 + 1.0) * (m2 + 1.0)) * 1e9);
  }
  /* probrej */
  { double plo[10], phi[10], th = 0.2; int nt = 1;
    k = 3; reps = 5000; t0 = now(); for (int i = 0; i < reps; i++) probrej(&k, &nt, &th, I3, a3, b3, plo, phi, &r); t1 = now();
    printf("probrej k=3 two-sided:  %8.2f us\n", (t1 - t0) / reps * 1e6);
    double am20[10] = {-20,-20,-20,-20,-20,-20,-20,-20,-20,-20};
    t0 = now(); for (int i = 0; i < reps; i++) probrej(&k, &nt, &th, I3, am20, b3, plo, phi, &r); t1 = now();
    printf("probrej k=3 one-sided:  %8.2f us\n", (t1 - t0) / reps * 1e6);
    k = 10; reps = 1000; t0 = now(); for (int i = 0; i < reps; i++) probrej(&k, &nt, &th, I10, a10, b10, plo, phi, &r); t1 = now();
    printf("probrej k=10 two-sided: %8.2f us\n", (t1 - t0) / reps * 1e6);
    t0 = now(); for (int i = 0; i < reps; i++) probrej(&k, &nt, &th, I10, am20, b10, plo, phi, &r); t1 = now();
    printf("probrej k=10 one-sided: %8.2f us\n", (t1 - t0) / reps * 1e6);
  }
  /* gsbound1 / gsbound */
  { double fp3[3] = {0.001, 0.0056, 0.0184}, fp10[10] = {0.00009,0.00028,0.00061,0.00110,0.00178,0.00266,0.00378,0.00516,0.00683,0.00871};
    double a[10], b[10], plo[10], tol = 1e-6, th = 0.; int ret, pe = 0;
    double am20[10] = {-20,-20,-20,-20,-20,-20,-20,-20,-20,-20};
    k = 3; reps = 2000; t0 = now(); for (int i = 0; i < reps; i++) gsbound1(&k, &th, I3, am20, b, plo, fp3, &tol, &r, &ret, &pe); t1 = now();
    printf("gsbound1 k=3 one-sided (a=-20):  %8.2f us  b=%.6f %.6f %.6f\n", (t1 - t0) / reps * 1e6, b[0], b[1], b[2]);
    for (int i = 0; i < 3; i++) a[i] = a3[i];
    t0 = now(); for (int i = 0; i < reps; i++) gsbound1(&k, &th, I3, a, b, plo, fp3, &tol, &r, &ret, &pe); t1 = now();
    printf("gsbound1 k=3 with lower bound:   %8.2f us  b=%.6f %.6f %.6f\n", (t1 - t0) / reps * 1e6, b[0], b[1], b[2]);
    k = 10; reps = 300; t0 = now(); for (int i = 0; i < reps; i++) gsbound1(&k, &th, I10, am20, b, plo, fp10, &tol, &r, &ret, &pe); t1 = now();
    printf("gsbound1 k=10 one-sided (a=-20): %8.2f us\n", (t1 - t0) / reps * 1e6);
    t0 = now(); for (int i = 0; i < reps; i++) gsbound(&k, I10, a, b, fp10, fp10, &tol, &r, &ret, &pe); t1 = now();
    printf("gsbound  k=10 symmetric:         %8.2f us\n", (t1 - t0) / reps * 1e6);
  }
  return 0;
}
