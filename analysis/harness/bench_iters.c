#include <stdio.h>
#include <math.h>
#include "R.h"
#include "Rmath.h"
#include "gsDesign.h"
extern int gs_newton_iters, gs_newton_calls;
static double I3[3]  = {0.3333333333333333, 0.6666666666666666, 1.0};
static double a3[3]  = {-0.2387178, 0.9403609, 2.0031906};
static double I10[10] = {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};
int main(void) {
  int r = 18, k, ret, pe = 0; double tol = 1e-6, th = 0., b[10], plo[10], a[10];
  double fp3[3] = {0.001, 0.0056, 0.0184}, fp10[10] = {0.00009,0.00028,0.00061,0.00110,0.00178,0.00266,0.00378,0.00516,0.00683,0.00871};
  double am20[10] = {-20,-20,-20,-20,-20,-20,-20,-20,-20,-20};
  k = 3; gs_newton_iters = gs_newton_calls = 0; gsbound1(&k, &th, I3, am20, b, plo, fp3, &tol, &r, &ret, &pe);
  printf("gsbound1 k=3 one-sided:   %d Newton iterations over %d analyses (b: %.7f %.7f %.7f)\n", gs_newton_iters, gs_newton_calls, b[0], b[1], b[2]);
  for (int i = 0; i < 3; i++) a[i] = a3[i];
  gs_newton_iters = gs_newton_calls = 0; gsbound1(&k, &th, I3, a, b, plo, fp3, &tol, &r, &ret, &pe);
  printf("gsbound1 k=3 lower bound: %d Newton iterations over %d analyses\n", gs_newton_iters, gs_newton_calls);
  k = 10; gs_newton_iters = gs_newton_calls = 0; gsbound1(&k, &th, I10, am20, b, plo, fp10, &tol, &r, &ret, &pe);
  printf("gsbound1 k=10 one-sided:  %d Newton iterations over %d analyses\n", gs_newton_iters, gs_newton_calls);
  double thd = 3.24; /* theta such that mu is moderate; beta-spending style call with theta=-delta */
  double bs10[10] = {0.005,0.008,0.010,0.011,0.012,0.012,0.012,0.011,0.010,0.009}; double neg = -thd; double nb[10]; for (int i = 0; i < 10; i++) nb[i] = -3.0;
  gs_newton_iters = gs_newton_calls = 0; gsbound1(&k, &neg, I10, nb, b, plo, bs10, &tol, &r, &ret, &pe);
  printf("gsbound1 k=10 theta=-3.24 (beta-spending style): %d Newton iterations over %d analyses, err=%d\n", gs_newton_iters, gs_newton_calls, ret);
  gs_newton_iters = gs_newton_calls = 0; gsbound(&k, I10, a, b, fp10, fp10, &tol, &r, &ret, &pe);
  printf("gsbound k=10 symmetric:   %d Newton iterations over %d analyses\n", gs_newton_iters, gs_newton_calls);
  /* zero spend at an interim -> how many iterations wasted? */
  double fpz[3] = {0.001, 0.0, 0.024};
  gs_newton_iters = gs_newton_calls = 0; k = 3; gsbound1(&k, &th, I3, am20, b, plo, fpz, &tol, &r, &ret, &pe);
  printf("gsbound1 k=3 zero spend at analysis 2: %d Newton iterations over %d analyses; b=%.4f %.4f %.4f\n", gs_newton_iters, gs_newton_calls, b[0], b[1], b[2]);

  /* erfc vs pnorm: absolute and relative (p>1e-30) differences */
  { double maxabs = 0, maxrel = 0;
    for (double xx = -12; xx < 12; xx += 0.0005) { double p1 = pnorm(xx, 0., 1., 1, 0), p2 = 0.5 * erfc(-xx * M_SQRT1_2);
      double d = fabs(p1 - p2); if (d > maxabs) maxabs = d; if (p1 > 1e-30 && d / p1 > maxrel) maxrel = d / p1; }
    printf("pnorm vs 0.5*erfc: max abs diff %.3e, max rel diff (p>1e-30) %.3e\n", maxabs, maxrel);
  }
  return 0;
}
