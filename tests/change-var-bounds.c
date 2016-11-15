#include <stdio.h>
#include "nope.h"
#include <assert.h>

/**
 * min f(x) = 0.5*(x_1^2 + x_2^2)
 * s.t     1 <= -2x_1 - 5x_2 <= 4.0
 *         5 <=  x_2  +  x_3 <= Inf
 *               x_1  +  x_3 = 0
 *           -1 <= x_1 <= 1
 *            0 <= x_2 <= 0
 *           -1 <= x_3 <= 7
 *
 */

#define UNUSED(x) (void)(x)

Nope *nope;
#include "nope_interface.h"

void core_cfn (int *st, const int *n, const int *m, const double *x, double *f, double *c) {
  UNUSED(st);
  UNUSED(n);
  UNUSED(m);
  *f = 0.5*(x[0]*x[0] + x[1]*x[1]);
  c[0] = -2*x[0] - 5*x[1];
  c[1] = x[1] + x[2];
  c[2] = x[0] + x[2];
}

void core_cofg (int *st, const int *n, const double *x, double *f, double *g, bool *grad) {
  UNUSED(st);
  UNUSED(n);
  *f = 0.5*(x[0]*x[0] + x[1]*x[1]);
  if (*grad) {
    g[0] = x[0];
    g[1] = x[1];
  }
}

void core_ccfsg (int *st, const int *n, const int *m, const double *x, double *c, int *nnzj, const int
    *jmax, double *Jval, int *Jvar, int *Jfun, const bool *grad) {
  UNUSED(st);
  UNUSED(n);
  UNUSED(m);
  UNUSED(jmax);
  c[0] = -2*x[0] - 5*x[1];
  c[1] = x[1] + x[2];
  c[2] = x[0] + x[2];
  if (!*grad)
    return;
  *nnzj = 6;
  Jval[0] = -2;
  Jval[1] = -5;
  Jval[2] = 1;
  Jval[3] = 1;
  Jval[4] = 1;
  Jval[5] = 1;
  Jvar[0] = 1;
  Jvar[1] = 2;
  Jvar[2] = 2;
  Jvar[3] = 3;
  Jvar[4] = 1;
  Jvar[5] = 3;
  Jfun[0] = 1;
  Jfun[1] = 1;
  Jfun[2] = 2;
  Jfun[3] = 2;
  Jfun[4] = 3;
  Jfun[5] = 3;
}

void core_cdimen (int *st, const int *input, int *n, int *m) {
  UNUSED(st);
  UNUSED(input);
  *n = 3;
  *m = 3;
}

void core_csetup (int *st, const int *input, const int *out, const int *io_buffer, int *n, int *m,
    double *x, double *bl, double *bu, double *y, double *cl, double *cu, bool
    *equatn, bool *linear, const int *e_order, const int *l_order, const int *v_order) {
  int i = 0;
  UNUSED(st);
  UNUSED(input);
  UNUSED(out);
  UNUSED(io_buffer);
  UNUSED(n);
  UNUSED(m);
  UNUSED(e_order);
  UNUSED(l_order);
  UNUSED(v_order);
  for (i = 0; i < 3; i++) {
    x[i] = 0;
    equatn[i] = false;
    linear[i] = true;
  }
  bl[0] = -1.0;
  bu[0] = 1.0;
  bl[1] = 0.0;
  bu[1] = 0.0;
  bl[2] = -1.0;
  bu[2] = 7.0;

  y[0] = 0;
  cl[0] =  1.0;
  cu[0] =  4.0;
  cl[1] = 5.0;
  cu[1] = 1e20;
  cl[2] = 0.0;
  cu[2] = 0.0;
  equatn[2] = true;
}

void core_cdimsj (int *st, int *nnzj) {
  UNUSED(st);
  *nnzj = 6;
}

int main () {
  int n, m;

  nope = initializeNope();

  setFuncs(nope, core_cdimen, 0, 0, 0, 0, 0, core_csetup, 0, core_cfn,
      core_cofg, 0, core_ccfsg, core_cdimsj);
  runNope(nope);

  ppDIMEN(nope, &n, &m);
  assert(n == 2);
  assert(m == 1);

  assert(nope->x[1] == 0.0);
  double x[n], bl[n], bu[n];
  double y[m], cl[m], cu[m];
  bool equatn[m], linear[m];
  int amax;
  runConSetup(nope, &n, x, bl, bu, &m, y, cl, cu, equatn, linear, &amax);
  assert(cl[0] == 0.0);
  assert(cu[0] == 0.0);
  assert(bl[0] == -1.0);
  assert(bu[0] == -0.5);
  assert(bl[1] == 5);
  assert(bu[1] == 7);

  destroyNope(nope);

  return 0;
}
