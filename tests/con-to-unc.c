#include <stdio.h>
#include "nope.h"
#include <assert.h>

/**
 * min f(x) = 0.5*dot(x,x)
 * s.t x_1 - x_2 = 0
 *     x_1 = -1
 *
 */

#define UNUSED(x) (void)(x)

int nvar = 10;

Nope *nope;
#include "nope_interface.c"

void core_cfn (int *st, const int *n, const int *m, const double *x, double *f, double *c) {
  int i;
  UNUSED(st);
  UNUSED(m);
  *f = 0.0;
  for (i = 0; i < *n; i++)
    *f += x[i]*x[i];
  *f /= 2;
  c[0] = x[0] - x[1];
}

void core_cofg (int *st, const int *n, const double *x, double *f, double *g, bool *grad) {
  int i;
  UNUSED(st);
  *f = 0.0;
  for (i = 0; i < *n; i++)
    *f += x[i]*x[i];
  *f /= 2;
  if (*grad) {
    for (i = 0; i < *n; i++)
      g[i] = x[i];
  }
}

void core_ccfsg (int *st, const int *n, const int *m, const double *x, double *c, int *nnzj, const int
    *jmax, double *Jval, int *Jvar, int *Jfun, const bool *grad) {
  UNUSED(st);
  UNUSED(n);
  UNUSED(m);
  UNUSED(jmax);
  c[0] = x[0] - x[1];
  if (!*grad)
    return;
  *nnzj = 2;
  Jval[0] = 1;
  Jvar[0] = 1;
  Jfun[0] = 1;
  Jval[1] = -1;
  Jvar[1] = 2;
  Jfun[1] = 1;
}

void core_cdimen (int *st, const int *input, int *n, int *m) {
  UNUSED(st);
  UNUSED(input);
  *n = nvar;
  *m = 1;
}

void core_csetup (int *st, const int *input, const int *out, const int *io_buffer, int *n, int *m,
    double *x, double *bl, double *bu, double *y, double *cl, double *cu, bool
    *equatn, bool *linear, const int *e_order, const int *l_order, const int *v_order) {
  int i = 0;
  UNUSED(st);
  UNUSED(input);
  UNUSED(out);
  UNUSED(io_buffer);
  UNUSED(m);
  UNUSED(e_order);
  UNUSED(l_order);
  UNUSED(v_order);
  for (i = 0; i < *n; i++) {
    x[i] = 2;
    bl[i] = -1e20;
    bu[i] = 1e20;
  }
  bl[0] = -1;
  bu[0] = -1;
  y[0] = 0;
  cl[0] = 0;
  cu[0] = 0;
  equatn[0] = true;
  linear[0] = true;
}

void core_cdimsj (int *st, int *nnzj) {
  UNUSED(st);
  *nnzj = 2;
}

int main () {
  int n, m;

  nope = initializeNope();

  setFuncs(nope, core_cdimen, 0, 0, 0, 0, 0, core_csetup, 0, core_cfn,
      core_cofg, 0, core_ccfsg, core_cdimsj);
  runNope(nope);

  assert(nope->x[0] == -1.0);
  assert(nope->x[1] == -1.0);

  ppDIMEN(nope, &n, &m);
  double x[n], bl[n], bu[n];

  runUncSetup(nope, &n, x, bl, bu);
  for (int i = 0; i < n; i++)
    assert(x[i] == 2.0);

  // Testing if ufn works
  int st = 0;
  double f = 0;
  ufn(&st, &n, x, &f);
  assert(f == 1.0 + 2*n);
  // Testing if uofg works
  double g[n];
  bool grad = true;
  uofg(&st, &n, x, &f, g, &grad);
  for (int i = 0; i < n; i++)
    assert(g[i] == 2.0);

  destroyNope(nope);

  if (n != nvar-2 || m != 0) {
    printf("Error: nvar = %d, ncon = %d\n", n, m);
    return 1;
  }

  return 0;
}
