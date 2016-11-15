#include <stdio.h>
#include "nope.h"

/**
 * min f(x) = 0.5*(x_1^2 + x_2^2)
 * s.t x_1 + x_2 - 1 = 0
 */

#define UNUSED(x) (void)(x)

Nope *nope;
#include "nope_interface.h"

void core_cfn (int *st, const int *n, const int *m, const double *x, double *f, double *c) {
  UNUSED(st);
  UNUSED(n);
  UNUSED(m);
  *f = 0.5*(x[0]*x[0] + x[1]*x[1]);
  c[0] = x[0] + x[1] - 1;
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
  c[0] = x[0] + x[1] - 1;
  if (!*grad)
    return;
  *nnzj = 2;
  Jval[0] = 1;
  Jval[1] = 1;
  Jvar[0] = 1;
  Jvar[1] = 2;
  Jfun[0] = 1;
  Jfun[1] = 1;
}

void core_cdimen (int *st, const int *input, int *n, int *m) {
  UNUSED(st);
  UNUSED(input);
  *n = 2;
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
  UNUSED(n);
  UNUSED(m);
  UNUSED(e_order);
  UNUSED(l_order);
  UNUSED(v_order);
  for (i = 0; i < 2; i++) {
    x[i] = 0;
    bl[i] = -1e20;
    bu[i] = 1e20;
  }
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

  ppDIMEN(nope, &n, &m);

  destroyNope(nope);

  if (n != 2 || m != 1)
    return 1;

  return 0;
}
