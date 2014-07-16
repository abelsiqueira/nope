#include <stdio.h>
#include "nope.h"

/**
 * min f(x) = 0.5*dot(x,x)
 * s.t sum(x) - 1 = 0
 *     x_1 - x_2 = 0
 *     x_1 = -1
 */

#define UNUSED(x) (void)(x)

int nvar = 10;

Nope *nope;
#include "nope_interface.h"

void core_cfn (int *st, int *n, int *m, double *x, double *f, double *c) {
  int i;
  UNUSED(st);
  UNUSED(m);
  *f = 0.0;
  c[0] = -1;
  for (i = 0; i < *n; i++) {
    *f = x[i]*x[i];
    c[0] += x[i];
  }
  *f /= 2;
  c[1] = x[0] - x[1];
}

void core_cofg (int *st, int *n, double *x, double *f, double *g, bool *grad) {
  int i;
  UNUSED(st);
  *f = 0.0;
  for (i = 0; i < *n; i++)
    *f = x[i]*x[i];
  *f /= 2;
  if (*grad) {
    for (i = 0; i < *n; i++)
      g[i] = x[i];
  }
}

void core_chprod (int *st, int *n, int *m, bool *goth, double *x, double *y,
    double *p, double *q) {
  int i;
  UNUSED(st);
  UNUSED(m);
  UNUSED(goth);
  UNUSED(x);
  UNUSED(y);
  for (i = 0; i < *n; i++)
    q[i] = p[i];
}

void core_ccfsg (int *st, int *n, int *m, double *x, double *c, int *nnzj, int
    *jmax, double *Jval, int *Jvar, int *Jfun, bool *grad) {
  int i;
  UNUSED(st);
  UNUSED(m);
  UNUSED(jmax);
  c[0] = -1;
  for (i = 0; i < *n; i++)
    c[0] += x[i];
  c[1] = x[0] - x[1];
  if (!*grad)
    return;
  *nnzj = *n+2;
  for (i = 0; i < *n; i++) {
    Jval[i] = 1;
    Jvar[i] = i+1;
    Jfun[i] = 1;
  }
  Jval[*n] = 1;
  Jval[*n+1] = -1;
  Jvar[*n] = 1;
  Jvar[*n+1] = 2;
  Jfun[*n] = 2;
  Jfun[*n+1] = 2;
}

void core_cdimen (int *st, int *input, int *n, int *m) {
  UNUSED(st);
  UNUSED(input);
  *n = nvar;
  *m = 2;
}

void core_csetup (int *st, int *input, int *out, int *io_buffer, int *n, int *m,
    double *x, double *bl, double *bu, double *y, double *cl, double *cu, bool
    *equatn, bool *linear, int *e_order, int *l_order, int *v_order) {
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
    x[i] = 0;
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
  y[1] = 0;
  cl[1] = 0;
  cu[1] = 0;
  equatn[1] = true;
  linear[1] = true;
}

void core_cdimsj (int *st, int *nnzj) {
  UNUSED(st);
  *nnzj = nvar;
}

int main () {
  int n, m;

  nope = initializeNope();

  setFuncs(nope, core_cdimen, 0, 0, 0, 0, core_csetup, core_cfn, core_cofg,
      core_chprod, core_ccfsg, core_cdimsj);
  runNope(nope);

  ppDIMEN(nope, &n, &m);

  destroyNope(nope);

  if (n != nvar-2 || m != 1)
    return 1;

  return 0;
}
