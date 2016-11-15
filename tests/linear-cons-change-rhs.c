#include <stdio.h>
#include "nope.h"
#include <assert.h>

/**
 * min f(x) = 0.5*dot(x,x)
 * s.t x_i - x_{i+1} + 1 = 0, i = 1,2,...,n-1
 *     1 <= x_1 <= 1
 *
 * Notice that this is not Ax = b, but Ax - b = 0, so the nope->rhs should be
 * corrected.
 */

#define UNUSED(x) (void)(x)

int nvar = 10;

Nope *nope;
#include "nope_interface.h"

void core_cfn (int *st, const int *n, const int *m, const double *x, double *f, double *c) {
  int i;
  UNUSED(st);
  UNUSED(m);
  *f = 0.0;
  for (i = 0; i < *n; i++)
    *f = x[i]*x[i];
  *f /= 2;
  for (i = 0; i < *m; i++)
    c[i] = x[i] - x[i+1] + 1;
}

void core_cofg (int *st, const int *n, const double *x, double *f, double *g, bool *grad) {
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

void core_ccfsg (int *st, const int *n, const int *m, const double *x, double *c, int *nnzj, const int
    *jmax, double *Jval, int *Jvar, int *Jfun, const bool *grad) {
  int i;
  UNUSED(st);
  UNUSED(n);
  UNUSED(m);
  UNUSED(jmax);
  for (i = 0; i < *m; i++)
    c[i] = x[i] - x[i+1] + 1;
  if (!*grad)
    return;
  *nnzj = 2*(*m);
  for (i = 0; i < *m; i++) {
    Jval[2*i] = 1;
    Jval[2*i+1] = -1;
    Jvar[2*i] = i+1;
    Jvar[2*i+1] = i+2;
    Jfun[2*i] = i+1;
    Jfun[2*i+1] = i+1;
  }
}

void core_cdimen (int *st, const int *input, int *n, int *m) {
  UNUSED(st);
  UNUSED(input);
  *n = nvar;
  *m = nvar-1;
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
    x[i] = 0;
    bl[i] = -1e20;
    bu[i] = 1e20;
  }
  bl[0] = 1;
  bu[0] = 1;
  for (i = 0; i < *m; i++) {
    y[i] = 0;
    cl[i] = 0;
    cu[i] = 0;
    equatn[i] = true;
    linear[i] = true;
  }
}

void core_cdimsj (int *st, int *nnzj) {
  UNUSED(st);
  *nnzj = 2*(nvar-1);
}

int main () {
  int n, m;

  nope = initializeNope();

  setFuncs(nope, core_cdimen, 0, 0, 0, 0, 0, core_csetup, 0, core_cfn,
      core_cofg, 0, core_ccfsg, core_cdimsj);
  runNope(nope);

  ppDIMEN(nope, &n, &m);

  // Test the solution if [1, 2, 3, ..., n]
  for (int i = 0; i < nope->nvar; i++) {
    assert(nope->x[i] == i+1);
  }
  // int st = 0;
  // double f = 1e20;
  // nope->origin_cfn(&st, &nope->nvar, &nope->ncon, nope->x, &f, nope->c);

  destroyNope(nope);

  if (n != 0 || m != 0)
    return 1;

  return 0;
}
