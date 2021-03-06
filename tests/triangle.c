#include <stdio.h>
#include "nope.h"
#include <assert.h>

/**
 * min f(x) = 0.5*dot(x,x)
 * s.t x_1 = 1
 *     x_1 + x_2 = 2
 *     ...
 *     x_1 + x_2 + ... + x_{n-1} = n-1
 *     n variables
 *
 */

#define UNUSED(x) (void)(x)

int nvar = 10;

Nope *nope;
#include "nope_interface.h"

void core_cfn (int *st, const int *n, const int *m, const double *x, double *f, double *c) {
  int i, j;
  UNUSED(st);
  UNUSED(m);
  *f = 0.0;
  for (i = 0; i < *n; i++)
    *f = x[i]*x[i];
  *f /= 2;
  for (i = 0; i < *n-1; i++) {
    c[i] = -i-1;
    for (j = 0; j <= i; j++)
      c[i] += x[j];
  }
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
  int i, j;
  UNUSED(st);
  UNUSED(m);
  UNUSED(jmax);
  for (i = 0; i < *n-1; i++) {
    c[i] = -i-1;
    for (j = 0; j <= i; j++)
      c[i] += x[j];
  }
  if (!*grad)
    return;
  *nnzj = 0;
  for (i = 0; i < *n-1; i++) {
    for (j = 0; j <= i; j++) {
      Jval[*nnzj] = 1;
      Jvar[*nnzj] = j+1;
      Jfun[*nnzj] = i+1;
      (*nnzj)++;
    }
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
  *nnzj = nvar*(nvar-1)/2;
}

int main () {
  int n, m;

  nope = initializeNope();

  setFuncs(nope, core_cdimen, 0, 0, 0, 0, 0, core_csetup, 0, core_cfn,
      core_cofg, 0, core_ccfsg, core_cdimsj);
  runNope(nope);

  ppDIMEN(nope, &n, &m);

  destroyNope(nope);

  for (int i = 0; i < nvar-1; i++)
    assert(nope->x[i] == 1.0);

  if (n != 1 || m != 0) {
    printf("Error: nvar = %d, ncon = %d\n", n, m);
    return 1;
  }

  return 0;
}
