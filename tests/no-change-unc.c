#include <stdio.h>
#include "nope.h"

/**
 * min f(x) = 0.5*(x_1^2 + x_2^2)
 * s.t x_1 + x_2 - 1 = 0
 */

#define UNUSED(x) (void)(x)

Nope *nope;
#include "nope_interface.h"

void core_ufn (int *st, const int *n, const double *x, double *f) {
  int i;
  UNUSED(st);
  *f = 0.0;
  for (i = 0; i < *n; i++)
    *f = x[i]*x[i];
  *f /= 2;
}

void core_uofg (int *st, const int *n, const double *x, double *f, double *g, const bool *grad) {
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

void core_cdimen (int *st, const int *input, int *n, int *m) {
  UNUSED(st);
  UNUSED(input);
  *n = 5;
  *m = 0;
}

void core_usetup (int *st, const int *input, const int *out, const int *io_buffer, int *n, double
    *x, double *bl, double *bu) {
  int i = 0;
  UNUSED(st);
  UNUSED(input);
  UNUSED(out);
  UNUSED(io_buffer);
  for (i = 0; i < *n; i++) {
    x[i] = 1;
    bl[i] = -1e20;
    bu[i] = 1e20;
  }
}

int main () {
  int n, m;

  nope = initializeNope();

  setFuncs(nope, core_cdimen, core_usetup, 0, core_ufn, core_uofg,
      0, 0, 0, 0, 0, 0, 0, 0);
  runNope(nope);

  ppDIMEN(nope, &n, &m);

  destroyNope(nope);

  if (n != 5 || m != 0)
    return 1;

  return 0;
}
