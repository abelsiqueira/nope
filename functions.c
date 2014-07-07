#include "preprocessor.h"

void ppDIMEN (Preprocessor *prep, int *nvar, int *ncon) {
  if (prep == 0 || prep->status != STATUS_OK)
    return;
  *nvar = prep->nvar - prep->nfix;
  *ncon = prep->ncon - prep->ntrivial;
}

void ppUFN (Preprocessor *prep, int * status, int * n, double * x, double * f) {
  if (prep == 0 || prep->status != STATUS_OK)
    return;
  (*prep->origin_ufn)(status, n, x, f);
}

void ppUOFG (Preprocessor *prep, int * status, int * n, double * x, double *
    f, double * g, _Bool * grad) {
  if (prep == 0 || prep->status != STATUS_OK)
    return;
  (*prep->origin_uofg)(status, n, x, f, g, grad);
}

void ppUHPROD (Preprocessor *prep, int * status, int * n, _Bool * goth,
    double * x, double * vector, double * result) {
  if (prep == 0 || prep->status != STATUS_OK)
    return;
  (*prep->origin_uhprod)(status, n, goth, x, vector, result);
}

void ppCFN (Preprocessor *prep, int * status, int * n, int * m, double * x,
    double * f, double * c) {
  int i = 0;
  if (prep == 0 || prep->status != STATUS_OK)
    return;
  for (i = 0; i < *n; i++)
    prep->x[prep->not_fixed_index[i]] = x[i];
  (*prep->origin_cfn)(status, &prep->nvar, m, prep->x, f, c);
}

void ppCOFG (Preprocessor *prep, int * status, int * n, double * x, double *
    f, double * g, _Bool * grad) {
  int i = 0;
  if (prep == 0 || prep->status != STATUS_OK)
    return;
  for (i = 0; i < *n; i++)
    prep->x[prep->not_fixed_index[i]] = x[i];
  (*prep->origin_cofg)(status, &prep->nvar, prep->x, f, prep->g, grad);
  for (i = 0; i < *n; i++)
    g[i] = prep->g[prep->not_fixed_index[i]];
}

void ppCHPROD (Preprocessor *prep, int * status, int * n, int * m, _Bool *
    goth, double * x, double * y, double * vector, double * result) {
  int i = 0;
  if (prep == 0 || prep->status != STATUS_OK)
    return;
  for (i = 0; i < prep->nvar; i++)
    prep->workspace1[i] = 0.0;
  for (i = 0; i < *n; i++) {
    prep->x[prep->not_fixed_index[i]] = x[i];
    prep->workspace1[prep->not_fixed_index[i]] = vector[i];
  }
  (*prep->origin_chprod)(status, &prep->nvar, m, goth, prep->x, y,
      prep->workspace1, prep->workspace2);
  for (i = 0; i < *n; i++)
    result[i] = prep->workspace2[prep->not_fixed_index[i]];
}

void ppCCFSG (Preprocessor *prep, int * status, int * n, int * m,
    double * x, double * c, int * nnzj, int * lj, double * Jval, int *
    Jvar, int * Jfun, _Bool * grad) {
  int i, k;
  if (prep == 0 || prep->status != STATUS_OK)
    return;
  for (i = 0; i < *n; i++)
    prep->x[prep->not_fixed_index[i]] = x[i];
  (*prep->origin_ccfsg)(status, &prep->nvar, m, prep->x, c, nnzj, lj, Jval,
      Jvar, Jfun, grad);
  printJacobian (prep->ncon, prep->nvar, *nnzj, Jval, Jvar, Jfun);
  // First the fixed variables are removed. This leaves some empty
  // columns
  for (k = 0; k < *nnzj; k++) {
    if (prep->is_fixed[Jvar[k]-1] || prep->is_trivial[Jfun[k]-1]) {
      Jval[k] = Jval[*nnzj-1];
      Jvar[k] = Jvar[*nnzj-1];
      Jfun[k] = Jfun[*nnzj-1];
      (*nnzj)--;
      k--;
    }
  }
  printJacobian (prep->ncon, prep->nvar, *nnzj, Jval, Jvar, Jfun);
  // Now we reorder the columns, reducing by one for each fixed
  // variable before that column
  for (i = prep->nfix-1; i >= 0; i--) {
    for (k = 0; k < *nnzj; k++) {
      if (Jvar[k]-1 > prep->fixed_index[i]) {
        Jvar[k]--;
      }
    }
  }

  for (i = prep->ntrivial-1; i >= 0; i--) {
    for (k = 0; k < *nnzj; k++) {
    }
  }
  printJacobian (prep->ncon, *n, *nnzj, Jval, Jvar, Jfun);
}
