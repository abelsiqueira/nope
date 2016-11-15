/**
 * NOPE
 * Copyright (C) 2014  Abel Soares Siqueira
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#include "nope.h"
#include <assert.h>

#define ASSERT_NOPE_OK assert( (nope != 0) && (nope->status == STATUS_OK) )

void ppDIMEN (Nope *nope, int *nvar, int *ncon) {
  ASSERT_NOPE_OK;
  *nvar = nope->nvar - nope->nfix;
  *ncon = nope->ncon - nope->ntrivial;
}

void ppUNAMES (Nope *nope, int * status, const int * n, char * pname, char * vnames) {
  char v[10*nope->nvar];
  UNUSED(n);
  UNUSED(vnames);
  // We, knowingly, decided to not pass the names of the variables
  if (nope->ncon == 0)
    (*nope->origin_unames)(status, &nope->nvar, pname, v);
  else {
    char c[10*nope->ncon];
    (*nope->origin_cnames)(status, &nope->nvar, &nope->ncon, pname, v, c);
  }
}

void ppUFN (Nope *nope, int * status, const int * n, const double * x, double * f) {
  int i = 0;
  ASSERT_NOPE_OK;
  for (i = 0; i < *n; i++) {
    nope->x[nope->not_fixed_index[i]] = x[i];
  }
  if (nope->ncon == 0)
    (*nope->origin_ufn)(status, &nope->nvar, nope->x, f);
  else {
    (*nope->origin_cfn)(status, &nope->nvar, &nope->ncon, nope->x, f, nope->c);
  }
}

void ppUOFG (Nope *nope, int * status, const int * n, const double * x, double *
    f, double * g, const _Bool * grad) {
  int i = 0;
  ASSERT_NOPE_OK;
  for (i = 0; i < *n; i++)
    nope->x[nope->not_fixed_index[i]] = x[i];
  if (nope->ncon == 0)
    (*nope->origin_uofg)(status, &nope->nvar, nope->x, f, nope->g, grad);
  else {
    // CUTEst uses cofg with non-const grad, but uofg with const grad
    bool grad_ = *grad;
    (*nope->origin_cofg)(status, &nope->nvar, nope->x, f, nope->g, &grad_);
  }
  if (*grad) {
    for (i = 0; i < *n; i++)
      g[i] = nope->g[nope->not_fixed_index[i]];
  }
}

void ppUHPROD (Nope *nope, int * status, const int * n, const _Bool * goth,
    const double * x, const double * vector, double * result) {
  int i = 0;
  ASSERT_NOPE_OK;
  for (i = 0; i < nope->nvar; i++)
    nope->workspace1[i] = 0.0;
  for (i = 0; i < *n; i++) {
    nope->x[nope->not_fixed_index[i]] = x[i];
    nope->workspace1[nope->not_fixed_index[i]] = vector[i];
  }
  if (nope->ncon == 0) {
    (*nope->origin_uhprod)(status, &nope->nvar, goth, nope->x,
        nope->workspace1, nope->workspace2);
  } else {
    // Remember that trivial constraints are linear, so they don't appear on the
    // Hessian
    (*nope->origin_chprod)(status, &nope->nvar, &nope->ncon, goth,
        nope->x, nope->y, nope->workspace1, nope->workspace2);
  }
  for (i = 0; i < *n; i++)
    result[i] = nope->workspace2[nope->not_fixed_index[i]];
}

void ppCNAMES (Nope *nope, int * status, const int * n, const int * m, char * pname, char *
    vnames, char * cnames) {
  char v[10*nope->nvar], c[10*nope->ncon];
  UNUSED(n);
  UNUSED(m);
  UNUSED(vnames);
  UNUSED(cnames);
  (*nope->origin_cnames)(status, &nope->nvar, &nope->ncon, pname, v, c);
  // We, knowingly, decided to not pass the names of the variables
}

void ppCFN (Nope *nope, int * status, const int * n, const int * m, const double * x,
    double * f, double * c) {
  int i = 0;
  ASSERT_NOPE_OK;
  for (i = 0; i < *n; i++)
    nope->x[nope->not_fixed_index[i]] = x[i];
  (*nope->origin_cfn)(status, &nope->nvar, &nope->ncon, nope->x, f, nope->c);
  for (i = 0; i < *m; i++)
    c[i] = nope->c[nope->not_trivial_index[i]];
}

void ppCOFG (Nope *nope, int * status, const int * n, const double * x, double *
    f, double * g, _Bool * grad) {
  int i = 0;
  ASSERT_NOPE_OK;
  for (i = 0; i < *n; i++)
    nope->x[nope->not_fixed_index[i]] = x[i];
  (*nope->origin_cofg)(status, &nope->nvar, nope->x, f, nope->g, grad);
  if (*grad) {
    for (i = 0; i < *n; i++)
      g[i] = nope->g[nope->not_fixed_index[i]];
  }
}

void ppCHPROD (Nope *nope, int * status, const int * n, const int * m, const _Bool *
    goth, const double * x, const double * y, double * vector, double * result) {
  int i = 0;
  ASSERT_NOPE_OK;
  for (i = 0; i < nope->nvar; i++)
    nope->workspace1[i] = 0.0;
  for (i = 0; i < *n; i++) {
    nope->x[nope->not_fixed_index[i]] = x[i];
    nope->workspace1[nope->not_fixed_index[i]] = vector[i];
  }
  // Remember that trivial constraints are linear, so they don't appear on the
  // Hessian
  for (i = 0; i < *m; i++) {
    nope->y[nope->not_trivial_index[i]] = y[i];
  }
  (*nope->origin_chprod)(status, &nope->nvar, &nope->ncon, goth,
      nope->x, nope->y, nope->workspace1, nope->workspace2);
  for (i = 0; i < *n; i++)
    result[i] = nope->workspace2[nope->not_fixed_index[i]];
}

void ppCCFSG (Nope *nope, int * status, const int * n, const int * m, const
    double * x, double * c, int * nnzj, const int * lj, double * Jval, int *
    Jvar, int * Jfun, const _Bool * grad) {
  int i, k;
  ASSERT_NOPE_OK;
  for (i = 0; i < *n; i++)
    nope->x[nope->not_fixed_index[i]] = x[i];
  (*nope->origin_ccfsg)(status, &nope->nvar, &nope->ncon, nope->x, nope->c,
      nnzj, lj, Jval, Jvar, Jfun, grad);
  for (i = 0; i < *m; i++)
    c[i] = nope->c[nope->not_trivial_index[i]];
  if (!(*grad))
    return;
  // First the fixed variables are removed. This leaves some empty
  // columns
  for (k = 0; k < *nnzj; k++) {
    if (nope->is_fixed[Jvar[k]-1] || nope->is_trivial[Jfun[k]-1]) {
      (*nnzj)--;
      if (k < *nnzj) {
        Jval[k] = Jval[*nnzj];
        Jvar[k] = Jvar[*nnzj];
        Jfun[k] = Jfun[*nnzj];
        k--;
      }
    }
  }
  // Now we reorder the columns, reducing by one for each fixed
  // variable before that column
  for (i = nope->nfix-1; i >= 0; i--) {
    for (k = 0; k < *nnzj; k++) {
      if (Jvar[k]-1 > nope->fixed_index[i]) {
        Jvar[k]--;
      }
    }
  }

  for (i = nope->ntrivial-1; i >= 0; i--) {
    for (k = 0; k < *nnzj; k++) {
      if (Jfun[k]-1 > nope->trivial_index[i]) {
        Jfun[k]--;
      }
    }
  }
}
