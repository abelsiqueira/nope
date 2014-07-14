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
#include "preprocessor.h"

void ppDIMEN (Preprocessor *prep, int *nvar, int *ncon) {
  if (prep == 0 || prep->status != STATUS_OK)
    return;
  *nvar = prep->nvar - prep->nfix;
  *ncon = prep->ncon - prep->ntrivial;
}

void ppUFN (Preprocessor *prep, int * status, int * n, double * x, double * f) {
  int i = 0;
  if (prep == 0 || prep->status != STATUS_OK)
    return;
  for (i = 0; i < *n; i++)
    prep->x[prep->not_fixed_index[i]] = x[i];
  (*prep->origin_ufn)(status, &prep->nvar, prep->x, f);
}

void ppUOFG (Preprocessor *prep, int * status, int * n, double * x, double *
    f, double * g, _Bool * grad) {
  int i = 0;
  if (prep == 0 || prep->status != STATUS_OK)
    return;
  for (i = 0; i < *n; i++)
    prep->x[prep->not_fixed_index[i]] = x[i];
  (*prep->origin_uofg)(status, &prep->nvar, prep->x, f, prep->g, grad);
  if (*grad) {
    for (i = 0; i < *n; i++)
      g[i] = prep->g[prep->not_fixed_index[i]];
  }
}

void ppUHPROD (Preprocessor *prep, int * status, int * n, _Bool * goth,
    double * x, double * vector, double * result) {
  int i = 0;
  if (prep == 0 || prep->status != STATUS_OK)
    return;
  for (i = 0; i < prep->nvar; i++)
    prep->workspace1[i] = 0.0;
  for (i = 0; i < *n; i++) {
    prep->x[prep->not_fixed_index[i]] = x[i];
    prep->workspace1[prep->not_fixed_index[i]] = vector[i];
  }
  (*prep->origin_uhprod)(status, &prep->nvar, goth, prep->x,
      prep->workspace1, prep->workspace2);
  for (i = 0; i < *n; i++)
    result[i] = prep->workspace2[prep->not_fixed_index[i]];
}

void ppCFN (Preprocessor *prep, int * status, int * n, int * m, double * x,
    double * f, double * c) {
  int i = 0;
  if (prep == 0 || prep->status != STATUS_OK)
    return;
  for (i = 0; i < *n; i++)
    prep->x[prep->not_fixed_index[i]] = x[i];
  (*prep->origin_cfn)(status, &prep->nvar, &prep->ncon, prep->x, f,
      prep->c);
  for (i = 0; i < *m; i++)
    c[i] = prep->c[prep->not_trivial_index[i]];
}

void ppCOFG (Preprocessor *prep, int * status, int * n, double * x, double *
    f, double * g, _Bool * grad) {
  int i = 0;
  if (prep == 0 || prep->status != STATUS_OK)
    return;
  for (i = 0; i < *n; i++)
    prep->x[prep->not_fixed_index[i]] = x[i];
  (*prep->origin_cofg)(status, &prep->nvar, prep->x, f, prep->g, grad);
  if (*grad) {
    for (i = 0; i < *n; i++)
      g[i] = prep->g[prep->not_fixed_index[i]];
  }
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
  for (i = 0; i < *m; i++) {
    prep->y[prep->not_trivial_index[i]] = y[i];
  }
  (*prep->origin_chprod)(status, &prep->nvar, &prep->ncon, goth,
      prep->x, prep->y, prep->workspace1, prep->workspace2);
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
  (*prep->origin_ccfsg)(status, &prep->nvar, &prep->ncon, prep->x, prep->c,
      nnzj, lj, Jval, Jvar, Jfun, grad);
  for (i = 0; i < *m; i++)
    c[i] = prep->c[prep->not_trivial_index[i]];
  if (!(*grad))
    return;
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
      if (Jfun[k]-1 > prep->trivial_index[i]) {
        Jfun[k]--;
      }
    }
  }
}
