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
#include <stdio.h>
#include <stdlib.h>

Nope * initializeNope () {
  Nope * nope = (Nope *) malloc (sizeof(Nope));

  nope->origin_usetup = 0;
  nope->origin_csetup = 0;
  nope->origin_ufn = 0;
  nope->origin_uofg = 0;
  nope->origin_uhprod = 0;
  nope->origin_cdimen = 0;
  nope->origin_cdimsj = 0;
  nope->origin_cfn = 0;
  nope->origin_cofg = 0;
  nope->origin_chprod = 0;
  nope->origin_ccfsg = 0;
  nope->status = STATUS_UNINITIALIZED;
  nope->constrained = 0;

  return nope;
}

void destroyNope (Nope *nope) {
  if (nope != 0)
    free(nope);
}

void setFuncs (Nope *nope, pcdimen cdimen, pusetup usetup,
    pufn ufn, puofg uofg, puhprod uhprod, pcsetup csetup, pcfn cfn,
    pcofg cofg, pchprod chprod, pccfsg ccfsg, pcdimsj cdimsj) {
  nope->origin_cdimen = cdimen;
  nope->origin_usetup = usetup;
  nope->origin_ufn = ufn;
  nope->origin_uofg = uofg;
  nope->origin_uhprod = uhprod;
  nope->origin_csetup = csetup;
  nope->origin_cfn = cfn;
  nope->origin_cofg = cofg;
  nope->origin_chprod = chprod;
  nope->origin_ccfsg = ccfsg;
  nope->origin_cdimsj = cdimsj;

  nope->status = STATUS_READY;
}

int runNope (Nope *nope) {
  int status = 0;
  int funit = 42, fout = 6, io_buffer = 11;
  int efirst = 0, lfirst = 0, nvfirst = 0;
  _Bool grad = true;
  int i, j;
  int knotfixed = 0, knottrivial = 0;

  if (nope == 0 || nope->status != STATUS_READY) {
    return 1;
  }

  (*nope->origin_cdimen)(&status, &funit, &nope->nvar, &nope->ncon);

  nope->nfix = 0;
  nope->ntrivial = 0;
  nope->fixed_index = (int *) malloc(nope->nvar*sizeof(int));
  nope->not_fixed_index = (int *) malloc(nope->nvar*sizeof(int));
  nope->is_fixed = (_Bool *) malloc(nope->nvar*sizeof(_Bool));
  nope->trivial_index = (int *) malloc(nope->ncon*sizeof(int));
  nope->not_trivial_index = (int *) malloc(nope->ncon*sizeof(int));
  nope->is_trivial = (_Bool *) malloc(nope->ncon*sizeof(_Bool));
  nope->x = (double *) malloc(nope->nvar*sizeof(double));
  nope->g = (double *) malloc(nope->nvar*sizeof(double));
  nope->bl = (double *) malloc(nope->nvar*sizeof(double));
  nope->bu = (double *) malloc(nope->nvar*sizeof(double));
  if (nope->ncon > 0) {
    nope->c = (double *) malloc(nope->ncon*sizeof(double));
    nope->y = (double *) malloc(nope->ncon*sizeof(double));
    nope->cl = (double *) malloc(nope->ncon*sizeof(double));
    nope->cu = (double *) malloc(nope->ncon*sizeof(double));
    nope->linbndl = (double *) malloc(nope->ncon*sizeof(double));
    nope->linbndu = (double *) malloc(nope->ncon*sizeof(double));
    nope->equatn = (_Bool *) malloc(nope->ncon*sizeof(_Bool));
    nope->linear = (_Bool *) malloc(nope->ncon*sizeof(_Bool));
  }
  nope->workspace1 = (double *) malloc(nope->nvar*sizeof(double));
  nope->workspace2 = (double *) malloc(nope->nvar*sizeof(double));

  if (nope->ncon == 0) {
    (nope->origin_usetup)(&status, &funit, &fout, &io_buffer,
        &nope->nvar, nope->x, nope->bl, nope->bu);
  } else {
    (nope->origin_csetup)(&status, &funit, &fout, &io_buffer,
        &nope->nvar, &nope->ncon, nope->x, nope->bl, nope->bu,
        nope->y, nope->cl, nope->cu, nope->equatn, nope->linear,
        &efirst, &lfirst, &nvfirst);
  }

  if (nope->ncon > 0)
    (*nope->origin_cdimsj)(&status, &nope->jmax);
  else
    nope->jmax = 0;
  if (nope->ncon > 0) {
    nope->Jval = (double *) malloc(nope->jmax*sizeof(double));
    nope->Jvar = (int *) malloc(nope->jmax*sizeof(int));
    nope->Jfun = (int *) malloc(nope->jmax*sizeof(int));
  }

  nope->status = STATUS_PROCESSING;
  while (nope->status == STATUS_PROCESSING) {
    (*nope->origin_ccfsg)(&status, &nope->nvar, &nope->ncon,
        nope->x, nope->c, &nope->nnzj, &nope->jmax, nope->Jval,
        nope->Jvar, nope->Jfun, &grad);
    for (i = 0; i < nope->ncon; i++) {
      nope->linbndl[i] = nope->cl[i] - nope->c[i];
      nope->linbndu[i] = nope->cu[i] - nope->c[i];
    }
    for (i = 0; i < nope->nnzj; i++) {
      j = nope->Jfun[i]-1;
      if (nope->linear[j]) {
        nope->linbndl[j] += nope->Jval[i]*nope->x[nope->Jvar[i]-1];
        nope->linbndu[j] += nope->Jval[i]*nope->x[nope->Jvar[i]-1];
      }
    }
    nope->status = STATUS_PROCESSED;
    findFixedVariables(nope);
    findTrivialConstraints(nope);
  }

  nope->nfix = 0;
  knotfixed = 0;
  nope->ntrivial = 0;
  knottrivial = 0;
  for (i = 0; i < nope->nvar; i++) {
    if (nope->is_fixed[i])
      nope->fixed_index[nope->nfix++] = i;
    else
      nope->not_fixed_index[knotfixed++] = i;
  }

  for (i = 0; i < nope->ncon; i++) {
    if (nope->is_trivial[i]) {
      nope->trivial_index[nope->ntrivial++] = i;
      nope->y[i] = 0;
    } else
      nope->not_trivial_index[knottrivial++] = i;
  }

  nope->status = STATUS_OK;
  return 0;
}

void runUncSetup (Nope *nope, int *nvar, double *x, double
    *bl, double *bu) {
  int i;
  for (i = 0; i < *nvar; i++) {
    x[i] = nope->x[nope->not_fixed_index[i]];
    bl[i] = nope->bl[nope->not_fixed_index[i]];
    bu[i] = nope->bu[nope->not_fixed_index[i]];
  }
}

void runConSetup (Nope *nope, int *nvar, double *x, double
    *bl, double *bu, int *ncon, double *y, double *cl, double *cu,
    _Bool *equatn, _Bool *linear, int *jmax) {
  int i, j;

  for (i = 0; i < *nvar; i++) {
    j = nope->not_fixed_index[i];
    x[i] = nope->x[j];
    bl[i] = nope->bl[j];
    bu[i] = nope->bu[j];
  }
  for (i = 0; i < *ncon; i++) {
    j = nope->not_trivial_index[i];
    y[i] = nope->y[j];
    cl[i] = nope->cl[j];
    cu[i] = nope->cu[j];
    equatn[i] = nope->equatn[j];
    linear[i] = nope->linear[j];
  }
  *jmax = nope->jmax;
}

void findFixedVariables (Nope *nope) {
  int i, status = 0;
  double f = 0;
  for (i = 0; i < nope->nvar; i++) {
    if (nope->is_fixed[i])
      continue;
    if (nope->bu[i] - nope->bl[i] < PPEPS) {
      nope->status = STATUS_PROCESSING;
      nope->is_fixed[i] = 1;
      nope->x[i] = (nope->bl[i] + nope->bu[i])/2;
    } else {
      nope->is_fixed[i] = 0;
    }
  }

  (*nope->origin_cfn)(&status, &nope->nvar, &nope->ncon, nope->x, &f, nope->c);

  for (i = 0; i < nope->nnzj; i++) {
    if (nope->is_fixed[nope->Jvar[i]-1]) {
      if (nope->linear[nope->Jfun[i]-1]) {
        nope->linbndl[nope->Jfun[i]-1] -= nope->Jval[i]*nope->x[nope->Jvar[i]-1];
        nope->linbndu[nope->Jfun[i]-1] -= nope->Jval[i]*nope->x[nope->Jvar[i]-1];
      }
      nope->Jval[i] = nope->Jval[nope->nnzj-1];
      nope->Jvar[i] = nope->Jvar[nope->nnzj-1];
      nope->Jfun[i] = nope->Jfun[nope->nnzj-1];
      nope->nnzj--;
      i--;
    }
  }
}

void findTrivialConstraints (Nope *nope) {
  int i, j, k;
  int nnzj = nope->nnzj;
  int ncon = nope->ncon;
  int nnzj_per_line[ncon];
  int index_of_nnzj[ncon];
  double newlim;

  for (i = 0; i < ncon; i++) {
    nnzj_per_line[i] = 0;
    index_of_nnzj[i] = -1;
  }
  for (i = 0; i < nnzj; i++) {
    nnzj_per_line[nope->Jfun[i]-1]++;
    index_of_nnzj[nope->Jfun[i]-1] = i;
  }

  for (i = 0; i < ncon; i++) {
    if (nope->is_trivial[i])
      continue;
    if (nope->linear[nope->not_trivial_index[i]]) {
      if (nnzj_per_line[i] == 0) {
        nope->status = STATUS_PROCESSING;
        nope->is_trivial[i] = true;
      } else if (nnzj_per_line[i] == 1) {
        nope->status = STATUS_PROCESSING;
        nope->is_trivial[i] = true;
        j = index_of_nnzj[i];
        k = nope->Jvar[j]-1;
        if (nope->equatn[i]) {
          nope->x[k] = nope->linbndl[nope->Jfun[j]-1]/nope->Jval[j];
          nope->is_fixed[k] = true;
        } else {
          // cl <= a x_j <= cu
          // cl/a <= x_j <= cu/a
          newlim = nope->linbndu[nope->Jfun[j]-1]/nope->Jval[j];
          if (newlim < nope->bu[k])
            nope->bu[k] = newlim;
          newlim = nope->linbndl[nope->Jfun[j]-1]/nope->Jval[j];
          if (newlim > nope->bl[k])
            nope->bl[k] = newlim;

          if (nope->bu[k] - nope->bl[k] < PPEPS) {
            nope->x[k] = (nope->bl[k] + nope->bu[k])/2;
            nope->is_fixed[k] = true;
          }
        }
      }
    }
  }
}

void printJacobian (int ncon, int nvar, int nnzj, double *Jval, int *Jvar,
    int *Jfun) {
  int i, j;
  int n = nvar*ncon;
  double M[n];

  for (i = 0; i < n; i++)
    M[i] = 0;

  for (i = 0; i < nnzj; i++)
    M[Jvar[i]-1 + (Jfun[i]-1)*nvar] = Jval[i];

  for (i = 0; i < ncon; i++) {
    for (j = 0; j < nvar; j++)
      printf("%3.2lf ", M[j + i*nvar]);
    printf("\n");
  }
}
