#include "preprocessor.h"
#include <stdio.h>
#include <stdlib.h>

Preprocessor * initializePreprocessor () {
  Preprocessor * prep = (Preprocessor *) malloc (sizeof(Preprocessor));
  prep->origin_usetup = 0;
  prep->origin_csetup = 0;
  prep->origin_ufn = 0;
  prep->origin_uofg = 0;
  prep->origin_uhprod = 0;
  prep->origin_cfn = 0;
  prep->origin_cofg = 0;
  prep->origin_chprod = 0;
  prep->origin_ccfsg = 0;
  prep->status = STATUS_UNINITIALIZED;
  prep->constrained = 0;

  return prep;
}

void destroyPreprocessor (Preprocessor *prep) {
  if (prep != 0)
    free(prep);
}

void setUncFuncs (Preprocessor *prep, pusetup usetup, pufn ufn, puofg uofg,
    puhprod uhprod) {
  prep->origin_usetup = usetup;
  prep->origin_ufn = ufn;
  prep->origin_uofg = uofg;
  prep->origin_uhprod = uhprod;
  prep->status = STATUS_OK;
  prep->constrained = 0;
}

void setConFuncs (Preprocessor *prep, pcsetup csetup, pcfn cfn, pcofg cofg,
    pchprod chprod, pccfsg ccfsg) {
  prep->origin_csetup = csetup;
  prep->origin_cfn = cfn;
  prep->origin_cofg = cofg;
  prep->origin_chprod = chprod;
  prep->origin_ccfsg = ccfsg;
  prep->status = STATUS_OK;
  prep->constrained = 1;
}

int process (Preprocessor *prep) {
  if (prep == 0 || prep->status != STATUS_OK) {
    return 1;
  }
  return 0;
}

void runUncSetup (Preprocessor *prep, int *nvar, double *x, double
    *bl, double *bu) {
  int status = 0;
  int funit = 42, fout = 6, io_buffer = 11;

  (prep->origin_usetup)(&status, &funit, &fout, &io_buffer, nvar, x, bl, bu);
  prep->nvar = *nvar;
  prep->ncon = 0;
  findFixedVariables(prep, *nvar, x, bl, bu);
  *nvar = *nvar - prep->nfix;
}

void runConSetup (Preprocessor *prep, int *nvar, double *x, double *bl,
    double *bu, int *ncon, double *y, double *cl, double *cu,
    _Bool *equatn, _Bool *linear) {
  int status = 0;
  int funit = 42, fout = 6, io_buffer = 11;
  int efirst = 0, lfirst = 0, nvfirst = 0;

  (prep->origin_csetup)(&status, &funit, &fout, &io_buffer, nvar, ncon,
      x, bl, bu, y, cl, cu, equatn, linear, &efirst, &lfirst,
      &nvfirst);
  prep->nvar = *nvar;
  prep->ncon = *ncon;
  findFixedVariables(prep, *nvar, x, bl, bu);
  *nvar = *nvar - prep->nfix;
}

void findFixedVariables (Preprocessor *prep, int nvar, double *x,
    double *bl, double *bu) {
  int i, kfixed = 0, knotfixed = 0;

  prep->nfix = 0;
  prep->fixed_index = (int *) malloc(nvar*sizeof(int));
  prep->not_fixed_index = (int *) malloc(nvar*sizeof(int));
  prep->is_fixed = (_Bool *) malloc(nvar*sizeof(_Bool));
  prep->x = (double *) malloc(nvar*sizeof(double));
  prep->g = (double *) malloc(nvar*sizeof(double));
  prep->workspace1 = (double *) malloc(nvar*sizeof(double));
  prep->workspace2 = (double *) malloc(nvar*sizeof(double));

  for (i = 0; i < nvar; i++) {
    if (bu[i] - bl[i] < PPEPS) {
      prep->is_fixed[i] = 1;
      prep->x[i] = (bl[i]+bu[i])/2;
      prep->nfix++;
      prep->fixed_index[kfixed++] = i;
    } else {
      prep->is_fixed[i] = 0;
      prep->x[i] = x[i];
      prep->not_fixed_index[knotfixed++] = i;
    }
  }

  for (i = 0; i < nvar - prep->nfix; i++) {
    x[i] = prep->x[prep->not_fixed_index[i]];
    bl[i] = bl[prep->not_fixed_index[i]];
    bu[i] = bu[prep->not_fixed_index[i]];
  }
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
    if (prep->is_fixed[Jvar[k]-1]) {
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
  printJacobian (prep->ncon, *n, *nnzj, Jval, Jvar, Jfun);
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
