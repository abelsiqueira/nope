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
  prep->origin_cdimen = 0;
  prep->origin_cdimsj = 0;
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

void setFuncs (Preprocessor *prep, pcdimen cdimen, pusetup usetup,
    pufn ufn, puofg uofg, puhprod uhprod, pcsetup csetup, pcfn cfn,
    pcofg cofg, pchprod chprod, pccfsg ccfsg, pcdimsj cdimsj) {
  prep->origin_cdimen = cdimen;
  prep->origin_usetup = usetup;
  prep->origin_ufn = ufn;
  prep->origin_uofg = uofg;
  prep->origin_uhprod = uhprod;
  prep->origin_csetup = csetup;
  prep->origin_cfn = cfn;
  prep->origin_cofg = cofg;
  prep->origin_chprod = chprod;
  prep->origin_ccfsg = ccfsg;
  prep->origin_cdimsj = cdimsj;

  prep->status = STATUS_READY;
}

int runPreprocessor (Preprocessor *prep) {
  int status = 0;
  int funit = 42, fout = 6, io_buffer = 11;
  int efirst = 0, lfirst = 0, nvfirst = 0;
  _Bool grad = true;
  int i;
  int knotfixed = 0, knottrivial = 0;

  if (prep == 0 || prep->status != STATUS_READY) {
    return 1;
  }

  (*prep->origin_cdimen)(&status, &funit, &prep->nvar, &prep->ncon);

  prep->nfix = 0;
  prep->ntrivial = 0;
  prep->fixed_index = (int *) malloc(prep->nvar*sizeof(int));
  prep->not_fixed_index = (int *) malloc(prep->nvar*sizeof(int));
  prep->is_fixed = (_Bool *) malloc(prep->nvar*sizeof(_Bool));
  prep->trivial_index = (int *) malloc(prep->ncon*sizeof(int));
  prep->not_trivial_index = (int *) malloc(prep->ncon*sizeof(int));
  prep->is_trivial = (_Bool *) malloc(prep->ncon*sizeof(_Bool));
  prep->x = (double *) malloc(prep->nvar*sizeof(double));
  prep->g = (double *) malloc(prep->nvar*sizeof(double));
  prep->bl = (double *) malloc(prep->nvar*sizeof(double));
  prep->bu = (double *) malloc(prep->nvar*sizeof(double));
  if (prep->ncon > 0) {
    prep->c = (double *) malloc(prep->ncon*sizeof(double));
    prep->y = (double *) malloc(prep->ncon*sizeof(double));
    prep->cl = (double *) malloc(prep->ncon*sizeof(double));
    prep->cu = (double *) malloc(prep->ncon*sizeof(double));
    prep->equatn = (_Bool *) malloc(prep->ncon*sizeof(_Bool));
    prep->linear = (_Bool *) malloc(prep->ncon*sizeof(_Bool));
  }
  prep->workspace1 = (double *) malloc(prep->nvar*sizeof(double));
  prep->workspace2 = (double *) malloc(prep->nvar*sizeof(double));

  if (prep->ncon == 0) {
    (prep->origin_usetup)(&status, &funit, &fout, &io_buffer,
        &prep->nvar, prep->x, prep->bl, prep->bu);
  } else {
    (prep->origin_csetup)(&status, &funit, &fout, &io_buffer,
        &prep->nvar, &prep->ncon, prep->x, prep->bl, prep->bu,
        prep->y, prep->cl, prep->cu, prep->equatn, prep->linear,
        &efirst, &lfirst, &nvfirst);
  }

  if (prep->ncon > 0)
    (*prep->origin_cdimsj)(&status, &prep->jmax);
  else
    prep->jmax = 0;
  if (prep->ncon > 0) {
    prep->Jval = (double *) malloc(prep->jmax*sizeof(double));
    prep->Jvar = (int *) malloc(prep->jmax*sizeof(int));
    prep->Jfun = (int *) malloc(prep->jmax*sizeof(int));
  }

  (*prep->origin_ccfsg)(&status, &prep->nvar, &prep->ncon,
      prep->x, prep->y, &prep->nnzj, &prep->jmax, prep->Jval,
      prep->Jvar, prep->Jfun, &grad);
  printJacobian (prep->ncon, prep->nvar, prep->nnzj, prep->Jval,
      prep->Jvar, prep->Jfun);
  findFixedVariables(prep);
  printJacobian (prep->ncon, prep->nvar, prep->nnzj, prep->Jval,
      prep->Jvar, prep->Jfun);
  findTrivialConstraints(prep);
  printJacobian (prep->ncon, prep->nvar, prep->nnzj, prep->Jval,
      prep->Jvar, prep->Jfun);

  for (i = 0; i < prep->nvar; i++) {
    if (prep->is_fixed[i])
      prep->fixed_index[prep->nfix++] = i;
    else
      prep->not_fixed_index[knotfixed++] = i;
  }

  for (i = 0; i < prep->ncon; i++) {
    if (prep->is_trivial[i])
      prep->trivial_index[prep->ntrivial++] = i;
    else
      prep->not_trivial_index[knottrivial++] = i;
  }

  prep->status = STATUS_OK;
  return 0;
}

void runUncSetup (Preprocessor *prep, int *nvar, double *x, double
    *bl, double *bu) {
  int i;
  for (i = 0; i < *nvar; i++) {
    x[i] = prep->x[prep->not_fixed_index[i]];
    bl[i] = prep->bl[prep->not_fixed_index[i]];
    bu[i] = prep->bu[prep->not_fixed_index[i]];
  }
}

void runConSetup (Preprocessor *prep, int *nvar, double *x, double
    *bl, double *bu, int *ncon, double *y, double *cl, double *cu,
    _Bool *equatn, _Bool *linear, int *jmax) {
  int i, j;

  for (i = 0; i < *nvar; i++) {
    j = prep->not_fixed_index[i];
    x[i] = prep->x[j];
    bl[i] = prep->bl[j];
    bu[i] = prep->bu[j];
  }
  for (i = 0; i < *ncon; i++) {
    j = prep->not_fixed_index[i];
    y[i] = prep->y[j];
    cl[i] = prep->cl[j];
    cu[i] = prep->cu[j];
    equatn[i] = prep->equatn[j];
    linear[i] = prep->linear[j];
  }
  *jmax = prep->jmax;
}

void findFixedVariables (Preprocessor *prep) {
  int i;
  for (i = 0; i < prep->nvar; i++) {
    if (prep->bu[i] - prep->bl[i] < PPEPS) {
      prep->is_fixed[i] = 1;
      prep->x[i] = (prep->bl[i] + prep->bu[i])/2;
    } else {
      prep->is_fixed[i] = 0;
    }
  }

  for (i = 0; i < prep->nnzj; i++) {
    if (prep->is_fixed[prep->Jvar[i]-1]) {
      if (prep->linear[prep->Jfun[i]-1]) {
        prep->c[prep->Jfun[i]-1] -= prep->Jval[i];
      }
      prep->Jval[i] = prep->Jval[prep->nnzj-1];
      prep->Jvar[i] = prep->Jvar[prep->nnzj-1];
      prep->Jfun[i] = prep->Jfun[prep->nnzj-1];
      prep->nnzj--;
      i--;
    }
  }
}

void findTrivialConstraints (Preprocessor *prep) {
  int i, j;
  int nnzj = prep->nnzj;
  int ncon = prep->ncon - prep->ntrivial;
  int nnzj_per_line[ncon];
  int index_of_nnzj[ncon];

  for (i = 0; i < ncon; i++) {
    nnzj_per_line[i] = 0;
    index_of_nnzj[i] = -1;
  }
  for (i = 0; i < nnzj; i++) {
    nnzj_per_line[prep->Jfun[i]-1]++;
    index_of_nnzj[prep->Jfun[i]-1] = i;
  }

  for (i = 0; i < ncon; i++) {
    if (prep->linear[prep->not_trivial_index[i]]) {
      if (nnzj_per_line[i] == 0)
        prep->is_trivial[i] = true;
      else if (nnzj_per_line[i] == 1) {
        prep->is_trivial[i] = true;
        j = index_of_nnzj[i];
        prep->x[prep->Jvar[j]-1] = prep->c[prep->Jfun[j]-1]/prep->Jval[j];
        prep->is_fixed[i] = true;
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
