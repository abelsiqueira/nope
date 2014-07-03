#include "preprocessor.h"
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
  prep->nvar = nvar;
  prep->ncon = 0;
  findFixedVariable(prep, nvar, x, bl, bu);
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
  prep->nvar = nvar;
  prep->ncon = ncon;
  findFixedVariable(prep, nvar, x, bl, bu);
}

void findFixedVariables (Preprocessor *prep, int nvar, double*x, double *bl)

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
  if (prep == 0 || prep->status != STATUS_OK)
    return;
  (*prep->origin_cfn)(status, n, m, x, f, c);
}

void ppCOFG (Preprocessor *prep, int * status, int * n, double * x, double *
    f, double * g, _Bool * grad) {
  if (prep == 0 || prep->status != STATUS_OK)
    return;
  (*prep->origin_cofg)(status, n, x, f, g, grad);
}

void ppCHPROD (Preprocessor *prep, int * status, int * n, int * m, _Bool *
    goth, double * x, double * y, double * vector, double * result) {
  if (prep == 0 || prep->status != STATUS_OK)
    return;
  (*prep->origin_chprod)(status, n, m, goth, x, y, vector, result);
}

void ppCCFSG (Preprocessor *prep, int * status, int * n, int * m, double *
    x, double * c, int * nnzj, int * lj, double * Jval, int * Jvar, int * Jfun,
    _Bool * grad) {
  if (prep == 0 || prep->status != STATUS_OK)
    return;
  (*prep->origin_ccfsg)(status, n, m, x, c, nnzj, lj, Jval, Jvar,
      Jfun, grad);
}
