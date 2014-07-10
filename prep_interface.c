#include "preprocessor.h"
extern Preprocessor *prep;

void ufn (int *status, int *n, double *x, double *f) {
  ppUFN(prep, status, n, x, f);
}

void uofg (int * status, int * n, double * x,
    double * f, double * g, _Bool * grad) {
  ppUOFG(prep, status, n, x, f, g, grad);
}

void uhprod (int * status, int * n, _Bool *
    goth, double * x, double * vector, double * result) {
  ppUHPROD(prep, status, n, goth, x, vector, result);
}

void cfn (int * status, int * n, int * m, double
    * x, double * f, double * c) {
  ppCFN(prep, status, n, m, x, f, c);
}

void cofg (int * status, int * n, double * x,
    double * f, double * g, _Bool * grad) {
  ppCOFG(prep, status, n, x, f, g, grad);
}

void chprod (int * status, int * n, int * m,
    _Bool * goth, double * x, double * y, double * vector,
    double * result) {
  ppCHPROD(prep, status, n, m, goth, x, y, vector, result);
}

void ccfsg (int * status, int * n, int * m, double
    * x, double * c, int * nnzj, int * lj, double * Jval, int * Jvar,
    int * Jfun, _Bool * grad) {
  ppCCFSG(prep, status, n, m, x, c, nnzj, lj, Jval, Jvar, Jfun, grad);
}

