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
extern Nope *nope;

void unames (int *status, const int *n, char *pname, char *vnames) {
  ppUNAMES(nope, status, n, pname, vnames);
}

void ufn (int *status, const int *n, const double *x, double *f) {
  ppUFN(nope, status, n, x, f);
}

void uofg (int * status, const int * n, const double * x,
    double * f, double * g, const _Bool * grad) {
  ppUOFG(nope, status, n, x, f, g, grad);
}

void uhprod (int * status, const int * n, const _Bool *
    goth, const double * x, const double * vector, double * result) {
  ppUHPROD(nope, status, n, goth, x, vector, result);
}

void cnames (int *status, const int *n, const int *m, char *pname, char *vnames, char
    *cnames) {
  ppCNAMES(nope, status, n, m, pname, vnames, cnames);
}

void cfn (int * status, const int * n, const int * m, const double
    * x, double * f, double * c) {
  ppCFN(nope, status, n, m, x, f, c);
}

void cofg (int * status, const int * n, const double * x,
    double * f, double * g, _Bool * grad) {
  ppCOFG(nope, status, n, x, f, g, grad);
}

void chprod (int * status, const int * n, const int * m,
    const _Bool * goth, const double * x, const double * y, double * vector,
    double * result) {
  ppCHPROD(nope, status, n, m, goth, x, y, vector, result);
}

void ccfsg (int * status, const int * n, const int * m, const double
    * x, double * c, int * nnzj, const int * lj, double * Jval, int * Jvar,
    int * Jfun, const _Bool * grad) {
  ppCCFSG(nope, status, n, m, x, c, nnzj, lj, Jval, Jvar, Jfun, grad);
}

