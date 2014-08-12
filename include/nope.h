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
#ifndef nope_h
#define nope_h

#include "nope_definitions.h"

typedef struct _nope {
  // Setup function pointers
  pusetup origin_usetup;
  pcsetup origin_csetup;
  // Unprocessed function pointers
  punames origin_unames;
  pufn    origin_ufn;
  puofg   origin_uofg;
  puhprod origin_uhprod;
  pcnames origin_cnames;
  pcdimen origin_cdimen;
  pcdimsj origin_cdimsj;
  pcfn    origin_cfn;
  pcofg   origin_cofg;
  pchprod origin_chprod;
  pccfsg  origin_ccfsg;
  // Information
  int status;
  _Bool constrained;
  int nvar, ncon, nfix, ntrivial, nlinear, jmax;
  // Data
  int *fixed_index, *trivial_index;
  int *not_fixed_index, *not_trivial_index;
  int nnzj;
  _Bool *is_fixed, *is_trivial;
  double *x, *g, *bl, *bu;
  double *c, *y, *cl, *cu;
  double *linbndl, *linbndu;
  double *Jval;
  int *Jvar, *Jfun;
  _Bool *equatn, *linear;
  double *workspace1, *workspace2;
} Nope;

// {Cons,Des}tructor
Nope * initializeNope ();
void destroyNope (Nope *);

// Setup functions
void setFuncs (Nope *, pcdimen, pusetup, punames, pufn, puofg, puhprod, pcsetup,
    pcnames, pcfn, pcofg, pchprod, pccfsg, pcdimsj);
void runUncSetup (Nope *, int *, double *, double *, double *);
void runConSetup (Nope *, int *, double *, double *, double *,
    int *, double *, double *, double *, _Bool *, _Bool *, int *);

// Debug Functions
void printJacobian (int, int, int, double *, int *, int *);

// Run the processor
int runNope (Nope *);
void findFixedVariables (Nope *);
void findTrivialConstraints (Nope *);

// Interface functions
void ppDIMEN (Nope *nope,
              int * nvar,
              int * ncon);

void ppUNAMES (Nope *nope,
               int * status,
               int * n,
               char * pname,
               char * vnames);

void ppUFN (Nope *nope,
            int * status,  // (out) Output status
            int * n,       // (in)  Number of variables
            double * x,      // (in)  Current estimate of the solution
            double * f);     // (out) Objective function evaluated at x

void ppUOFG (Nope *nope,
             int * status,
             int * n,
             double * x,
             double * f,
             double * g,
             _Bool * grad);

void ppUHPROD (Nope *nope,
               int * status,
               int * n,
               _Bool * goth,
               double * x,
               double * vector,
               double * result);

void ppCNAMES (Nope *nope,
               int * status,
               int * n,
               int * m,
               char * pname,
               char * vnames,
               char * cnames);

void ppCFN (Nope *nope,
            int * status,
            int * n,
            int * m,
            double * x,
            double * f,
            double * c);

void ppCOFG (Nope *nope,
             int * status,
             int * n,
             double * x,
             double * f,
             double * g,
             _Bool * grad);

void ppCHPROD (Nope *nope,
               int * status,
               int * n,
               int * m,
               _Bool * goth,
               double * x,
               double * y,
               double * vector,
               double * result);

void ppCCFSG (Nope *nope,
              int * status,
              int * n,
              int * m,
              double * x,
              double * c,
              int * nnzj,
              int * lj,
              double * Jval,
              int * Jvar,
              int * Jfun,
              _Bool * grad);
#endif
