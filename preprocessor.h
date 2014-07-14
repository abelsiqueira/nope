NOPE
Copyright (C) 2014  Abel Soares Siqueira

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
#ifndef preprocessor_h
#define preprocessor_h

#include "preprocessor_definitions.h"

typedef struct _preprocessor {
  // Setup function pointers
  pusetup origin_usetup;
  pcsetup origin_csetup;
  // Unprocessed function pointers
  pufn    origin_ufn;
  puofg   origin_uofg;
  puhprod origin_uhprod;
  pcdimen origin_cdimen;
  pcdimsj origin_cdimsj;
  pcfn    origin_cfn;
  pcofg   origin_cofg;
  pchprod origin_chprod;
  pccfsg  origin_ccfsg;
  // Information
  int status;
  _Bool constrained;
  int nvar, ncon, nfix, ntrivial, jmax;
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
} Preprocessor;

// {Cons,Des}tructor
Preprocessor * initializePreprocessor ();
void destroyPreprocessor (Preprocessor *);

// Setup functions
void setFuncs (Preprocessor *, pcdimen, pusetup, pufn, puofg, puhprod,
    pcsetup, pcfn, pcofg, pchprod, pccfsg, pcdimsj);
void runUncSetup (Preprocessor *, int *, double *, double *, double *);
void runConSetup (Preprocessor *, int *, double *, double *, double *,
    int *, double *, double *, double *, _Bool *, _Bool *, int *);

// Debug Functions
void printJacobian (int, int, int, double *, int *, int *);

// Run the processor
int runPreprocessor (Preprocessor *);
void findFixedVariables (Preprocessor *);
void findTrivialConstraints (Preprocessor *);

// Interface functions
void ppDIMEN (Preprocessor *prep,
              int * nvar,
              int * ncon);

void ppUFN (Preprocessor *prep,
            int * status,  // (out) Output status
            int * n,       // (in)  Number of variables
            double * x,      // (in)  Current estimate of the solution
            double * f);     // (out) Objective function evaluated at x

void ppUOFG (Preprocessor *prep,
             int * status,
             int * n,
             double * x,
             double * f,
             double * g,
             _Bool * grad);

void ppUHPROD (Preprocessor *prep,
               int * status,
               int * n,
               _Bool * goth,
               double * x,
               double * vector,
               double * result);

void ppCFN (Preprocessor *prep,
            int * status,
            int * n,
            int * m,
            double * x,
            double * f,
            double * c);

void ppCOFG (Preprocessor *prep,
             int * status,
             int * n,
             double * x,
             double * f,
             double * g,
             _Bool * grad);

void ppCHPROD (Preprocessor *prep,
               int * status,
               int * n,
               int * m,
               _Bool * goth,
               double * x,
               double * y,
               double * vector,
               double * result);

void ppCCFSG (Preprocessor *prep,
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
