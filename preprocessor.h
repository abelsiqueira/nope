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
  pcfn    origin_cfn;
  pcofg   origin_cofg;
  pchprod origin_chprod;
  pccfsg  origin_ccfsg;
  // Information
  int status;
  _Bool constrained;
  int nvar, ncon, nfix;
  // Data
  int *fixed_index;
  int *not_fixed_index;
  _Bool *is_fixed;
  double *x, *g;
  double *workspace1, *workspace2;
} Preprocessor;

// {Cons,Des}tructor
Preprocessor * initializePreprocessor ();
void destroyPreprocessor (Preprocessor *);

// Setup functions
void setUncFuncs (Preprocessor *, pusetup, pufn, puofg, puhprod);
void setConFuncs (Preprocessor *, pcsetup, pcfn, pcofg, pchprod, pccfsg);
void runUncSetup (Preprocessor *, int *, double *, double *, double *);
void runConSetup (Preprocessor *, int *, double *, double *, double *,
    int *, double *, double *, double *, _Bool *, _Bool *);
void findFixedVariables (Preprocessor *, int, double *, double *, double *);

// Debug Functions
void printJacobian (int, int, int, double *, int *, int *);

// Run the processor
int process (Preprocessor *);

// Interface functions
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
