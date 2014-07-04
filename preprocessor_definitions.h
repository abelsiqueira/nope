#ifndef preprocessor_definitions_h
#define preprocessor_definitions_h

#define UNUSED(x) (void)(x)

#include <stdbool.h>

#define STATUS_UNINITIALIZED -1
#define STATUS_OK 0
#define STATUS_ERROR 1

#define PPEPS 1e-12

// Unconstrained
typedef void (*pusetup) (int *, int *, int *, int *, int *, double *, double *,
                         double *);
typedef void (*pufn)    (int *, int *, double *, double *);
typedef void (*puofg)   (int *, int *, double *, double *, double *, _Bool *);
typedef void (*puhprod) (int *, int *, _Bool *, double *, double *, double *);
// Constrained
typedef void (*pcsetup) (int *, int *, int *, int *, int *, int *, double *,
                         double *, double *, double *, double *, double *, _Bool *, _Bool *,
                         int *, int *, int *);
typedef void (*pcfn)    (int *, int *, int *, double *, double *, double *);
typedef void (*pcofg)   (int *, int *, double *, double *, double *, _Bool *);
typedef void (*pccfsg)  (int *, int *, int *, double *, double *, int *, int *,
                         double *, int *, int *, _Bool *);
typedef void (*pchprod) (int *, int *, int *, _Bool *, double *, double *, double *,
                         double *);

#endif
