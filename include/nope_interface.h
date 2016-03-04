#ifndef nope_interface_h
#define nope_interface_h

// Unconstrained
typedef void (*pusetup) (int *, const int *, const int *, const int *, int *, double *, double *,
                         double *);
typedef void (*punames) (int *, const int *, char *, char *);
typedef void (*pufn)    (int *, const int *, const double *, double *);
typedef void (*puofg)   (int *, const int *, const double *, double *, double *, const _Bool *);
typedef void (*puhprod) (int *, const int *, const _Bool *, const double *, const double *, double *);
// Constrained
typedef void (*pcdimen) (int *, const int *, int *, int *);
typedef void (*pcdimsj) (int *, int *);
typedef void (*pcsetup) (int *, const int *, const int *, const int *, int *, int *, double *,
                         double *, double *, double *, double *, double *, _Bool
                         *, _Bool *, const int *, const int *, const int *);
typedef void (*pcnames) (int *, const int *, const int *, char *, char *, char *);
typedef void (*pcfn)    (int *, const int *, const int *, const double *, double *, double *);
typedef void (*pcofg)   (int *, const int *, const double *, double *, double *, _Bool *);
typedef void (*pccfsg)  (int *, const int *, const int *, const double *, double *, int *, const int *,
                         double *, int *, int *, const _Bool *);
typedef void (*pchprod) (int *, const int *, const int *, const _Bool *, const double *, const double *,
                         double *, double *);

#endif
