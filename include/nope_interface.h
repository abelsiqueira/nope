#ifndef nope_interface_h
#define nope_interface_h

// Unconstrained
typedef void (*pusetup) (int *, int *, int *, int *, int *, double *, double *,
                         double *);
typedef void (*pufn)    (int *, int *, double *, double *);
typedef void (*puofg)   (int *, int *, double *, double *, double *, _Bool *);
typedef void (*puhprod) (int *, int *, _Bool *, double *, double *, double *);
// Constrained
typedef void (*pcdimen) (int *, int *, int *, int *);
typedef void (*pcdimsj) (int *, int *);
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
