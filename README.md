NOPE - Nonlinear Optimization Preprocessing tEchniques
==================

**Warning**
This library tries to bring preprocessing to nonlinear programming problems.
Everything here is experimental, even the name, use at your own risk.

In this library, we try to implement a preprocessor to use before calling a
nonlinear optimization method. Since much varies within methods, some guidelines
are necessary:

  - We use [CUTEst](http://ccpforge.cse.rl.ac.uk/gf/project/cutest/wiki) as the
    matching syntax;
  - We want compatibility with Fortran (because other people may want to use
    this);
  - We want this to work with our method
    ([DCICPP](https://github.com/abelsiqueira/dcicpp.git)).

Note that we have not tested this in any way with Fortran, so we are currently
not following our own guidelines. If you want to help here, we welcome you.

LICENSE
-------

This work is licensed under GPLv3. See [LICENSE](LICENSE).

INSTALL
-------

Currently, only

    $ make

then link the library to your executable (in CUTEst, add in your prototype
file).
