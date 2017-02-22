# polyoptim

Optimal Trajectory Generation in Table Tennis

This repository investigates the nonlinear equality constrained optimization
necessary to generate table tennis striking (and returning) 3rd order
polynomials in joint space.

We use the NLOPT toolbox to run the constrained nonlinear optimization:
http://ab-initio.mit.edu/wiki/index.php/NLopt

We investigate the following problems:

1. Testing for faster optimization - can we speed up the optimization with more optimized structures and gradient-based optimizers?

2. Can we use automatic differentiation to supply gradients ?

3. Can we try altogether different optimal control strategies in table tennis? Can we use NLOPT to solve them sufficiently well?

4. Should we try different optimization toolboxes?

# carma

Includes ARMADILLO matrix library in C++. 
Make sure to run tcsh and then check the following:

```echo $LD_LIBRARY_PATH```

If carma directory is not included add it, e.g. :

```setenv LD_LIBRARY_PATH ~/robolab/barrett/carma```

