# optimpoly

Optimal Trajectory Generation in Table Tennis

This repository investigates the nonlinear equality constrained optimization
necessary to generate table tennis striking (and returning) 3rd order
polynomials in joint space.

We use the NLOPT toolbox to run the constrained nonlinear optimization:
http://ab-initio.mit.edu/wiki/index.php/NLopt

Notes and Questions:

1. Should we have the repository completely independent from SL files (code and header)?
2. Stick to either 1 or 0-based indexing (but not both) ?
3. Include Eigen and switch to C++? 
