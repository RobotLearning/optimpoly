/*
 * example.c
 *
 * Examples testing various libraries:
 *
 * 1. NLOPT Example showing nonlinear optimization with derivatives and nonlinear constraints.
 * One can run the following example: http://ab-initio.mit.edu/wiki/index.php/NLopt_Tutorial
 * by including nlopt_example_run() in main()
 *
 * Here we can test various optimization algorithms provided by NLOPT
 * Especially LN (local, no-derivative) vs. LD (local, derivative-supplied)
 *
 * Created on: Jun 2, 2016
 *      Author: okoc
 */

#include <stdio.h>
#include <stdlib.h>
#include "sys/time.h"

// optimization and math libraries
#include <math.h>
#include <nlopt.h>

// example functions
double myfunc(unsigned n, const double *x, double *grad, void *my_func_data);
double myconstraint(unsigned n, const double *x, double *grad, void *data);
void nlopt_example_run();
long get_time();

int main() {
	nlopt_example_run();
	return 1;
}

/*
 * NLOPT nonlinear optimization example
 * http://ab-initio.mit.edu/wiki/index.php/NLopt_Tutorial
 *
 */
void nlopt_example_run() {

	double x[3] = { 1.234, 5.678, 2.456 };  /* some initial guess */
	double minf; /* the minimum objective value, upon return */
	double tol = 1e-4;
	double lb[3] = { -HUGE_VAL, -10, -20 }; /* lower bounds */
	double init_time;
	nlopt_opt opt;

	opt = nlopt_create(NLOPT_LD_MMA, 3);
	//opt = nlopt_create(NLOPT_LN_COBYLA, 3); /* algorithm and dimensionality */
	nlopt_set_lower_bounds(opt, lb);
	nlopt_set_min_objective(opt, myfunc, NULL);
	nlopt_add_inequality_constraint(opt, myconstraint, NULL, tol);
	nlopt_set_xtol_rel(opt, 1e-2);

	init_time = get_time();

	if (nlopt_optimize(opt, x, &minf) < 0) {
	    printf("NLOPT failed!\n");
	}
	else {
		printf("NLOPT took %f ms\n", (get_time() - init_time)/1e3);
	    printf("Found minimum at f(%g,%g,%g) = %0.10g\n", x[0], x[1], x[2], minf);
	}
	nlopt_destroy(opt);

}

/*
 * Cost function
 * Calculates also the gradient if grad is TRUE
 *
 */
double myfunc(unsigned n, const double *x, double *grad, void *my_func_data) {

    if (grad) {
        grad[0] = 2*x[0];
        grad[1] = 2*x[1];
        grad[2] = 2*x[2];
    }
    return (x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
}

/*
 * Constraint function that returns the amount of constraint violation
 * Calculates also the gradient of the constraint function if grad is TRUE
 */
double myconstraint(unsigned n, const double *x, double *grad, void *data) {

    if (grad) {
        grad[0] = -x[1];
        grad[1] = -x[0];
        grad[2] = 0.0;
    }
    return 1 - x[0]*x[1];
}

/*
 * Return time of day as micro seconds
 */
long get_time() {
	struct timeval tv;
	if (gettimeofday(&tv, (struct timezone *)0) == 0)
		return (tv.tv_sec*1000*1000 + tv.tv_usec);  //us

	return 0.;
}


