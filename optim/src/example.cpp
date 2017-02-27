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

typedef struct {
    double a, b;
} my_constraint_data;

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

	double lb[2] = { -HUGE_VAL, 0.0 }; /* lower bounds */
	//nlopt_opt local_opt;
	//local_opt = nlopt_create(NLOPT_LD_TNEWTON_PRECOND, 2);
	//nlopt_set_xtol_rel(local_opt, 1e-4);
	nlopt_opt opt;
	opt = nlopt_create(NLOPT_LN_COBYLA, 2); /* algorithm and dimensionality */
	//nlopt_set_local_optimizer(opt, local_opt);
	//nlopt_destroy(local_opt);
	nlopt_set_lower_bounds(opt, lb);
	nlopt_set_min_objective(opt, myfunc, NULL);
	my_constraint_data data[2] = { {2,0}, {-1,1} };
	double tol = 1e-8;
	nlopt_add_inequality_constraint(opt, myconstraint, &data[0], tol);
	nlopt_add_inequality_constraint(opt, myconstraint, &data[1], tol);
	nlopt_set_xtol_rel(opt, 1e-4);

	double x[2] = { 1.234, 5.678 };  /* some initial guess */
	double minf; /* the minimum objective value, upon return */
	double init_time = get_time();

	if (nlopt_optimize(opt, x, &minf) < 0) {
	    printf("NLOPT failed!\n");
	}
	else {
		printf("NLOPT took %f ms\n", (get_time() - init_time)/1e3);
	    printf("Found minimum at f(%g,%g) = %0.10g\n", x[0], x[1], minf);
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
        grad[0] = 0.0;
        grad[1] = 0.5 / sqrt(x[1]);
    }
    return sqrt(x[1]);
}

/*
 * Constraint function that returns the amount of constraint violation
 * Calculates also the gradient of the constraint function if grad is TRUE
 */
double myconstraint(unsigned n, const double *x, double *grad, void *data) {

    my_constraint_data *d = (my_constraint_data *) data;
    double a = d->a, b = d->b;
    if (grad) {
        grad[0] = 3 * a * (a*x[0] + b) * (a*x[0] + b);
        grad[1] = -1.0;
    }
    return ((a*x[0] + b) * (a*x[0] + b) * (a*x[0] + b) - x[1]);
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


