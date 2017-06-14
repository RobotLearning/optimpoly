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

#include <adolc/adouble.h>            // use of active doubles
#include <adolc/drivers/drivers.h>    // use of "Easy to Use" drivers
#include <adolc/taping.h>             // use of taping

// optimization and math libraries
#include <math.h>
#include <nlopt.h>

long get_time();

// example functions
template <class type>
type myfunc_exact(unsigned n, const type *x, double *grad, void *my_func_data);
template <class type>
type myfunc(unsigned n, const type *x, void *my_func_data);
template <class type>
type myconstraint(unsigned n, const type *x, double *grad, void *data);
void nlopt_example_run();

// auto diff versions
double myfunc_auto(unsigned n, const double *x, double *grad, void *my_func_data);
void nlopt_autodiff_run();


#define DIM 4
#define DIMY 2

typedef struct {
    adouble *x_auto;
    adouble y_auto;
} auto_diff_data;

int main() {
	printf("Testing toy problem with EXACT derivatives...\n");
	nlopt_example_run();
	printf("Testing toy problem with AUTO derivatives...\n");
	nlopt_autodiff_run();
	return 1;
}

/*
 * NLOPT + autodiff example
 */
void nlopt_autodiff_run() {

	auto_diff_data data = {new adouble[DIM], 1.0};
	double x[4] = { 1.234, 5.678, 2.456, -4.124 };  /* some initial guess */
	double minf; /* the minimum objective value, upon return */
	double tol = 1e-4;
	double lb[4] = { -HUGE_VAL, -10, -10, -10 }; /* lower bounds */
	double init_time;
	nlopt_opt opt;

	opt = nlopt_create(NLOPT_LD_MMA, 4);
	//opt = nlopt_create(NLOPT_LN_COBYLA, 4); /* algorithm and dimensionality */
	nlopt_set_lower_bounds(opt, lb);
	nlopt_set_min_objective(opt, myfunc_auto, &data);
	//nlopt_add_inequality_constraint(opt, myconstraint, NULL, tol);
	nlopt_set_xtol_rel(opt, 1e-2);

	init_time = get_time();

	if (nlopt_optimize(opt, x, &minf) < 0) {
	    printf("NLOPT failed!\n");
	}
	else {
		printf("NLOPT took %f ms\n", (get_time() - init_time)/1e3);
	    printf("Found minimum at f(%g,%g,%g,%g) = %0.10g\n", x[0], x[1], x[2], x[3], minf);
	}
	nlopt_destroy(opt);
	delete[] data.x_auto;
}

/*
 * NLOPT nonlinear optimization example
 */
void nlopt_example_run() {

	double x[4] = { 1.234, 5.678, 2.456, -4.124 };  /* some initial guess */
	double minf; /* the minimum objective value, upon return */
	double tol = 1e-4;
	double lb[4] = { -HUGE_VAL, -10, -10, -10 }; /* lower bounds */
	double init_time;
	nlopt_opt opt;

	opt = nlopt_create(NLOPT_LD_MMA, 4);
	//opt = nlopt_create(NLOPT_LN_COBYLA, 4); /* algorithm and dimensionality */
	nlopt_set_lower_bounds(opt, lb);
	nlopt_set_min_objective(opt, myfunc_exact, NULL);
	nlopt_add_inequality_constraint(opt, myconstraint, NULL, tol);
	nlopt_set_xtol_rel(opt, 1e-2);

	init_time = get_time();

	if (nlopt_optimize(opt, x, &minf) < 0) {
	    printf("NLOPT failed!\n");
	}
	else {
		printf("NLOPT took %f ms\n", (get_time() - init_time)/1e3);
	    printf("Found minimum at f(%g,%g,%g,%g) = %0.10g\n", x[0], x[1], x[2], x[3], minf);
	}
	nlopt_destroy(opt);

}

/*
 * Autodiffentiates the cost function
 * Provides value and also the auto-derivative values
 *
 */
double myfunc_auto(unsigned n, const double *x, double *grad, void *my_func_data) {

	if (grad) {
	    trace_on(1); // tag = 1, keep = 0 by default
	    auto_diff_data *data = (auto_diff_data*)my_func_data;
	    adouble *x_auto = data->x_auto;
	    adouble y_auto = 1;
	    double y;
	    for(int i = 0; i < DIM; i++) {
	        x_auto[i] <<= x[i];
	    }
	    y_auto = myfunc(n, x_auto, my_func_data);
	    y_auto >>= y;
	    trace_off(1);
	    gradient(1,DIM,x,grad);
	    return y;
	}
	else
		return myfunc(n, x, my_func_data);
}

/*
 * Cost function
 * Calculates also the gradient if grad is TRUE
 *
 */
template <class type>
type myfunc(unsigned n, const type *x, void *my_func_data) {

    return (x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[3]*x[3]);
}


/*
 * Cost function
 * Calculates also the gradient if grad is TRUE
 *
 */
template <class type>
type myfunc_exact(unsigned n, const type *x, double *grad, void *my_func_data) {

    if (grad) {
        grad[0] = 2*x[0];
        grad[1] = 2*x[1];
        grad[2] = 2*x[2];
        grad[3] = 2*x[3];
    }
    return (x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[3]*x[3]);
}

/*
 * Constraint function that returns the amount of constraint violation
 * Calculates also the gradient of the constraint function if grad is TRUE
 */
template <class type>
type myconstraint(unsigned n, const type *x, double *grad, void *data) {

    if (grad) {
        grad[0] = -x[1]*x[2]*x[3];
        grad[1] = -x[0]*x[2]*x[3];
        grad[2] = -x[0]*x[1]*x[3];
        grad[3] = -x[0]*x[1]*x[2];
    }
    return 1 - x[0]*x[1]*x[2]*x[3];
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


