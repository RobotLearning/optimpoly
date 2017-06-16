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

#define DIM 4 // dimension of the problem
#define CONSTR_DIM 2 // number of constraints
#define TAG_OBJ 1
#define TAG_CONSTR 2

typedef struct {
    double ** jac;
} auto_diff_data;

// template of objective that both auto-diff and exact versions share
template <class type>
type myfunc(unsigned n, const type *x, void *my_func_data);
template<class type>
void myconstraint(unsigned m, type *result, unsigned n, const type *x, void *data);

// exact objective and constraint functions
double myfunc_exact(unsigned n, const double *x, double *grad, void *my_func_data);
void myconstraint_exact(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data);
void nlopt_example_run();

// auto diff versions of objective and constraints
double myfunc_auto(unsigned n, const double *x, double *grad, void *my_func_data);
void myconstraint_auto(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data);
void nlopt_autodiff_run();
void generate_tape(auto_diff_data & data);
void get_starting_point(int n, double *x);
bool check_optim_result(const int res);
long get_time();

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

	double x[DIM];
	double minf; /* the minimum objective value, upon return */
	double tol_ineq[CONSTR_DIM] = {1e-4, 1e-4};
	double lb[4] = { -HUGE_VAL, -10, -10, -10 }; /* lower bounds */
	double init_time;
	int res; // error code
	nlopt_opt opt;
	auto_diff_data data;
	generate_tape(data);

	opt = nlopt_create(NLOPT_LD_MMA, 4);
	//opt = nlopt_create(NLOPT_LN_COBYLA, 4); /* algorithm and dimensionality */
	nlopt_set_lower_bounds(opt, lb);
	nlopt_set_min_objective(opt, myfunc_auto, &data);
	nlopt_add_inequality_mconstraint(opt, CONSTR_DIM, myconstraint_auto, &data, tol_ineq);
	nlopt_set_xtol_rel(opt, 1e-2);
	get_starting_point(DIM,x);

	init_time = get_time();
	if ((res = nlopt_optimize(opt, x, &minf)) < 0) {
	    printf("NLOPT failed!\n");
	}
	else {
		printf("NLOPT took %f ms\n", (get_time() - init_time)/1e3);
	    printf("Found minimum at f(%g,%g,%g,%g) = %0.10g\n", x[0], x[1], x[2], x[3], minf);
	}
	check_optim_result(res);
	nlopt_destroy(opt);
	for (int i = 0; i < CONSTR_DIM; i++)
	    delete[] data.jac[i];
	delete[] data.jac;
}

/*
 * NLOPT nonlinear optimization example
 */
void nlopt_example_run() {

	double x[4] = { 1.234, 5.678, 2.456, -4.124 };  /* some initial guess */
	double minf; /* the minimum objective value, upon return */
	double tol_ineq[CONSTR_DIM] = {1e-4, 1e-4};
	double lb[4] = { -HUGE_VAL, -10, -10, -10 }; /* lower bounds */
	double init_time;
	int res; // error code
	nlopt_opt opt;

	opt = nlopt_create(NLOPT_LD_MMA, 4);
	//opt = nlopt_create(NLOPT_LN_COBYLA, 4); /* algorithm and dimensionality */
	nlopt_set_lower_bounds(opt, lb);
	nlopt_set_min_objective(opt, myfunc_exact, NULL);
	nlopt_add_inequality_mconstraint(opt, CONSTR_DIM, myconstraint_exact, NULL, tol_ineq);
	nlopt_set_xtol_rel(opt, 1e-2);

	init_time = get_time();

	if ((res = nlopt_optimize(opt, x, &minf)) < 0) {
	    printf("NLOPT failed!\n");
	}
	else {
		printf("NLOPT took %f ms\n", (get_time() - init_time)/1e3);
	    printf("Found minimum at f(%g,%g,%g,%g) = %0.10g\n", x[0], x[1], x[2], x[3], minf);
	}
	check_optim_result(res);
	nlopt_destroy(opt);
}

/*
 * Cost function template
 *
 */
template <class type>
type myfunc(unsigned n, const type *x, void *my_func_data) {

    return (x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[3]*x[3]);
}

/*
 * Autodiffentiates the cost function
 * Provides value and also the auto-derivative values
 *
 */
double myfunc_auto(unsigned n, const double *x, double *grad, void *my_func_data) {

	if (grad) {
	    gradient(TAG_OBJ,DIM,x,grad);
	}
	return myfunc(n, x, my_func_data);
}

/*
 * Cost function
 * Calculates also the gradient if grad is TRUE
 *
 */
double myfunc_exact(unsigned n, const double *x, double *grad, void *my_func_data) {

    if (grad) {
        grad[0] = 2*x[0];
        grad[1] = 2*x[1];
        grad[2] = 2*x[2];
        grad[3] = 2*x[3];
    }
    return (x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[3]*x[3]);
}

/*
 * Constraint function template
 */
template<class type>
void myconstraint(unsigned m, type *result, unsigned n, const type *x, void *data) {

    result[0] = 1 - x[0]*x[1]*x[2]*x[3];
    result[1] = x[1]*x[2];
}

/*
 * Constraint function that returns the amount of constraint violation
 * Calculates also the gradient of the constraint function if grad is TRUE
 */
void myconstraint_auto(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data) {

    if (grad) {
    	auto_diff_data * ddata = (auto_diff_data*)data;
        jacobian(TAG_CONSTR,m,n,x,ddata->jac);
        int idx = 0;
        for(int i = 0; i < m; i++)
        	for(int j=0; j < n; j++)
        		grad[idx++] = ddata->jac[i][j];
    }
    result[0] = 1 - x[0]*x[1]*x[2]*x[3];
    result[1] = x[1]*x[2];
}

/*
 * Constraint function that returns the amount of constraint violation
 * Calculates also the gradient of the constraint function if grad is TRUE
 */
void myconstraint_exact(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data)  {

    if (grad) {
        grad[0] = -x[1]*x[2]*x[3];
        grad[1] = -x[0]*x[2]*x[3];
        grad[2] = -x[0]*x[1]*x[3];
        grad[3] = -x[0]*x[1]*x[2];
        grad[4] = 0.0;
        grad[5] = x[2];
        grad[6] = x[1];
        grad[7] = 0.0;
    }
    result[0] = 1 - x[0]*x[1]*x[2]*x[3];
    result[1] = x[1]*x[2];
}

void get_starting_point(int n, double *x) {

	double x0[4] = { 1.234, 5.678, 2.456, -4.124 };  /* some initial guess */
	for (int i = 0; i < n; i++)
		x[i] = x0[i];
}

/*
 * Initialize ADOLC operations
 */
void generate_tape(auto_diff_data & data) {

	int n = DIM;
	int m = CONSTR_DIM;

	adouble *x_auto = new adouble[n];
	adouble *g_auto = new adouble[m];
	adouble f_auto;
	double *x = new double[n];
	double sig;
	double dummy;
	data.jac = new double*[m];
	for (int i = 0; i < m; i++)
		data.jac[i] = new double[n];

	get_starting_point(n, x);

	trace_on(TAG_OBJ);
	for (int i = 0; i < n; i++)
		x_auto[i] <<= x[i];
	f_auto = myfunc(n,x_auto,NULL);
	f_auto >>= dummy;
	trace_off();

	trace_on(TAG_CONSTR);
	for (int i = 0; i < n; i++)
		x_auto[i] <<= x[i];

	myconstraint(m,g_auto,n,x_auto,NULL);
	for(int i = 0; i < m; i++)
		g_auto[i] >>= dummy;
	trace_off();

	delete[] x_auto;
	delete[] x;
	delete[] g_auto;
}

/*
 * Give info about the optimization after termination
 *
 */
bool check_optim_result(const int res) {

	bool flag = false;
	switch (res) {
	case NLOPT_SUCCESS:
		printf("Success!\n");
		flag = true;
		break;
	case NLOPT_STOPVAL_REACHED:
		printf("Optimization stopped because stopval (above) was reached.\n");
		flag = true;
		break;
	case NLOPT_FTOL_REACHED:
		printf("Optimization stopped because ftol_rel "
				"or ftol_abs (above) was reached.\n");
		flag = true;
		break;
	case NLOPT_XTOL_REACHED:
		flag = true;
		printf("Optimization stopped because xtol_rel or xtol_abs (above) was reached.\n");
		break;
	case NLOPT_MAXEVAL_REACHED:
		flag = true;
		printf("Optimization stopped because maxeval (above) was reached.\n");
		break;
	case NLOPT_MAXTIME_REACHED:
		flag = true;
		printf("Optimization stopped because maxtime (above) was reached.\n");
		break;
	case NLOPT_FAILURE:
		printf("Epic fail!\n");
		break;
	case NLOPT_INVALID_ARGS:
		printf("Invalid arguments (e.g. lower bounds are bigger than "
				"upper bounds, an unknown algorithm was specified, etcetera).\n");
		break;
	case NLOPT_OUT_OF_MEMORY:
		printf("Ran out of memory!\n");
		break;
	case NLOPT_ROUNDOFF_LIMITED:
		printf("Halted because roundoff errors limited progress."
			"(In this case, the optimization still typically returns a useful result.\n");
		break;
	case NLOPT_FORCED_STOP:
		printf("Halted because of a forced termination: "
				"the user called nlopt_force_stop(opt)"
				"on the optimization’s nlopt_opt object "
				"opt from the user’s objective function or constraints.\n");
		break;

	}
	return flag;
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
