/*
 ============================================================================
 Name        : optimpoly.c
 Author      : Okan
 Version     :
 Date        : 30/05/2016
 Description : NLOPT example in C
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>

// optimization and math libraries
#include <math.h>
#include <nlopt.h>

// SL variables and kinematics
#include "SL.h"
#include "SL_user.h"
#include "SL_kinematics.h"

// defines
#define DOF 7
#define OPTIM_DIM 2*DOF+1

typedef struct {
    double a, b;
} my_constraint_data;

// example functions
double myfunc(unsigned n, const double *x, double *grad, void *my_func_data);
double myconstraint(unsigned n, const double *x, double *grad, void *data);
void nlopt_example_run();

// utility
long getTime();

// optimization methods
double costfunc(unsigned n, const double *x, double *grad, void *my_func_data);
void vec_minus(double *a1, const double *a2);
void vec_plus(double *a1, const double *a2);
double inner_prod(const double *a1, const double *a2);
void const_vec(const int n, const double val, double * vec);
void calc_poly_coeff(double *a1, double *a2, const double *q0, const double *x);

// global variables
static double q0[DOF];

int main(void) {

	// initialize the variables
	//q0 = [1.0; -0.2; -0.1; 1.8; -1.57; 0.1; 0.3];
	q0[0] = 1.0;
	q0[1] = -0.2;
	q0[2] = -0.1;
	q0[3] = 1.8;
	q0[4] = -1.57;
	q0[5] = 0.1;
	q0[6] = 0.3;

	double initTime = getTime();
	nlopt_example_run();
	printf("NLOPT took %f ms\n", (getTime() - initTime)/1e3);

	return TRUE;
}

/*
 * NLOPT optimization routine for table tennis traj gen
 */
void optim_poly_nlopt_run() {

	double val = 3.0;
	double lb[OPTIM_DIM]; /* lower bounds */
	double ub[OPTIM_DIM]; /* upper bounds */

	const_vec(OPTIM_DIM,-val,lb);
	const_vec(OPTIM_DIM,val,ub);

	nlopt_opt opt;
	opt = nlopt_create(NLOPT_LD_MMA, OPTIM_DIM); /* algorithm and dimensionality */
	nlopt_set_lower_bounds(opt, lb);
	nlopt_set_upper_bounds(opt, ub);
	nlopt_set_min_objective(opt, costfunc, NULL);
	double tol = 1e-8;
	nlopt_set_xtol_rel(opt, 1e-4);

	double x[OPTIM_DIM] = {0};  /* some initial guess */
	double minf; /* the minimum objective value, upon return */

	double init_time = getTime();
	if (nlopt_optimize(opt, x, &minf) < 0) {
	    printf("NLOPT failed!\n");
	}
	else {
	    printf("Found minimum at f = %0.10g\n", minf);
	    printf("Optim took %f ms\n", (getTime() - init_time)/1e3);
	}
	nlopt_destroy(opt);
}

/*
 * Returns constant vector of val value from 1 to n
 */
void const_vec(const int n, const double val, double * vec) {

	int i;
	for (i = 0; i < n; i++) {
		vec[i] = val;
	}
}

/*
 * Calculates the cost function for table tennis trajectory generation optimization
 * to find spline (3rd order strike+return) polynomials
 */
double costfunc(unsigned n, const double *x, double *grad, void *my_func_data) {

	int i;
	double a1[DOF];
	double a2[DOF];
	double q0dot[DOF]; // all zeros
	double qfdot[DOF]; // all zeros
	double qf[DOF]; // opt value
	double T = x[2*DOF];

	if (grad) {
		//TODO:

		for (i = 0; i < DOF; i++) {
			qf[i] = x[i];
			qfdot[i] = x[i+DOF];
			grad[i] = (6/pow(T,3))*(qf[i]-q0[i]) - (3/(T*T))*(q0dot[i]+qfdot[i]);
			grad[i+DOF] = (-3/(T*T))*(qf[i]-q0[i]) + (1/T)*(2*qfdot[i]+q0dot[i]);
		}
		//time derivative of cost J
		vec_minus(qf,q0);
		//assuming q0dot = 0 here;
		grad[2*DOF] = (-9/pow(T,4))*inner_prod(qf,qf) +
					  (6/pow(T,3))*inner_prod(qf,qfdot) +
					  (3/(T*T))*inner_prod(q0dot,qfdot) -
					  (1/(T*T))*inner_prod(qfdot,qfdot);
	}



	// calculate the polynomial coeffs which are used in the cost calculation
	calc_poly_coeff(a1,a2,q0,x);

	return T * (3*T*T*inner_prod(a1,a1) +
			3*T*inner_prod(a1,a2) + inner_prod(a2,a2));
}

//
//% supply gradient if algorithm supports it
//if nargout > 1
//    G = [(6/(T^3))*(qf-q0) - (3/T^2)*(q0dot + qfdot);
//         (-3/T^2)*(qf-q0) + (1/T)*(2*qfdot + q0dot)];
//    G(end+1) = (-9/(T^4))*((qf-q0)'*(qf-q0)) + ...
//               (6/(T^3))*((qf-q0)'*(q0dot+qfdot)) +...
//               (3/(T^2))*(q0dot'*(q0dot + qfdot)) +...
//               -(1/T^2)*((qfdot+2*q0dot)'*(qfdot+2*q0dot));
//end


/*
 * Returns a1 - a2 vector into a1, assuming both have dof = 7 length
 */
void vec_minus(double *a1, const double *a2) {

	int i;
	for (i = 0; i < DOF; i++) {
		a1[i] = a1[i] - a2[i];
	}
}

/*
 * Returns a1 + a2 vector into a1, assuming both have dof = 7 length
 */
void vec_plus(double *a1, const double *a2) {

	int i;
	for (i = 0; i < DOF; i++) {
		a1[i] = a1[i] + a2[i];
	}
}

/*
 * Returns the inner product between two vectors of size N_DOFS
 */
double inner_prod(const double *a1, const double *a2) {

	int i;
	double val = 0.0;
	for (i = 0; i < DOF; i++) {
		val += a1[i]*a2[i];
	}

	return val;
}

/*
 * Calculate the polynomial coefficients from the optimized variables qf,qfdot,T
 * p(t) = a1*t^3 + a2*t^2 + a3*t + a4
 */
void calc_poly_coeff(double *a1, double *a2, const double *q0, const double *x) {

	// variables to be optimized
	static double q0dot[DOF];
	int i;
	double T = x[2*DOF];

	for (i = 0; i < DOF; i++) {
		a1[i] = (2/pow(T,3))*(q0[i]-x[i]) + (1/(T*T))*(q0dot[i] + x[i+DOF]);
		a2[i] = (3/(T*T))*(x[i]-q0[i]) - (1/T)*(x[i+DOF] + 2*q0dot[i]);
	}

	return;
}

/*
 * NLOPT nonlinear optimization example
 * http://ab-initio.mit.edu/wiki/index.php/NLopt_Tutorial
 *
 */
void nlopt_example_run() {

	double lb[2] = { -HUGE_VAL, 0 }; /* lower bounds */
	nlopt_opt opt;
	opt = nlopt_create(NLOPT_LD_MMA, 2); /* algorithm and dimensionality */
	nlopt_set_lower_bounds(opt, lb);
	nlopt_set_min_objective(opt, myfunc, NULL);
	my_constraint_data data[2] = { {2,0}, {-1,1} };
	double tol = 1e-8;
	nlopt_add_inequality_constraint(opt, myconstraint, &data[0], tol);
	nlopt_add_inequality_constraint(opt, myconstraint, &data[1], tol);
	nlopt_set_xtol_rel(opt, 1e-4);

	double x[2] = { 1.234, 5.678 };  /* some initial guess */
	double minf; /* the minimum objective value, upon return */

	if (nlopt_optimize(opt, x, &minf) < 0) {
	    printf("NLOPT failed!\n");
	}
	else {
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
long getTime() {
	struct timeval tv;
	if (gettimeofday(&tv, (struct timezone *)0) == 0)
		return (tv.tv_sec*1000*1000 + tv.tv_usec);  //us

	return 0.;
}
