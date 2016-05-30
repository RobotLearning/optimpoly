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

typedef struct {
    double a, b;
} my_constraint_data;

double myfunc(unsigned n, const double *x, double *grad, void *my_func_data);
double myconstraint(unsigned n, const double *x, double *grad, void *data);
void nlopt_example_run();
long getTime();

// global variables
int dof = 7;
double q0[dof];

int main(void) {

	SL_DJstate joint_des_state;

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
 * Calculates the cost function for table tennis trajectory generation optimization
 * to find spline (3rd order strike+return) polynomials
 */
double costfunc(unsigned n, const double *x, double *grad, void *my_func_data) {

	int i;
	const int dof = 7;
	double a1[dof];
	double a2[dof];
	double q0dot[dof]; // all zeros
	double qfdot[dof]; // all zeros
	double qf[dof]; // opt value
	double qfdot[dof];
	double T = x[2*dof];

	if (grad) {
		//TODO:

		for (i = 0; i < dof; i++) {
			qf[i] = x[i];
			qfdot[i] = x[i+dof];
			grad[i] = (6/pow(T,3))*(qf[i]-q0[i]) - (3/(T*T))*(q0dot[i]+qfdot[i]);
			grad[i+dof] = (-3/(T*T))*(qf[i]-q0[i]) + (1/T)*(2*qfdot[i]+q0dot[i]);
		}
		grad[2*dof] = (-9/pow(T,4))*inner_prod()//time derivative of cost J
	}



	// calculate the polynomial coeffs which are used in the cost calculation
	calc_poly_coeff(a1,a2,q0,x);

	return T * (3*T^2*inner_prod(a1,a1) + 3*T*inner_prod(a1,a2) + inner_prod(a2,a2));
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
 * Returns the inner product between two vectors of size N_DOFS
 */
double inner_prod(const double *a1, const double *a2) {

	int i;
	double val = 0.0;
	for (i = 0; i < dof; i++) {
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
	static double q0dot[dof];
	int i;
	double T = x[2*dof];

	for (i = 0; i < dof; i++) {
		a1[i] = (2/pow(T,3))*(q0[i]-x[i]) + (1/(T*T))*(q0dot[i] + x[i+dof]);
		a2[i] = (3/(T*T))*(x[i]-q0[i]) - (1/T)*(x[i+dof] + 2*q0dot[i]);
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
