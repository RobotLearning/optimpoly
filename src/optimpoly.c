/*
 ============================================================================
 Name        : optimpoly.c
 Author      : Okan
 Version     :
 Date        : 30/05/2016
 Description : Nonlinear optimization in C using the NLOPT library
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
#include "SL_user_common.h"
#include "SL_common.h"
#include "SL_kinematics_body.h"

// defines
#define DOF 7
#define OPTIM_DIM 2*DOF+1
#define MAX_VEL 200

// global variables
char joint_names[][20]= {
  {"dummy"},
  {"R_SFE"},
  {"R_SAA"},
  {"R_HR"},
  {"R_EB"},
  {"R_WR"},
  {"R_WFE"},
  {"R_WAA"}
};

// initialization needs to be done for this mapping
int  link2endeffmap[] = {0,PALM};

static double q0[DOF];

// utility method
long get_time();
void vec_minus(double *a1, const double *a2);
void vec_plus(double *a1, const double *a2);
double inner_prod(const double *a1, const double *a2);
void const_vec(const int n, const double val, double * vec);

// optimization related methods
double costfunc(unsigned n, const double *x, double *grad, void *my_func_data);
void kinematics_constr(unsigned m, double *result, unsigned n, const double *x, double *grad, void *f_data);
void calc_poly_coeff(double *a1, double *a2, const double *q0, const double *x);
void optim_poly_nlopt_run();
void guesstimate_soln(double * x);
void set_bounds(double *lb, double *ub);
void load_joint_limits();
int read_sensor_offsets(char *fname);


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

	double initTime = get_time();
	//nlopt_example_run();
	optim_poly_nlopt_run();
	printf("NLOPT took %f ms\n", (get_time() - initTime)/1e3);

	return TRUE;
}

/*
 * NLOPT optimization routine for table tennis traj gen
 */
void optim_poly_nlopt_run() {

	double val = 3.0;
	static double tol[OPTIM_DIM];
	static double lb[OPTIM_DIM]; /* lower bounds */
	static double ub[OPTIM_DIM]; /* upper bounds */
	static double x[OPTIM_DIM]; /* initial guess */

	load_joint_limits();
	set_bounds(lb,ub);
	const_vec(OPTIM_DIM,1e-2,tol);

	nlopt_opt opt;
	opt = nlopt_create(NLOPT_LD_SLSQP, OPTIM_DIM); /* algorithm and dimensionality */
	nlopt_set_lower_bounds(opt, lb);
	nlopt_set_upper_bounds(opt, ub);
	nlopt_set_min_objective(opt, costfunc, NULL);
	nlopt_add_equality_mconstraint(opt, 3, kinematics_constr, NULL, tol);

	nlopt_set_xtol_rel(opt, 1e-2);

	//double maxtime = 1e-3;
	//nlopt_set_maxtime(opt, maxtime);

	guesstimate_soln(x);

	double minf; /* the minimum objective value, upon return */

	if (nlopt_optimize(opt, x, &minf) < 0) {
	    printf("NLOPT failed!\n");
	}
	else {
	    printf("Found minimum at f = %0.10g\n", minf);
	}
	nlopt_destroy(opt);
}

/*
 * Set upper and lower bounds on the optimization
 */
void set_bounds(double *lb, double *ub) {

	// lower bounds and upper bounds for qf are the joint limits
	int i;
	for (i = 1; i <= DOF; i++) {
		ub[i-1] = joint_range[i][MAX_THETA];
		lb[i-1] = joint_range[i][MIN_THETA];
		ub[i-1+DOF] = MAX_VEL;
		lb[i-1+DOF] = -MAX_VEL;
	}
	// constraints on final time
	ub[2*DOF] = 1.0;
	lb[2*DOF] = 0.0;
}

/*
 * Load the joint limits from config/ file into joint_range array
 *
 */
void load_joint_limits() {

	char *fname = "SensorOffset.cf";
	read_sensor_offsets(fname);

}

/*
 * Copied from SL_common.h. Dont want to call the function from SL because
 * we have to include a lot of extra SL files
 *
 */
int read_sensor_offsets(char *fname) {

  int j,i,rc;
  char   string[100];
  FILE  *in;

  /* get the max, min, and offsets of the position sensors */

  sprintf(string,"%s/robolab/barrett/%s%s",getenv("HOME"),CONFIG,fname);
  in = fopen(string,"r");
  if (in == NULL) {
    printf("ERROR: Cannot open file >%s<!\n",string);
    return FALSE;
  }

  /* find all joint variables and read them into the appropriate array */

  for (i=1; i<= n_dofs; ++i) {
    if (!find_keyword(in, &(joint_names[i][0]))) {
      printf("ERROR: Cannot find offset for >%s<!\n",joint_names[i]);
      fclose(in);
      return FALSE;
    }
    rc=fscanf(in,"%lf %lf %lf %lf %lf %lf",
	&joint_range[i][MIN_THETA], &joint_range[i][MAX_THETA],
	   &(joint_default_state[i].th),
	   &(joint_opt_state[i].th),
	   &(joint_opt_state[i].w),
	   &joint_range[i][THETA_OFFSET]);
    joint_default_state[i].thd = 0;
    joint_default_state[i].uff = 0;
  }

  fclose(in);

  return TRUE;

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
 * Estimate an initial solution for NLOPT
 * 2*dof + 1 dimensional problem
 *
 * TODO: for now initializing to x = [q0,zeros,0.5]
 */
void guesstimate_soln(double * x) {

	// initialize first dof entries to q0
	int i;
	for (i = 0; i < DOF; i++) {
		x[i] = q0[i];
		x[i+DOF] = 0.0;
	}
	x[2*DOF] = 0.5;
}

/*
 * Calculates the cost function for table tennis trajectory generation optimization
 * to find spline (3rd order strike+return) polynomials
 */
double costfunc(unsigned n, const double *x, double *grad, void *my_func_data) {

	int i;
	double a1[DOF];
	double a2[DOF];
	static double q0dot[DOF]; // all zeros
	static double qfdot[DOF]; // all zeros
	static double qf[DOF]; // opt value
	double T = x[2*DOF];

	if (grad) {

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

/*
 * This is the constraint that makes sure we touch the ball
 */
void kinematics_constr(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data) {

	int i;

	double dt = 0.02;
	double **ballPred = (double **) data;
	double T = x[2*DOF];
	int N = T/dt;

	if (grad) {
		// compute gradient of kinematics = jacobian
		//TODO:
		grad[0] = 0.0;
		grad[1] = 0.0;
		grad[2] = 0.0;
	}



	/* compute the desired link positions */
	linkInformationDes(joint_des_state, &base_state, &base_orient,endeff,
			           joint_cog_mpos_des, joint_axis_pos_des, joint_origin_pos_des,
			           link_pos_des, Alink_des);

	/* the desired endeffector information */
	for (i = 1; i <= 3; ++i) {
		result[i-1] = link_pos_des[6][i] - ballPred[i-1][N];
	}

}


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
 * Return time of day as micro seconds
 */
long get_time() {
	struct timeval tv;
	if (gettimeofday(&tv, (struct timezone *)0) == 0)
		return (tv.tv_sec*1000*1000 + tv.tv_usec);  //us

	return 0.;
}
