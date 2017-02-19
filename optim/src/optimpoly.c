/*
 ============================================================================
 Name        : optimpoly.c
 Author      : Okan
 Version     :
 Date        : 30/05/2016
 Description : Nonlinear optimization in C using the NLOPT library
 ============================================================================
 */

#include "SL.h"
#include "constants.h"
#include "optimpoly.h"
#include "utils.h"
#include "kinematics.h"
#include "stdlib.h"
#include "math.h"

/*
 * NLOPT optimization routine for table tennis traj gen
 *
 */
void nlopt_optim_poly_run(coptim params,
						  cracket racket) {

	static double x[OPTIM_DIM];
	static double tol[EQ_CONSTR_DIM];
	const_vec(EQ_CONSTR_DIM,1e-2,tol);
	init_soln_to_rest_posture(params,x); //parameters are the initial joint positions q0
	// set tolerances equal to second argument //

	nlopt_opt opt;
	opt = nlopt_create(NLOPT_LN_COBYLA, OPTIM_DIM);
	// LN = does not require gradients //
	nlopt_set_lower_bounds(opt, params.lb);
	nlopt_set_upper_bounds(opt, params.ub);
	nlopt_set_min_objective(opt, costfunc, &params);
	nlopt_add_inequality_mconstraint(opt, INEQ_CONSTR_DIM, joint_limits_ineq_constr,
			                         &params, tol);
	nlopt_add_equality_mconstraint(opt, EQ_CONSTR_DIM, kinematics_eq_constr,
			                         &racket, tol);
	nlopt_set_xtol_rel(opt, 1e-2);

	double initTime = get_time();
	double minf; // the minimum objective value, upon return //
	int res; // error code

	if ((res = nlopt_optimize(opt, x, &minf)) < 0) {
	    printf("NLOPT failed with exit code %d!\n", res);
	    //test_optim(x,TRUE);
	}
	else {
		//nlopt_example_run();
		printf("NLOPT success with exit code %d!\n", res);
		printf("NLOPT took %f ms\n", (get_time() - initTime)/1e3);
	    printf("Found minimum at f = %0.10g\n", minf);
	    //test_optim(x,TRUE);
	}
	nlopt_destroy(opt);
}

/*
 * Debug by testing the constraint violation of the solution vector
 *
 */
double test_optim(double *x, int info) {

	// give info on constraint violation
	double *grad = FALSE;
	static double kin_violation[EQ_CONSTR_DIM];
	static double lim_violation[INEQ_CONSTR_DIM]; // joint limit violations on strike and return
	kinematics_eq_constr(EQ_CONSTR_DIM, kin_violation, OPTIM_DIM, x, grad, NULL);
	joint_limits_ineq_constr(INEQ_CONSTR_DIM, lim_violation,
			                 OPTIM_DIM, x, grad, NULL);
	double cost = costfunc(OPTIM_DIM, x, grad, NULL);

	if (info) {
		// give info on solution vector
		print_optim_vec(x);
		printf("f = %.2f\n",cost);
		printf("Position constraint violation: [%.2f %.2f %.2f]\n",kin_violation[0],kin_violation[1],kin_violation[2]);
		printf("Velocity constraint violation: [%.2f %.2f %.2f]\n",kin_violation[3],kin_violation[4],kin_violation[5]);
		printf("Normal constraint violation: [%.2f %.2f %.2f]\n",kin_violation[6],kin_violation[7],kin_violation[8]);
	}

	for (int i = 0; i < INEQ_CONSTR_DIM; i++) {
		if (lim_violation[i] > 0.0)
			printf("Joint limit violated by %.2f on joint %d\n", lim_violation[i], i % DOF + 1);
	}

	return fmax(max_abs_array(kin_violation,EQ_CONSTR_DIM),
			    max_array(lim_violation,INEQ_CONSTR_DIM));
}

/*
 * Calculates the cost function for table tennis trajectory generation optimization
 * to find spline (3rd order strike+return) polynomials
 */
static double const_costfunc(unsigned n, const double *x, double *grad, void *my_func_params) {

	return 1.0;
}

/*
 * Calculates the cost function for table tennis trajectory generation optimization
 * to find spline (3rd order strike+return) polynomials
 */
static double costfunc(unsigned n, const double *x, double *grad, void *my_func_params) {

	double a1[DOF];
	double a2[DOF];
	static double *q0dot; // initial joint velocity
	static double *q0; // initial joint pos
	static int firsttime = TRUE;
	double T = x[2*DOF];

	if (firsttime) {
		firsttime = FALSE;
		coptim* optim_data = (coptim*)my_func_params;
		q0 = optim_data->q0;
		q0dot = optim_data->q0dot;
	}

	// calculate the polynomial coeffs which are used in the cost calculation
	calc_strike_poly_coeff(q0,q0dot,x,a1,a2);

	return T * (3*T*T*inner_prod(a1,a1) +
			3*T*inner_prod(a1,a2) + inner_prod(a2,a2));
}

/*
 * This is the inequality constraint that makes sure we never exceed the
 * joint limits during the striking and returning motion
 *
 */
static void joint_limits_ineq_constr(unsigned m, double *result,
		unsigned n, const double *x, double *grad, void *my_func_params) {

	static double a1[DOF];
	static double a2[DOF];
	static double a1ret[DOF]; // coefficients for the returning polynomials
	static double a2ret[DOF];
	static double joint_strike_max_cand[DOF];
	static double joint_strike_min_cand[DOF];
	static double joint_return_max_cand[DOF];
	static double joint_return_min_cand[DOF];
	static double *q0;
	static double *q0dot;
	static double *ub;
	static double *lb;
	static double Tret;
	static int firsttime = TRUE;

	if (firsttime) {
		firsttime = FALSE;
		coptim* optim_data = (coptim*)my_func_params;
		q0 = optim_data->q0;
		q0dot = optim_data->q0dot;
		ub = optim_data->ub;
		lb = optim_data->lb;
		Tret = optim_data->time2return;
	}

	// calculate the polynomial coeffs which are used for checking joint limits
	calc_strike_poly_coeff(q0,q0dot,x,a1,a2);
	calc_return_poly_coeff(q0,q0dot,x,Tret,a1ret,a2ret);
	// calculate the candidate extrema both for strike and return
	calc_strike_extrema_cand(a1,a2,x[2*DOF],q0,q0dot,
			joint_strike_max_cand,joint_strike_min_cand);
	calc_return_extrema_cand(a1ret,a2ret,x,Tret,joint_return_max_cand,joint_return_min_cand);

	/* deviations from joint min and max */
	for (int i = 0; i < DOF; i++) {
		result[i] = joint_strike_max_cand[i] - ub[i];
		result[i+DOF] = lb[i] - joint_strike_min_cand[i];
		result[i+2*DOF] = joint_return_max_cand[i] - ub[i];
		result[i+3*DOF] = lb[i] - joint_return_min_cand[i];
		//printf("%f %f %f %f\n", result[i],result[i+DOF],result[i+2*DOF],result[i+3*DOF]);
	}

}

/*
 * This is the constraint that makes sure we hit the ball
 */
static void kinematics_eq_constr(unsigned m, double *result, unsigned n,
		                  const double *x, double *grad, void *my_function_data) {

	static double racket_des_pos[CART];
	static double racket_des_vel[CART];
	static double racket_des_normal[CART];
	static Matrix link_pos_des;
	static Matrix joint_origin_pos_des;
	static Matrix joint_axis_pos_des;
	static Matrix Alink_des[N_LINKS+1];
	static Matrix racketTransform;
	static Matrix Jacobi;
	static Vector qfdot;
	static Vector xfdot;
	static double normal[CART];
	static int firsttime = TRUE;

	static double q[DOF];
	static double base_orient[4]; // quat
	static double base_state[3]; // pos
	static double eff_angles[3]; // initial euler angles for racket
	static double eff_pos[3]; // initial racket positions
	static cracket* racket;

	double T = x[2*DOF];
	int i;

	/* initialization of static variables */
	if (firsttime) {
		firsttime = FALSE;

		link_pos_des         = my_matrix(1,N_LINKS,1,3);
		joint_origin_pos_des = my_matrix(1,N_DOFS,1,3);
		joint_axis_pos_des   = my_matrix(1,N_DOFS,1,3);
		Jacobi               = my_matrix(1,2*CART,1,N_DOFS);
		racketTransform      = my_matrix(1,4,1,4);
		qfdot                = my_vector(1,DOF);
		xfdot                = my_vector(1,2*CART);

		for (i = 0; i <= N_LINKS; ++i)
			Alink_des[i] = my_matrix(1,4,1,4);

		base_orient[1] = 1.0;
		eff_pos[Z] = 0.3; // set racket

		// homogeneous transform instead of using quaternions
		racketTransform[1][1] = 1;
		racketTransform[2][3] = 1;
		racketTransform[3][2] = -1;
		racketTransform[4][4] = 1;

		racket = (cracket*) my_function_data;
	}

	if (grad) {
		// compute gradient of kinematics = jacobian
		//TODO:
		grad[0] = 0.0;
		grad[1] = 0.0;
		grad[2] = 0.0;
	}

	// interpolate at time T to get the desired racket parameters
	first_order_hold(racket,T,racket_des_pos,racket_des_vel,racket_des_normal);

	// extract state information from array to joint_des_state structure
	for (i = 0; i < DOF; i++) {
		q[i] = x[i];
		qfdot[i+1] = x[i+DOF];
	}

	// compute the desired link positions
	kinematics(q, base_state, base_orient,
			      eff_angles, eff_pos,
			      joint_axis_pos_des, joint_origin_pos_des,
			      link_pos_des, Alink_des);

	/*
	// compute the racket normal
	/*mat_mult(Alink_des[6],racketTransform,Alink_des[6]);
	for (i = 0; i < CART; i++) {
		normal[i] = Alink_des[6][i+1][3];
	}

	// compute the jacobian //
	jacobian(link_pos_des, joint_origin_pos_des, joint_axis_pos_des, Jacobi);
	mat_vec_mult(Jacobi, qfdot, xfdot);

	// deviations from the desired racket frame
	for (i = 0; i < CART; i++) {
		result[i] = link_pos_des[6][i+1] - racket_des_pos[i];
		result[i + CART] = xfdot[i+1] - racket_des_vel[i];
		result[i + 2*CART] = normal[i] - racket_des_normal[i];
	}*/

}

/*
 * First order hold to interpolate linearly at time T
 * between racket pos,vel,normal entries
 *
 * IF T is nan, racket variables are assigned to zero-element of
 * relevant racket entries
 *
 */
static void first_order_hold(const cracket* racket, const double T, double racket_pos[NCART],
		               double racket_vel[NCART], double racket_n[NCART]) {

	if (isnan(T)) {
		printf("Warning: T value is nan!\n");

		for(int i = 0; i < CART; i++) {
			racket_pos[i] = racket->pos[i][0];
			racket_vel[i] = racket->vel[i][0];
			racket_n[i] = racket->normal[i][0];
		}
	}
	else {
		int N = (int) (T/dt);
		double Tdiff = T - N*dt;
		int Nmax = racket->Nmax;

		for (int i = 0; i < CART; i++) {
			if (N < Nmax) {
				racket_pos[i] = racket->pos[i][N] +
						(Tdiff/dt) * (racket->pos[i][N+1] - racket->pos[i][N]);
				racket_vel[i] = racket->vel[i][N] +
						(Tdiff/dt) * (racket->vel[i][N+1] - racket->vel[i][N]);
				racket_n[i] = racket->normal[i][N] +
						(Tdiff/dt) * (racket->normal[i][N+1] - racket->normal[i][N]);
			}
			else {
				racket_pos[i] = racket->pos[i][N];
				racket_vel[i] = racket->vel[i][N];
				racket_n[i] = racket->normal[i][N];
			}
		}
	}
}

/*
 * Calculate the polynomial coefficients from the optimized variables qf,qfdot,T
 * p(t) = a1*t^3 + a2*t^2 + a3*t + a4
 */
static void calc_strike_poly_coeff(const double *q0, const double *q0dot, const double *x,
		                    double *a1, double *a2) {

	double T = x[2*DOF];

	for (int i = 0; i < DOF; i++) {
		a1[i] = (2/pow(T,3))*(q0[i]-x[i]) + (1/(T*T))*(q0dot[i] + x[i+DOF]);
		a2[i] = (3/(T*T))*(x[i]-q0[i]) - (1/T)*(x[i+DOF] + 2*q0dot[i]);
	}

	return;
}

/*
 * Calculate the returning polynomial coefficients from the optimized variables qf,qfdot
 * and time to return constant T
 * p(t) = a1*t^3 + a2*t^2 + a3*t + a4
 */
static void calc_return_poly_coeff(const double *q0, const double *q0dot,
		                    const double *x, const double T,
		                    double *a1, double *a2) {

	for (int i = 0; i < DOF; i++) {
		a1[i] = (2/pow(T,3))*(x[i]-q0[i]) + (1/(T*T))*(q0dot[i] + x[i+DOF]);
		a2[i] = (3/(T*T))*(q0[i]-x[i]) - (1/T)*(2*x[i+DOF] + q0dot[i]);
	}

}

/*
 * Calculate the extrema candidates for each joint (2*7 candidates in total)
 * For the striking polynomial
 * Clamp to [0,T]
 *
 */
static void calc_strike_extrema_cand(const double *a1, const double *a2, const double T,
		                      const double *q0, const double *q0dot,
		                      double *joint_max_cand, double *joint_min_cand) {

	static double cand1, cand2;

	for (int i = 0; i < DOF; i++) {
		cand1 = fmin(T,fmax(0,(-a2[i] + sqrt(a2[i]*a2[i] - 3*a1[i]*q0dot[i]))/(3*a1[i])));
		cand2 =  fmin(T,fmax(0,(-a2[i] - sqrt(a2[i]*a2[i] - 3*a1[i]*q0dot[i]))/(3*a1[i])));
		cand1 = a1[i]*pow(cand1,3) + a2[i]*pow(cand1,2) + q0dot[i]*cand1 + q0[i];
		cand2 = a1[i]*pow(cand2,3) + a2[i]*pow(cand2,2) + q0dot[i]*cand2 + q0[i];
		joint_max_cand[i] = fmax(cand1,cand2);
		joint_min_cand[i] = fmin(cand1,cand2);
	}
}

/*
 * Calculate the extrema candidates for each joint (2*7 candidates in total)
 * For the return polynomial
 * Clamp to [0,TIME2RETURN]
 *
 */
static void calc_return_extrema_cand(const double *a1, const double *a2,
		                      const double *x, const double Tret,
		                      double *joint_max_cand, double *joint_min_cand) {

	static double cand1, cand2;

	for (int i = 0; i < DOF; i++) {
		cand1 = fmin(Tret, fmax(0,(-a2[i] + sqrt(a2[i]*a2[i] - 3*a1[i]*x[i+DOF]))/(3*a1[i])));
		cand2 =  fmin(Tret, fmax(0,(-a2[i] - sqrt(a2[i]*a2[i] - 3*a1[i]*x[i+DOF]))/(3*a1[i])));
		cand1 = a1[i]*pow(cand1,3) + a2[i]*pow(cand1,2) + x[i+DOF]*cand1 + x[i];
		cand2 = a1[i]*pow(cand2,3) + a2[i]*pow(cand2,2) + x[i+DOF]*cand2 + x[i];
		joint_max_cand[i] = fmax(cand1,cand2);
		joint_min_cand[i] = fmin(cand1,cand2);
	}
}

/*
 * Estimate an initial solution for NLOPT
 * 2*dof + 1 dimensional problem
 *
 * The closer to the optimum it is the faster alg should converge
 */
void init_soln_to_rest_posture(coptim params, double x[OPTIM_DIM]) {

	// initialize first dof entries to q0
	for (int i = 0; i < NDOF; i++) {
		x[i] = params.qrest[i];
		x[i+NDOF] = 0.0;
	}
	x[2*NDOF] = 0.50;
}
