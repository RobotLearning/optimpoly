/**
 * @file invkin.cpp
 *
 * @brief Inverse Kinematics optimization to find qf, qfdot joint angles
 * and velocities at Virtual Hitting Plane (VHP).
 *
 *  Created on: Mar 5, 2017
 *      Author: okoc
 */

#include "constants.h"
#include "utils.h"
#include "stdlib.h"
#include "math.h"
#include "kinematics.h"
#include "optim.h"

// firsttime checking
static bool firsttime[2]; //! TODO: remove this global var!

static double penalize_dist_to_limits(unsigned n, const double *x,
		                     double *grad, void *my_func_params);
static double const_costfunc(unsigned n, const double *x,
		                     double *grad, void *my_func_params);
static void init_last_soln(const optim * params, double x[2*NDOF]);
static void init_rest_soln(const coptim * params, double x[2*NDOF]);
static void kinematics_eq_constr(unsigned m, double *result, unsigned n,
		                  const double *x, double *grad, void *my_function_data);
static double test_optim(const double *x, const double T,
		          coptim * coparams,
				  racketdes * racketdata);
static void joint_limits_ineq_constr(unsigned m, double *result,
		                      unsigned n, const double *x, double *grad, void *data);

static void calc_strike_poly_coeff(const double *q0, const double *q0dot, const double *x,
		                    double *a1, double *a2);
static void calc_return_poly_coeff(const double *q0, const double *q0dot,
		                    const double *x, const double time2return,
		                    double *a1, double *a2);
static void calc_strike_extrema_cand(const double *a1, const double *a2, const double T,
		                      const double *q0, const double *q0dot,
							  double *joint_max_cand, double *joint_min_cand);
static void calc_return_extrema_cand(const double *a1, const double *a2,
		                      const double *x, const double time2return,
							  double *joint_max_cand, double *joint_min_cand);
static void finalize_soln(const double* x, optim * params,
		                  bool detach, double time_elapsed);

void VHP::run() {

	predict();
	optim();
}

void VHP::optim() {

}

/**
 * @brief inverse kinematics for table tennis trajectory generation
 *        based on a fixed Virtual Hitting Plane (VHP).
 *
 * Multi-threading entry point for the NLOPT optimization.
 * The optimization problem is solved online using COBYLA (see NLOPT).
 * Cost function in this case is the sum of squared distances from
 * joint limits.
 * VHP is held fixed at a constant (see constants.h) y-location.
 *
 * @param coparams Co-optimization parameters held fixed during optimization.
 * @param racketdata Predicted racket position,velocity and normals.
 * @param params Optimization parameters updated if the solution found is FEASIBLE.
 * @return the maximum of violations
 *         Maximum of :
 * 		   1. kinematics equality constraint violations
 *         2. joint limit violations throughout trajectory
 */
double nlopt_vhp_run(coptim *coparams,
					 racketdes *racketdata,
					 optim *params) {

	firsttime[0] = true;
	firsttime[1] = true;

	//print_input_structs(coparams, racketdata, params);
	params->update = false;
	params->running = true;
	nlopt_opt opt;
	double x[2*NDOF];
	double tol_eq[EQ_CONSTR_DIM];
	const_vec(EQ_CONSTR_DIM,1e-2,tol_eq);

	if (coparams->moving)
		init_last_soln(params,x);
	else
		init_rest_soln(coparams,x);

	// LN = does not require gradients //
	opt = nlopt_create(NLOPT_LN_COBYLA, 2*NDOF);
	nlopt_set_xtol_rel(opt, 1e-2);
	nlopt_set_lower_bounds(opt, coparams->lb);
	nlopt_set_upper_bounds(opt, coparams->ub);
	nlopt_set_min_objective(opt, penalize_dist_to_limits, coparams);
	nlopt_add_equality_mconstraint(opt, EQ_CONSTR_DIM, kinematics_eq_constr,
			                         racketdata, tol_eq);

	double init_time = get_time();
	double past_time = 0.0;
	double minf; // the minimum objective value, upon return //
	int res; // error code
	double max_violation;

	if ((res = nlopt_optimize(opt, x, &minf)) < 0) {
		if (coparams->verbose)
			printf("NLOPT failed with exit code %d!\n", res);
	    past_time = (get_time() - init_time)/1e3;
	    max_violation = 100.0;
	}
	else {
		past_time = (get_time() - init_time)/1e3;
		if (coparams->verbose) {
			printf("NLOPT success with exit code %d!\n", res);
			printf("NLOPT took %f ms\n", past_time);
			printf("Found minimum at f = %0.10g\n", minf);
		}
	    max_violation = test_optim(x,params->T,coparams,racketdata);
	    if (max_violation < 1e-2)
	    	finalize_soln(x,params,coparams->detach,past_time);
	}
	params->running = false;
	//nlopt_destroy(opt);
	return max_violation;
}

/*
 * Penalize (unweighted) squared distance (in joint space) to joint limit averages
 *
 * Adds also joint velocity penalties
 *
 */
static double penalize_dist_to_limits(unsigned n, const double *x, double *grad, void *my_func_params) {

	static double *ub;
	static double *lb;
	static double limit_avg[NDOF];
	double cost = 0.0;

	if (firsttime[0]) {
		firsttime[0] = false;
		coptim* optim_data = (coptim*)my_func_params;
		ub = optim_data->ub;
		lb = optim_data->lb;
		for (int i = 0; i < NDOF; i++) {
			limit_avg[i] = (ub[i] + lb[i])/2.0;
		}
	}

	for (int i = 0; i < NDOF; i++) {
		cost += pow(x[i] - limit_avg[i],2);
		cost += pow(x[i + NDOF], 2);
	}
	return cost;
}


/*
 * Constant function for simple inverse kinematics
 */
static double const_costfunc(unsigned n, const double *x, double *grad, void *my_func_params) {

	return 1.0;
}

/*
 * Initialize solution for NLOPT optimizer based
 * 2*NDOF dimensional inverse kinematics using last solution found
 *
 */
static void init_last_soln(const optim * params, double x[2*NDOF]) {

	// initialize first dof entries to q0
	for (int i = 0; i < NDOF; i++) {
		x[i] = params->qf[i];
		x[i+NDOF] = params->qfdot[i];
	}
}

/*
 * Initialize solution for NLOPT optimizer based
 * 2*NDOF dimensional inverse kinematics using rest posture
 *
 * The closer to the optimum it is the faster alg should converge
 */
static void init_rest_soln(const coptim * params, double x[2*NDOF]) {

	// initialize first dof entries to q0
	for (int i = 0; i < NDOF; i++) {
		x[i] = params->qrest[i];
		x[i+NDOF] = 0.0;
	}
}

/*
 * This is the constraint that makes sure we hit the ball
 */
static void kinematics_eq_constr(unsigned m, double *result, unsigned n,
		                  const double *x, double *grad, void *my_function_data) {

	static double pos[NCART];
	static double qfdot[NDOF];
	static double vel[NCART];
	static double normal[NCART];
	static double q[NDOF];
	static racketdes* racket_data;

	/* initialization of static variables */
	if (firsttime[1]) {
		firsttime[1] = false;
		racket_data = (racketdes*) my_function_data;
	}

	// extract state information from optimization variables
	for (int i = 0; i < NDOF; i++) {
		q[i] = x[i];
		qfdot[i] = x[i+NDOF];
	}

	// compute the actual racket pos,vel and normal
	calc_racket_state(q,qfdot,pos,vel,normal);

	// deviations from the desired racket frame
	for (int i = 0; i < NCART; i++) {
		result[i] = pos[i] - racket_data->pos[i][0];
		result[i + NCART] = vel[i] - racket_data->vel[i][0];
		result[i + 2*NCART] = normal[i] - racket_data->normal[i][0];
	}

}

/*
 * Debug by testing the constraint violation of the solution vector
 *
 */
static double test_optim(const double *x, const double T,
		          coptim * coparams,
				  racketdes * racketdata) {

	double x_[2*NDOF+1];
	for (int i = 0; i < 2*NDOF; i++)
		x_[i] = x[i];
	x_[2*NDOF] = T;

	// give info on constraint violation
	double *grad = 0;
	static double kin_violation[EQ_CONSTR_DIM];
	static double lim_violation[INEQ_CONSTR_DIM]; // joint limit violations on strike and return
	kinematics_eq_constr(EQ_CONSTR_DIM, kin_violation,
			             2*NDOF, x, grad, racketdata);
	joint_limits_ineq_constr(INEQ_CONSTR_DIM, lim_violation,
			                 2*NDOF, x_, grad, coparams);
	//double cost = costfunc(OPTIM_DIM, x, grad, coparams);

	if (coparams->verbose) {
		// give info on solution vector
		print_optim_vec(x_);
		//printf("f = %.2f\n",cost);
		printf("Position constraint violation: [%.2f %.2f %.2f]\n",kin_violation[0],kin_violation[1],kin_violation[2]);
		printf("Velocity constraint violation: [%.2f %.2f %.2f]\n",kin_violation[3],kin_violation[4],kin_violation[5]);
		printf("Normal constraint violation: [%.2f %.2f %.2f]\n",kin_violation[6],kin_violation[7],kin_violation[8]);
		for (int i = 0; i < INEQ_CONSTR_DIM; i++) {
			if (lim_violation[i] > 0.0)
				printf("Joint limit violated by %.2f on joint %d\n", lim_violation[i], i % NDOF + 1);
		}
	}

	return fmax(max_abs_array(kin_violation,EQ_CONSTR_DIM),
			    max_array(lim_violation,INEQ_CONSTR_DIM));
}

/*
 * This is the inequality constraint that makes sure we never exceed the
 * joint limits during the striking and returning motion
 *
 */
static void joint_limits_ineq_constr(unsigned m, double *result,
		unsigned n, const double *x, double *grad, void *my_func_params) {

	static double a1[NDOF];
	static double a2[NDOF];
	static double a1ret[NDOF]; // coefficients for the returning polynomials
	static double a2ret[NDOF];
	static double qdot_rest[NDOF];
	static double joint_strike_max_cand[NDOF];
	static double joint_strike_min_cand[NDOF];
	static double joint_return_max_cand[NDOF];
	static double joint_return_min_cand[NDOF];
	static double *q0;
	static double *q0dot;
	static double *qrest;
	static double *ub;
	static double *lb;
	static double Tret;
	static bool firsttime = true;

	if (firsttime) {
		firsttime = false;
		coptim* optim_data = (coptim*)my_func_params;
		q0 = optim_data->q0;
		q0dot = optim_data->q0dot;
		qrest = optim_data->qrest;
		ub = optim_data->ub;
		lb = optim_data->lb;
		Tret = optim_data->time2return;
	}

	// calculate the polynomial coeffs which are used for checking joint limits
	calc_strike_poly_coeff(q0,q0dot,x,a1,a2);
	calc_return_poly_coeff(qrest,qdot_rest,x,Tret,a1ret,a2ret);
	// calculate the candidate extrema both for strike and return
	calc_strike_extrema_cand(a1,a2,x[2*NDOF],q0,q0dot,
			joint_strike_max_cand,joint_strike_min_cand);
	calc_return_extrema_cand(a1ret,a2ret,x,Tret,joint_return_max_cand,joint_return_min_cand);

	/* deviations from joint min and max */
	for (int i = 0; i < NDOF; i++) {
		result[i] = joint_strike_max_cand[i] - ub[i];
		result[i+NDOF] = lb[i] - joint_strike_min_cand[i];
		result[i+2*NDOF] = joint_return_max_cand[i] - ub[i];
		result[i+3*NDOF] = lb[i] - joint_return_min_cand[i];
		//printf("%f %f %f %f\n", result[i],result[i+DOF],result[i+2*DOF],result[i+3*DOF]);
	}

}

/*
 * Calculate the polynomial coefficients from the optimized variables qf,qfdot,T
 * p(t) = a1*t^3 + a2*t^2 + a3*t + a4
 */
static void calc_strike_poly_coeff(const double *q0, const double *q0dot, const double *x,
		                    double *a1, double *a2) {

	double T = x[2*NDOF];

	for (int i = 0; i < NDOF; i++) {
		a1[i] = (2/pow(T,3))*(q0[i]-x[i]) + (1/(T*T))*(q0dot[i] + x[i+NDOF]);
		a2[i] = (3/(T*T))*(x[i]-q0[i]) - (1/T)*(x[i+NDOF] + 2*q0dot[i]);
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

	for (int i = 0; i < NDOF; i++) {
		a1[i] = (2/pow(T,3))*(x[i]-q0[i]) + (1/(T*T))*(q0dot[i] + x[i+NDOF]);
		a2[i] = (3/(T*T))*(q0[i]-x[i]) - (1/T)*(2*x[i+NDOF] + q0dot[i]);
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

	for (int i = 0; i < NDOF; i++) {
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

	for (int i = 0; i < NDOF; i++) {
		cand1 = fmin(Tret, fmax(0,(-a2[i] + sqrt(a2[i]*a2[i] - 3*a1[i]*x[i+NDOF]))/(3*a1[i])));
		cand2 =  fmin(Tret, fmax(0,(-a2[i] - sqrt(a2[i]*a2[i] - 3*a1[i]*x[i+NDOF]))/(3*a1[i])));
		cand1 = a1[i]*pow(cand1,3) + a2[i]*pow(cand1,2) + x[i+NDOF]*cand1 + x[i];
		cand2 = a1[i]*pow(cand2,3) + a2[i]*pow(cand2,2) + x[i+NDOF]*cand2 + x[i];
		joint_max_cand[i] = fmax(cand1,cand2);
		joint_min_cand[i] = fmin(cand1,cand2);
	}
}

/*
 * Finalize the solution and update target SL structure and hitTime value
 */
static void finalize_soln(const double* x, optim * params, bool detach, double time_elapsed) {

	// initialize first dof entries to q0
	for (int i = 0; i < NDOF; i++) {
		params->qf[i] = x[i];
		params->qfdot[i] = x[i+NDOF];
	}
	if (detach)
		params->T -= (time_elapsed/1e3);
	params->update = true;
}
