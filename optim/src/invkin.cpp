/*
 * invkin.c
 *
 * Inverse Kinematics optimization to find qf, qfdot joint angles
 * and velocities at Virtual Hitting Plane (VHP)
 *
 *  Created on: Mar 5, 2017
 *      Author: okoc
 */

#include "constants.h"
#include "utils.h"
#include "kinematics.h"
#include "stdlib.h"
#include "math.h"
#include "optim.h"

static double penalize_dist_to_limits(unsigned n, const double *x,
		                     double *grad, void *my_func_params);
static double const_costfunc(unsigned n, const double *x,
		                     double *grad, void *my_func_params);
static void init_invkin_soln(const optim * params, double x[2*NDOF]);
static void kinematics_eq_constr(unsigned m, double *result, unsigned n,
		                  const double *x, double *grad, void *my_function_data);
static double test_invkin(const double *x, InvKin * invkin, int info);
static void finalize_soln(const double* x, optim * params, double time_elapsed);

/*
 * Constructor for the Inverse Kinematics class
 * Sets the structures and calls the NLOPT optimization
 *
 */
InvKin::InvKin(coptim *coparams_,
		       racketdes *racketdata_,
		       optim *params_) :
	coparams(coparams_), racketdata(racketdata_), params(params_) {

	double tol_eq[EQ_CONSTR_DIM];
	const_vec(EQ_CONSTR_DIM,1e-2,tol_eq);
	// set tolerances equal to second argument

	// LN = does not require gradients //
	opt = nlopt_create(NLOPT_LN_COBYLA, 2*NDOF);
	nlopt_set_xtol_rel(opt, 1e-2);
	nlopt_set_lower_bounds(opt, coparams->lb);
	nlopt_set_upper_bounds(opt, coparams->ub);
	nlopt_set_min_objective(opt, penalize_dist_to_limits, this);
	nlopt_add_equality_mconstraint(opt, EQ_CONSTR_DIM,
			kinematics_eq_constr, this, tol_eq);

	for (int i = 0; i < NDOF; i++) {
		limit_avg[i] = (coparams->ub[i] + coparams->lb[i])/2.0;
	}
}

/*
 *
 * Run Inverse Kinematics to find a solution for Virtual Hitting Plane (VHP)
 *
 */
double InvKin::operator()() {

	params->update = FALSE;
	params->running = TRUE;

	double x[2*NDOF];
	init_invkin_soln(params,x); //parameters are the initial joint positions q0*/
	double init_time = get_time();
	double past_time = 0.0;
	double minf; // the minimum objective value, upon return //
	int res; // error code
	double max_violation;

	if ((res = nlopt_optimize(opt, x, &minf)) < 0) {
		printf("NLOPT failed with exit code %d!\n", res);
		past_time = (get_time() - init_time)/1e3;
		max_violation = 100.0;
	}
	else {
		past_time = (get_time() - init_time)/1e3;
		printf("NLOPT success with exit code %d!\n", res);
		printf("NLOPT took %f ms\n", past_time);
		printf("Found minimum at f = %0.10g\n", minf);
		max_violation = test_invkin(x,this,true);
		if (max_violation < 1e-2)
			finalize_soln(x,params,past_time);
	}
	params->running = FALSE;
	check_optim_result(res);
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

	double cost = 0.0;
	InvKin * invkin = (InvKin*)my_func_params;

	for (int i = 0; i < NDOF; i++) {
		cost += pow(x[i] - invkin->limit_avg[i],2);
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
 * Estimate an initial solution for NLOPT optimizer based
 * 2*NDOF dimensional inverse kinematics
 *
 * The closer to the optimum it is the faster alg should converge
 */
static void init_invkin_soln(const optim * params, double x[2*NDOF]) {

	// initialize first dof entries to q0
	for (int i = 0; i < NDOF; i++) {
		x[i] = params->qf[i];
		x[i+NDOF] = params->qfdot[i];
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

	InvKin * invkin = (InvKin*)my_function_data;

	// extract state information from optimization variables
	for (int i = 0; i < NDOF; i++) {
		q[i] = x[i];
		qfdot[i] = x[i+NDOF];
	}

	// compute the actual racket pos,vel and normal
	calc_racket_state(q,qfdot,pos,vel,normal);

	// deviations from the desired racket frame
	for (int i = 0; i < NCART; i++) {
		result[i] = pos[i] - invkin->racketdata->pos[i][0];
		result[i + NCART] = vel[i] - invkin->racketdata->vel[i][0];
		result[i + 2*NCART] = normal[i] - invkin->racketdata->normal[i][0];
	}

}

/*
 * Finalize the solution and update target SL structure and hitTime value
 */
static void finalize_soln(const double* x, optim * params, double time_elapsed) {

	// initialize first dof entries to q0
	for (int i = 0; i < NDOF; i++) {
		params->qf[i] = x[i];
		params->qfdot[i] = x[i+NDOF];
	}
	params->T = params->T - (time_elapsed/1e3);
	params->update = TRUE;
}

/*
 * Call test optim function (set hitting time T fixed)
 */
static double test_invkin(const double *x, InvKin * invkin, int info) {

	double x_[2*NDOF+1];
	for (int i = 0; i < 2*NDOF; i++) {
		x_[i] = x[i];
	}
	x_[2*NDOF] = invkin->params->T;
	return test_optim(x_,invkin,TRUE);
}
