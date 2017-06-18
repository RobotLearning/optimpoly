/**
 * @file invkin.cpp
 *
 * @brief Inverse Kinematics optimization to find qf, qfdot joint angles
 * and velocities at Virtual Hitting Plane (VHP).
 *
 *  Created on: Mar 5, 2017
 *      Author: okoc
 */

#include <armadillo>
#include "constants.h"
#include "utils.h"
#include "stdlib.h"
#include "math.h"
#include "kinematics.h"
#include "optim.h"
#include <adolc/adouble.h>            // use of active doubles
#include <adolc/drivers/drivers.h>    // use of "Easy to Use" drivers
#include <adolc/taping.h>             // use of taping

static double penalize_dist_to_limits(unsigned n, const double *x,
		                     double *grad, void *my_func_params);
static double const_costfunc(unsigned n, const double *x,
		                     double *grad, void *my_func_params);
static void kinematics_eq_constr(unsigned m, double *result, unsigned n,
		                  const double *x, double *grad, void *my_function_data);
static void kinematics_eq_constr(unsigned m, adouble *result, unsigned n,
		                  const adouble *x, double *grad, void *my_function_data);

/*
 * This is actually only used by Virtual Hitting Plane!
 */
void Optim::fix_hitting_time(double time_pred) {
	if (time_pred > 0.05)
		T = time_pred;
}

HittingPlane::HittingPlane(double qrest_[], double lb_[],
						   double ub_[], bool grad_) : GRAD_BASED_OPT(grad_) {

	double tol_eq[EQ_CONSTR_DIM];
	const_vec(EQ_CONSTR_DIM,1e-2,tol_eq);
	// set tolerances equal to second argument

	if (GRAD_BASED_OPT) {
		opt = nlopt_create(NLOPT_LD_MMA, OPTIM_DIM);
		param_des = new optim_des();
		generate_tape();
	}
	else {	// LN = does not require gradients //
		param_des = nullptr;
		opt = nlopt_create(NLOPT_LN_COBYLA, OPTIM_DIM);
	}
	nlopt_set_xtol_rel(opt, 1e-2);
	nlopt_set_lower_bounds(opt, lb_);
	nlopt_set_upper_bounds(opt, ub_);
	nlopt_set_min_objective(opt, penalize_dist_to_limits, this);
	nlopt_add_equality_mconstraint(opt,EQ_CONSTR_DIM,kinematics_eq_constr,this,tol_eq);

	for (int i = 0; i < NDOF; i++) {
		qrest[i] = qrest_[i];
		ub[i] = ub_[i];
		lb[i] = lb_[i];
		limit_avg[i] = (ub[i] + lb[i])/2.0;
	}
}

/*
 * Release auto-diff tape (jacobian double-pointer)
 */
HittingPlane::~HittingPlane() {

	//nlopt_destroy(opt);
	if (GRAD_BASED_OPT) {
		for (int i = 0; i < EQ_CONSTR_DIM; i++)
		    delete[] jac[i];
		delete[] jac;
	}
}

/*
 * Initialize ADOLC operations
 * Objective has easily computed gradients so no need to auto-diff
 */
void HittingPlane::generate_tape() {

	int n = OPTIM_DIM;
	int m = EQ_CONSTR_DIM;

	adouble *x_auto = new adouble[n];
	adouble *g_auto = new adouble[m];
	double *x = new double[n];
	double dummy;
	jac = new double*[m];
	for (int i = 0; i < m; i++)
		jac[i] = new double[n];

	init_rest_soln(x);

	trace_on(2);
	for (int i = 0; i < n; i++)
		x_auto[i] <<= x[i];

	kinematics_eq_constr(m,g_auto,n,x_auto,nullptr,this);
	for(int i = 0; i < m; i++)
		g_auto[i] >>= dummy;
	trace_off();

	delete[] x_auto;
	delete[] x;
	delete[] g_auto;
}

void HittingPlane::init_last_soln(double x[2*NDOF]) const {

	// initialize first dof entries to q0
	for (int i = 0; i < NDOF; i++) {
		x[i] = qf[i];
		x[i+NDOF] = qfdot[i];
	}
}

void HittingPlane::init_rest_soln(double x[2*NDOF]) const {

	// initialize first dof entries to q0
	for (int i = 0; i < NDOF; i++) {
		x[i] = qrest[i];
		x[i+NDOF] = 0.0;
	}
}

void HittingPlane::finalize_soln(const double x[2*NDOF], double time_elapsed) {

	if (T > 0.05) {
		// initialize first dof entries to q0
		for (int i = 0; i < NDOF; i++) {
			qf[i] = x[i];
			qfdot[i] = x[i+NDOF];
		}
		if (detach) {
			T -= (time_elapsed/1e3);
		}
		update = true;
	}
}

double HittingPlane::test_soln(const double x[]) {

	double x_[2*NDOF+1];
	for (int i = 0; i < 2*NDOF; i++)
		x_[i] = x[i];
	x_[2*NDOF] = T;

	// give info on constraint violation
	double *grad = 0;
	double kin_violation[EQ_CONSTR_DIM];
	double lim_violation[INEQ_CONSTR_DIM]; // joint limit violations on strike and return
	bool flag_temp = false;
	if (GRAD_BASED_OPT) {
		GRAD_BASED_OPT = false;
		flag_temp = true;
	}
	kinematics_eq_constr(EQ_CONSTR_DIM, kin_violation, 2*NDOF,
			             x, grad, (void*)this);
	if (flag_temp)
		GRAD_BASED_OPT = true;

	joint_limits_ineq_constr(INEQ_CONSTR_DIM, lim_violation,
			                 2*NDOF, x_, grad, (void*)this);
	//double cost = costfunc(OPTIM_DIM, x, grad, coparams);

	if (verbose) {
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
 * Penalize (unweighted) squared distance (in joint space) to joint limit averages
 *
 * Adds also joint velocity penalties
 *
 */
static double penalize_dist_to_limits(unsigned n, const double *x, double *grad, void *my_func_params) {

	double cost = 0.0;
	HittingPlane * vhp = (HittingPlane*)my_func_params;

	if (grad) {
		for (int i = 0; i < NDOF; i++) {
			grad[i] = 2 * (x[i] - vhp->limit_avg[i]);
			grad[i+NDOF] = 2 * x[i+NDOF];
		}
	}

	for (int i = 0; i < NDOF; i++) {
		cost += pow(x[i] - vhp->limit_avg[i],2);
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
 * This is the constraint that makes sure we hit the ball
 * AUTO version includes the gradients computed through auto-diff
 */
static void kinematics_eq_constr(unsigned m, double *result, unsigned n,
		                  const double *x, double *grad, void *my_function_data) {

	static double pos[NCART];
	static double vel[NCART];
	static double normal[NCART];
	static double qfdot[NDOF];
	static double qf[NDOF];

	HittingPlane * vhp = (HittingPlane*)my_function_data;
	if (vhp->GRAD_BASED_OPT) {
        jacobian(2,m,n,x,vhp->jac);
        int idx = 0;
        for(int i = 0; i < m; i++) {
        	for(int j = 0; j < n; j++) {
        		grad[idx++] = vhp->jac[i][j];
        		cout << vhp->jac[i][j] << "\t";
        	}
        	cout << endl;
        }
        cout << endl;
	}

	// extract state information from optimization variables
	for (int i = 0; i < NDOF; i++) {
		qf[i] = x[i];
		qfdot[i] = x[i+NDOF];
	}

	// compute the actual racket pos,vel and normal
	calc_racket_state(qf,qfdot,pos,vel,normal);

	// deviations from the desired racket frame
	for (int i = 0; i < NCART; i++) {
		result[i] = pos[i] - vhp->param_des->racket_pos(i);
		result[i + NCART] = vel[i] - vhp->param_des->racket_vel(i);
		result[i + 2*NCART] = normal[i] - vhp->param_des->racket_normal(i);
	}
}


/*
 * This is the constraint that makes sure we hit the ball
 */
static void kinematics_eq_constr(unsigned m, adouble *result, unsigned n,
		                  const adouble *x, double *grad, void *my_function_data) {

	static adouble pos[NCART];
	static adouble vel[NCART];
	static adouble normal[NCART];
	static adouble qfdot[NDOF];
	static adouble qf[NDOF];

	HittingPlane * vhp = (HittingPlane*)my_function_data;

	// extract state information from optimization variables
	for (int i = 0; i < NDOF; i++) {
		qf[i] = x[i];
		qfdot[i] = x[i+NDOF];
	}

	// compute the actual racket pos,vel and normal
	calc_racket_state(qf,qfdot,pos,vel,normal);

	// deviations from the desired racket frame
	for (int i = 0; i < NCART; i++) {
		result[i] = pos[i] - vhp->param_des->racket_pos(i);
		result[i + NCART] = vel[i] - vhp->param_des->racket_vel(i);
		result[i + 2*NCART] = normal[i] - vhp->param_des->racket_normal(i);
	}
}