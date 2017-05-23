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
#include <thread>
#include "constants.h"
#include "utils.h"
#include "stdlib.h"
#include "math.h"
#include "kinematics.h"
#include "tabletennis.h"
#include "kalman.h"
#include "optim.h"

using namespace arma;

static void calc_racket_strategy(const vec6 & ball_predicted,
		                  const vec2 & ball_land_des,
						  const double time_land_des,
						  racket_des & racket);
static bool predict_hitting_point(vec6 & ball_pred, double & time_pred,
		                   EKF & filter, game & game_state);
static void check_legal_bounce(const vec6 & ball_est, game & game_state);
static bool check_legal_ball(const vec6 & ball_est, const mat & balls_predicted, game & game_state);
static void predict_ball(const double & time_pred, mat & balls_pred, EKF & filter);

static double penalize_dist_to_limits(unsigned n, const double *x,
		                     double *grad, void *my_func_params);
static double const_costfunc(unsigned n, const double *x,
		                     double *grad, void *my_func_params);
static void kinematics_eq_constr(unsigned m, double *result, unsigned n,
		                  const double *x, double *grad, void *my_function_data);
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


HittingPlane::HittingPlane(double qrest_[NDOF], double lb[NDOF], double ub[NDOF]) {

	time_land_des = 0.8;
	ball_land_des[X] = 0.0;
	ball_land_des[Y] = dist_to_table - 3*table_length/4;
	double tol_eq[EQ_CONSTR_DIM];
	const_vec(EQ_CONSTR_DIM,1e-2,tol_eq);
	// set tolerances equal to second argument

	// LN = does not require gradients //
	opt = nlopt_create(NLOPT_LN_COBYLA, 2*NDOF);
	nlopt_set_xtol_rel(opt, 1e-2);
	nlopt_set_lower_bounds(opt, lb);
	nlopt_set_upper_bounds(opt, ub);
	nlopt_set_min_objective(opt, penalize_dist_to_limits, this);
	nlopt_add_equality_mconstraint(opt, EQ_CONSTR_DIM,
			kinematics_eq_constr, this, tol_eq);

	for (int i = 0; i < NDOF; i++) {
		qrest[i] = qrest_[i];
		limit_avg[i] = (ub[i] + lb[i])/2.0;
	}
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

	// initialize first dof entries to q0
	for (int i = 0; i < NDOF; i++) {
		qf[i] = x[i];
		qfdot[i] = x[i+NDOF];
	}
	if (detach)
		T -= (time_elapsed/1e3);
	update = true;
}

double HittingPlane::test_soln(const double x[2*NDOF]) const {

	double x_[2*NDOF+1];
	for (int i = 0; i < 2*NDOF; i++)
		x_[i] = x[i];
	x_[2*NDOF] = T;

	// give info on constraint violation
	double *grad = 0;
	static double kin_violation[EQ_CONSTR_DIM];
	static double lim_violation[INEQ_CONSTR_DIM]; // joint limit violations on strike and return
	kinematics_eq_constr(EQ_CONSTR_DIM, kin_violation, 2*NDOF,
			             x, grad, (void*)this);
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

bool HittingPlane::predict(EKF & filter) {

	bool flag = false;
	vec6 ball_pred;
	double time_pred;
	if (predict_hitting_point(ball_pred,time_pred,filter,game_state)) {
		calc_racket_strategy(ball_pred,ball_land_des,time_land_des,racket);
		T = time_pred;
		flag = true;
	}
	return flag;
}

void HittingPlane::optim() {

	update = false;
	running = true;
	double x[2*NDOF];
	double tol_eq[EQ_CONSTR_DIM];
	const_vec(EQ_CONSTR_DIM,1e-2,tol_eq);

	if (moving)
		init_last_soln(x);
	else
		init_rest_soln(x);

	double init_time = get_time();
	double past_time = 0.0;
	double minf; // the minimum objective value, upon return //
	int res; // error code
	double max_violation;

	if ((res = nlopt_optimize(opt, x, &minf)) < 0) {
		if (verbose)
			printf("NLOPT failed with exit code %d!\n", res);
	    past_time = (get_time() - init_time)/1e3;
	    max_violation = 100.0;
	}
	else {
		past_time = (get_time() - init_time)/1e3;
		if (verbose) {
			printf("NLOPT success with exit code %d!\n", res);
			printf("NLOPT took %f ms\n", past_time);
			printf("Found minimum at f = %0.10g\n", minf);
		}
	    max_violation = test_soln(x);
	    if (max_violation < 1e-2)
	    	finalize_soln(x,past_time);
	}
	running = false;
}

static void calc_racket_strategy(const vec6 & ball_predicted,
		                  const vec2 & ball_land_des,
						  const double time_land_des,
						  racket_des & racket) {

	TableTennis tennis = TableTennis(false,false);

	vec3 balls_out_vel;
	vec3 racket_des_normal;
	vec3 racket_des_vel;
	tennis.calc_des_ball_out_vel(ball_land_des,time_land_des,ball_predicted,balls_out_vel);
	tennis.calc_des_racket_normal(ball_predicted.rows(DX,DZ),balls_out_vel,racket_des_normal);
	tennis.calc_des_racket_vel(ball_predicted.rows(DX,DZ),balls_out_vel,racket_des_normal,racket_des_vel);

	// place racket centre on the predicted ball

	for (int j = 0; j < NCART; j++) {
		racket.pos[j] = ball_predicted(j);
		racket.vel[j] = racket_des_vel(j);
		racket.normal[j] = racket_des_normal(j);
	}
}

/*
 * Predict hitting point on the Virtual Hitting Plane (VHP)
 * if the ball is legal (legal detected bounce or legal predicted bounce)
 * and there is enough time (50 ms threshold)
 *
 * The location of the VHP is defined as a constant (constants.h)
 *
 */
static bool predict_hitting_point(vec6 & ball_pred, double & time_pred,
		                   EKF & filter, game & game_state) {

	const double time_min = 0.05;
	mat balls_path;
	bool valid_hp = false;
	predict_ball(time_pred,balls_path,filter);
	uvec vhp_index;
	unsigned idx;

	if (check_legal_ball(filter.get_mean(),balls_path,game_state)) {
		vhp_index = find(balls_path.row(Y) >= VHPY, 1);
		if (vhp_index.n_elem == 1) {
			idx = as_scalar(vhp_index);
			ball_pred = balls_path.col(idx);
			time_pred = DT * (idx + 1);
			if (time_pred > time_min)
				valid_hp = true;
		}
	}

	return valid_hp;
}

/*
 * Predict ball with the models fed into the filter
 *
 * Number of prediction steps is given by Nmax in racket
 * parameters
 *
 */
static void predict_ball(const double & time_pred, mat & balls_pred, EKF & filter) {

	int N = (int)(time_pred/DT);
	balls_pred = filter.predict_path(DT,N);
}

/*
 *
 * Checks for legal bounce
 * If an incoming ball has bounced before or bounces on opponents' court
 * it is declared ILLEGAL (legal_bounce as DATA MEMBER of player class)
 *
 * Bounce variable is static variable of estimate_ball_state method of player class
 * which is reset each time an incoming ball from ball gun is detected.
 *
 * TODO: also consider detecting HIT by robot racket
 *
 */
static void check_legal_bounce(const vec6 & ball_est, game & game_state) {


	static double last_z_vel = 0.0;
	bool ball_is_incoming = ball_est(DY) > 0.0;
	bool on_opp_court = (ball_est(Y) < (dist_to_table - (table_length/2.0)));

	if (last_z_vel < 0.0 && ball_est(DZ) > 0.0) { // bounce must have occurred
		if (ball_is_incoming) {
			if (game_state == LEGAL || on_opp_court) {
				game_state = ILLEGAL;
			}
			else {
				//cout << "Legal bounce occurred!" << endl;
				game_state = LEGAL;
			}
		}
	}

	last_z_vel = ball_est(DZ);
}

/*
 * Check if the table tennis trial is LEGAL (hence motion planning can be started).
 *
 * If it is exactly one ball when in awaiting mode
 * (before an actual legal bounce was detected) trial will be legal!
 *
 * If ball has bounced legally bounce, then there should be no more bounces.
 *
 * TODO: no need to check after ball passes table
 *
 */
static bool check_legal_ball(const vec6 & ball_est, const mat & balls_predicted, game & game_state) {

	int num_bounces = 0;
	int N = balls_predicted.n_cols;

	check_legal_bounce(ball_est, game_state);

	// if sign of z-velocity changes then the ball bounces
	for (int i = 0; i < N-1; i++) {
		if (balls_predicted(DZ,i) < 0 && balls_predicted(DZ,i+1) > 0) {
			num_bounces++;
		}
	}

	// one bounce is predicted before an actual bounce has happened
	if (game_state == AWAITING && num_bounces == 1) {
		return true;
	}
	// no bounce is predicted
	if (game_state == LEGAL && num_bounces == 0) {
		return true;
	}

	return false;
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
 */
static void kinematics_eq_constr(unsigned m, double *result, unsigned n,
		                  const double *x, double *grad, void *my_function_data) {

	static double pos[NCART];
	static double qfdot[NDOF];
	static double vel[NCART];
	static double normal[NCART];
	static double q[NDOF];

	HittingPlane * vhp = (HittingPlane*)my_function_data;

	// extract state information from optimization variables
	for (int i = 0; i < NDOF; i++) {
		q[i] = x[i];
		qfdot[i] = x[i+NDOF];
	}

	// compute the actual racket pos,vel and normal
	calc_racket_state(q,qfdot,pos,vel,normal);

	// deviations from the desired racket frame
	for (int i = 0; i < NCART; i++) {
		result[i] = pos[i] - vhp->racket.pos[i];
		result[i + NCART] = vel[i] - vhp->racket.vel[i];
		result[i + 2*NCART] = normal[i] - vhp->racket.normal[i];
	}
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
