/**
 * @file lazyoptim.cpp
 *
 * @brief NLOPT polynomial optimization functions for LAZY PLAYER are stored here.
 *
 * Lazy Player makes use of FIXED PLAYER to seed good trajectories.
 *
 *  Created on: Sep 11, 2016
 *      Author: okan
 */

#include <armadillo>
#include <iostream>
#include "constants.h"
#include "table.h"
#include "utils.h"
#include "stdlib.h"
#include "math.h"
#include "kinematics.h"
#include "optim.h"

#define INEQ_LAND_CONSTR_DIM 11
#define INEQ_JOINT_CONSTR_DIM 2*NDOF + 2*NDOF

static double costfunc(unsigned n, const double *x, double *grad, void *my_func_params);
static double punish_land_robot(const double *xland,
								const double xnet,
								const double *Rland,
								const double Rnet);
static void land_ineq_constr(unsigned m, double *result, unsigned n, const double *x, double *grad,
		                     void *my_func_params);
static void modify_ball_outgoing_vel(double* ballVel);
static void calc_times(const LazyOptim* opt,
		               const double *x,
					   double *netTime,
					   double *landTime,
					   double *xnet,
					   double *xland,
					   double *ball_norm, // ball to racket normal distance
					   double *ball_proj);
static void racket_contact_model(const double* racketVel,
		                         const double* racketNormal, double* ballVel);
static void interp_ball(const optim_des *params, const double T, double *ballpos, double *ballvel);


LazyOptim::LazyOptim(double qrest_[NDOF], double lb_[NDOF], double ub_[NDOF])
                          : FocusedOptim() {

	const_vec(NDOF,1.0,w.R_strike);
	const_vec(NCART,2e4,w.R_hit);
	const_vec(2,1e3,w.R_land);
	w.R_net = 1e2;

	double tol_ineq_land[INEQ_LAND_CONSTR_DIM];
	double tol_ineq_joint[INEQ_JOINT_CONSTR_DIM];
	const_vec(INEQ_LAND_CONSTR_DIM,1e-2,tol_ineq_land);
	const_vec(INEQ_JOINT_CONSTR_DIM,1e-3,tol_ineq_joint);

	// LN = does not require gradients //
	opt = nlopt_create(NLOPT_LN_COBYLA, 2*NDOF+1);
	nlopt_set_lower_bounds(opt, lb_);
	nlopt_set_upper_bounds(opt, ub_);
	nlopt_set_min_objective(opt, costfunc, this);
	nlopt_add_inequality_mconstraint(opt, INEQ_LAND_CONSTR_DIM,
			land_ineq_constr, this, tol_ineq_land);
	nlopt_add_inequality_mconstraint(opt, INEQ_JOINT_CONSTR_DIM,
			joint_limits_ineq_constr, this, tol_ineq_joint);
	nlopt_set_xtol_rel(opt, 1e-2);

	for (int i = 0; i < NDOF; i++) {
		qrest[i] = qrest_[i];
		ub[i] = ub_[i];
		lb[i] = lb_[i];
	}
}


double LazyOptim::test_soln(const double x[]) const {

	// give info on constraint violation
	static int count = 0;
	double *grad = 0;
	static double land_violation[INEQ_LAND_CONSTR_DIM];
	static double lim_violation[INEQ_JOINT_CONSTR_DIM]; // joint limit violations on strike and return
	joint_limits_ineq_constr(INEQ_JOINT_CONSTR_DIM, lim_violation, OPTIM_DIM, x, grad, (void*)this);
	land_ineq_constr(INEQ_LAND_CONSTR_DIM, land_violation, OPTIM_DIM, x, grad, (void*)this);
	double cost = costfunc(OPTIM_DIM, x, grad, (void*)this);

	if (verbose) {
		// give info on solution vector
		printf("Optim count: %d\n", (++count));
		print_optim_vec(x);
		printf("f = %.2f\n",cost);
		printf("Hitting constraints (b2r):\n");
		printf("Distance along normal: %.2f\n",land_violation[0]);
		printf("Distance along racket: %.2f\n", land_violation[2]);
		printf("Landing constraints:\n");
		printf("NetTime: %f\n", -land_violation[3]);
	    printf("LandTime: %f\n", -land_violation[6]-land_violation[3]);
		printf("Below wall by : %f\n", -land_violation[4]);
		printf("Above net by: %f\n", -land_violation[5]);
		printf("X: between table limits by [%f, %f]\n", -land_violation[7], -land_violation[8]);
		printf("Y: between table limits by [%f, %f]\n", -land_violation[9], -land_violation[10]);
		for (int i = 0; i < INEQ_JOINT_CONSTR_DIM; i++) {
			if (lim_violation[i] > 0.0) {
				printf("Joint limit violated by %.2f on joint %d\n", lim_violation[i], i%NDOF + 1);
			}
		}
	}

	return fmax(max_array(lim_violation,INEQ_JOINT_CONSTR_DIM),
			max_array(land_violation,INEQ_LAND_CONSTR_DIM));
}

/*
 * Calculates the cost function for table tennis Lazy Player (LP)
 * to find spline (3rd order strike+return) polynomials
 */
static double costfunc(unsigned n, const double *x, double *grad, void *my_func_params) {

	static double balldist;
	static double ballproj[NCART];
	static double J1, Jhit, Jland;
	static double a1[NDOF], a2[NDOF];
	static double netTime, landTime, xnet, xland[2];
	double T = x[2*NDOF];

	LazyOptim* opt = (LazyOptim*) my_func_params;
	double *q0 = opt->q0;
	double *q0dot = opt->q0dot;
	weights w = opt->w;

	// calculate the polynomial coeffs which are used in the cost calculation
	calc_strike_poly_coeff(q0,q0dot,x,a1,a2);

	// calculate the landing time
	calc_times(opt, x, &netTime, &landTime, &xnet, xland, &balldist, ballproj);

	J1 = T * (3*T*T*inner_winv_prod(NDOF,w.R_strike,a1,a1) +
			3*T*inner_winv_prod(NDOF,w.R_strike,a1,a2) +
			inner_winv_prod(NDOF,w.R_strike,a2,a2));

	//Jhit = inner_w_prod(NCART,w.R_hit,ballproj,ballproj);
	//Jland = punish_land_robot(xland,xnet,w.R_land, w.R_net);

	//std::cout << J1 << "\t" << Jhit << "\t" << Jland << std::endl;

	return J1 + Jhit + Jland;
}

/*
 * Punish the robot sufficiently as to induce proper landing behaviour
 *
 * The desired landing location chosen is the center of the table
 */
static double punish_land_robot(const double *xland,
								const double xnet,
								const double *Rland,
								const double Rnet) {
	// desired landing locations
	static double x_des_land = 0.0;
	static double y_des_land = dist_to_table - 3*table_length/4;
	static double z_des_net = floor_level - table_height + net_height + 1.0;

	return sqr(xland[X] - x_des_land)*Rland[X] +
			sqr(xland[Y] - y_des_land)*Rland[Y] +
			sqr(xnet - z_des_net)*Rnet;

}

/*
 * This is the constraint that makes sure we land the ball
 *
 */
static void land_ineq_constr(unsigned m, double *result, unsigned n, const double *x, double *grad,
		                     void *my_func_params) {

	static double xnet;
	static double xland[2];
	static double balldist;
	static double ballproj[NCART];
	static double netTime;
	static double landTime;
	static double table_xmax = table_width/2.0;
	static double table_ymax = dist_to_table - table_length;
	static double wall_z = 1.0;
	static double net_y = dist_to_table - table_length/2.0;
	static double net_z = floor_level - table_height + net_height;

	calc_times((LazyOptim*)my_func_params, x, &netTime, &landTime, &xnet, xland, &balldist, ballproj);

	result[0] = -balldist;
	result[1] = balldist - ball_radius;
	result[2] = inner_prod(NCART,ballproj,ballproj) - sqr(racket_radius);
	result[3] = -netTime;
	result[4] = xnet - wall_z;
	result[5] = -xnet + net_z;
	result[6] = netTime - landTime;
	result[7] = xland[X] - table_xmax;
	result[8] = -xland[X] - table_xmax;
	result[9] = xland[Y] - net_y;
	result[10] = -xland[Y] + table_ymax;
}

/*
 * Calculate the net hitting time and the landing time
 * Assuming that there was an interaction (hitting) between racket and ball.
 *
 * Since we call this function multiple times within one optimization step
 * we store the optimization variable x,
 * landTime and netTime variables and return these when the optimization
 * variable is the same (meaning optimization algorithm did not update it yet).
 *
 * TODO: simplify code! If landtime is not real, should we set it to -1.0?
 *
 */
static void calc_times(const LazyOptim *opt,
		               const double *x,
					   double *netTime,
					   double *landTime,
					   double *xnet,
					   double *xland,
					   double *ball_norm, // ball to racket normal distance
					   double *ball_proj) { // ball projected to racket plane

	static double x_[2*NDOF+1];
	static double b2rn_;
	static double b2rp_[NCART];
	static double xnet_;
	static double xland_[2];
	static double landTime_;
	static double netTime_; // local static variables for saving state

	static double qf[NDOF];
	static double qfdot[NDOF];
	static double vel[NCART];
	static double normal[NCART];
	static double pos[NCART];
	static double ballpos[NCART];
	static double ballvel[NCART];
	static double table_z = floor_level - table_height + ball_radius;
	static double g = 9.8;
	static double net_y = dist_to_table - table_length/2.0;
	double distBall2TableZ;

	if (vec_is_equal(2*NDOF+1,x,x_)) {
		*netTime = netTime_;
		*landTime = landTime_;
		*xnet = xnet_;
		make_equal(2,xland_,xland);
		*ball_norm = b2rn_;
		make_equal(NCART,b2rp_,ball_proj);
	}
	else {
		// extract state information from optimization variables
		for (int i = 0; i < NDOF; i++) {
			qf[i] = x[i];
			qfdot[i] = x[i + NDOF];
		}
		calc_racket_state(qf, qfdot, pos, vel, normal);
		interp_ball(opt->param_des, x[2*NDOF], ballpos, ballvel);
		racket_contact_model(vel, normal, ballvel);
		modify_ball_outgoing_vel(ballvel);
		distBall2TableZ = ballpos[Z] - table_z;

		if (sqr(ballvel[Z]) > -2*g*distBall2TableZ) {
			landTime_ = *landTime =
			fmax(0.10,(ballvel[Z] + sqrt(sqr(ballvel[Z]) + 2*g*distBall2TableZ))/g);
		}
		else {
			// landTime is not real!
			landTime_ = *landTime = 1.0;
		}
		netTime_ = *netTime = fmax(0.10,(net_y - ballpos[Y])/ballvel[Y]);

		xnet_ = *xnet = ballpos[Z] + netTime_ * ballvel[Z] - 0.5*g*netTime_*netTime_;
		xland_[X] = xland[X] = ballpos[X] + landTime_ * ballvel[X];
		xland_[Y] = xland[Y] = ballpos[Y] + landTime_ * ballvel[Y];

		// calculate deviation of ball to racket - hitting constraints
		make_equal(NCART,ballpos,ball_proj);
		vec_minus(NCART,pos,ball_proj);
		b2rn_ = *ball_norm = inner_prod(NCART,normal,ball_proj);
		for (int i = 0; i < NCART; i++) {
			normal[i] *= b2rn_;
		}
		vec_minus(NCART,normal,ball_proj); // we get (e - nn'e) where e is ballpos - racketpos
		make_equal(NCART,ball_proj,b2rp_);
		make_equal(2*NDOF+1,x,x_);
	}
}

/*
 * Modify ball outgoing velocity by multiplying with a constant vector (less than one)
 *
 * Necessary since the landing and net hitting locations are calculated using projectile motion
 *
 */
static void modify_ball_outgoing_vel(double* ballVel) {

	static double multiplier_x, multiplier_y, multiplier_z;
	static bool firsttime = true;

	if (firsttime) {
		firsttime = false;
		multiplier_x = 0.9;
		multiplier_y = 0.8;
		multiplier_z = 0.83;
	}

	ballVel[X] *= multiplier_x;
	ballVel[Y] *= multiplier_y;
	ballVel[Z] *= multiplier_z;
}

/*
 * Update the incoming ball velocity with outgoing ball velocity using MIRROR LAW
 *
 * The racket contact model in vector form is O = I + (1 + eps_R)*N*N'*(V - I)
 * where I is the incoming ball velocity
 *       N is the racket normal
 *       V is the racket velocity
 *       eps_R is the coefficient of restitution of the racket
 *
 *
 */
static void racket_contact_model(const double* racketVel,
		const double* racketNormal, double* ballVel) {

	static const double racket_param = 0.78;
	static double diffVel[NCART];
	static double normalMultSpeed[NCART];
	double speed;

	for (int i = 0; i < NCART; i++)
		diffVel[i] = racketVel[i] - ballVel[i];

	speed = (1 + racket_param) * inner_prod(NCART, racketNormal, diffVel);

	for (int i = 0; i < NCART; i++) {
		normalMultSpeed[i] = speed * racketNormal[i];
	}
	vec_plus(NCART,normalMultSpeed,ballVel);
}

/*
 * First order hold to interpolate linearly at time T
 * between racket pos,vel,normal entries
 *
 * IF T is nan, racket variables are assigned to zero-element of
 * relevant racket entries
 *
 */
static void interp_ball(const optim_des *data, const double T, double *ballpos, double *ballvel) {

    const double dt = data->dt;
	const int Nmax = data->Nmax;
	int N = (int) (T/dt);
	double Tdiff = T - N*dt;

	for (int i = 0; i < NCART; i++) {
		if (N < Nmax) {
			ballpos[i] = data->ball_pos(i,N) +
					(Tdiff/dt) * (data->ball_pos(i,N+1) - data->ball_pos(i,N));
			ballvel[i] = data->ball_vel(i,N) +
					(Tdiff/dt) * (data->ball_vel(i,N+1) - data->ball_vel(i,N));
		}
		else {
			ballpos[i] = data->ball_pos(i,N);
			ballvel[i] = data->ball_vel(i,N);
		}
	}
}
