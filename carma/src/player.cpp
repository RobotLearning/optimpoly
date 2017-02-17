/*
 * player.cpp
 *
 * Table Tennis Player Class
 *
 *  Created on: Feb 8, 2017
 *      Author: okoc
 */

#include <armadillo>
#include "constants.h"
#include "tabletennis.h"
#include "kalman.h"
#include "player.hpp"
#include "optim.hpp"
#include <thread>

using namespace std;
using namespace arma;

/*
 * Initialize Centred Table Tennis Player (CP)
 *
 * Tries to return the ball to the centre of the opponents court
 *
 */
Player::Player(const vec7 & q0, const EKF & filter) : filter(filter) {

	time_land_des = 0.8;
	time2return = 1.0;

	ball_land_des(X) = 0.0;
	ball_land_des(Y) = dist_to_table - 3*table_length/4;

	q_rest_des = q0;
	optim_params = {q0, zeros<vec>(7), 1.0, false};
	racket_params = {zeros<mat>(3,1), zeros<mat>(3,1), zeros<mat>(3,1)};
	moving = false;

}

/*
 * Empirical Bayes procedure to estimate prior given
 * matrix of observations Y of column length N, each row
 * corresponding to one
 *
 * If observations arrive as column vectors then we take
 * transpose of it.
 *
 *
 */
void Player::estimate_prior(const mat & observations, const vec & times) {

	vec6 x; mat66 P;
	int num_samples = times.n_elem;
	mat M = zeros<mat>(num_samples,3);

	// and create the data matrix
	for (int i = 0; i < num_samples; i++) {
		M(i,0) = 1.0;
		M(i,1) = times(i);
		M(i,2) = times(i) * times(i);
	}
	// solving for the parameters
	//cout << "Data matrix:" << endl << M << endl;
	mat Beta = solve(M,observations.t());
	//cout << "Parameters:" << endl << Beta << endl;
	x = join_horiz(Beta.row(0),Beta.row(1)).t(); //vectorise(Beta.rows(0,1));
	P.eye(6,6);
	P *= 0.1;
	filter.set_prior(x,P);
	filter.update(observations.col(0));

	double dt;
	for (int i = 1; i < times.n_elem; i++) {
		dt = times(i) - times(i-1);
		filter.predict(dt);
		filter.update(observations.col(i));
	}
}

/*
 * Filter the blob information with a kalman Filter and save the result in filtState.
 * KF is only used in simulation mode.
 *
 * Returns the bounce variable.
 * It will be reset to FALSE when filter is re-initialized.
 * Bounce variable is used for legal ball detection
 *
 * FIXME: Add automatic resetting!
 *
 */
void Player::estimate_ball_state(const vec3 & obs) {

	static wall_clock timer;
	static bool firsttime = true;
	static const int min_obs = 5;
	static int num_obs = 0;
	// observation matrix
	static mat OBS = zeros<mat>(3,min_obs);
	static vec TIMES = zeros<vec>(min_obs);
	static double t_cum;

	if (firsttime) {
		firsttime = false;
		timer.tic();
	}

	// get elapsed time
	double dt = timer.toc();
	t_cum += dt;

	/*if (*reset) {
		num_obs = 0;
		*reset = false;
		t_cum = 0.0; // t_cumulative
	}*/

	bool new_ball = check_new_obs(obs);
	if (num_obs < min_obs) {
		if (new_ball) {
			TIMES(num_obs) = t_cum;
			OBS.col(num_obs) = obs;
			num_obs++;
			if (num_obs == min_obs) {
				estimate_prior(OBS,TIMES);
			}
		}
	}
	else { // comes here if there are enough balls to start filter
		filter.predict(dt);
		if (new_ball)
			filter.update(obs);
	}
	timer.tic();
}

/*
 * Public interface, suitable for testing
 *
 */
vec6 Player::filt_ball_state(const vec3 & obs) {

	estimate_ball_state(obs);
	return filter.get_mean();
}

/*
 * Play Table Tennis
 */
joint Player::play(const joint & qact, const vec3 & obs) {

	joint qdes = {q_rest_des, zeros<vec>(NDOF), zeros<vec>(NDOF)};

	estimate_ball_state(obs);

	// initialize optimization and get the hitting parameters
	calc_optim_param();

	// generate movement or calculate next desired step
	calc_next_state(qact, qdes);

	return qdes;

}

/*
 * Calculate the optimization parameters using an NLOPT nonlinear optimization algorithm
 * in another thread
 *
 * The optimized parameters are: qf, qf_dot, T
 * assuming T_return and q0 are fixed
 *
 *
 */
void Player::calc_optim_param() {

	int N = 100;
	double dt = 0.02;
	vec6 state_est = filter.get_mean();

	// if ball is fast enough and robot is not moving consider optimization
	if (!moving && state_est(Y) > (dist_to_table - table_length/2) && state_est(DY) > 1.0) {

		mat balls_pred = filter.predict_path(dt,N);

		if (check_legal_ball(balls_pred)) { // ball is legal
			moving = true;
			calc_racket_strategy(balls_pred);
			opt minimizer = opt(LN_COBYLA, OPTIM_DIM);
			PolyOptim polyopt = PolyOptim(q_rest_des,racket_params,
					                  time2return,optim_params,minimizer);
			polyopt.setup();
			// run optimization in another thread
			thread my_thread(polyopt);
		}
	}
}

/*
 * Unfold the next desired state of the 3rd order polynomials in joint space
 * If movement finishes then the desired state velocities and accelerations are zeroed.
 *
 * Multithreading : if after initial lookup, the computation in the other thread terminates, then we
 * synchonize the values between the threads and compute the next desired states given these new polynomial
 * parameters (qf, qf_dot and T_hit)
 *
 */
void Player::calc_next_state(const joint & qact, joint & qdes) {

	static int idx = 0;
	static mat Q_des, Qd_des, Qdd_des;

	// this should be only for MPC?
	if (optim_params.update) {
		moving = true;
		optim_params.update = false;
		generate_strike(qact,Q_des,Qd_des,Qdd_des);
		// call polynomial generation
	}

	// make sure we update after optim finished
	if (moving) {
		qdes.q = Q_des.col(idx);
		qdes.qd = Qd_des.col(idx);
		qdes.qdd = Qdd_des.col(idx);
		idx++;
		if (idx == Q_des.n_cols) {
			// hitting process will finish
			moving = false;
			qdes.q = zeros<vec>(7);
			qdes.qd = zeros<vec>(7);
			qdes.qdd = zeros<vec>(7);
		}
	}

}

/*
 * Create batch hitting and returning joint state 3rd degree polynomials
 *
 */
void Player::generate_strike(const joint & qact, mat & Q, mat & Qd, mat & Qdd) const {

	// first create hitting polynomials
	vec7 a2, a3;
	vec7 b2, b3; // for returning
	double T = optim_params.t;
	vec7 qf = optim_params.q;
	vec7 qfdot = optim_params.qdot;
	vec7 qnow = qact.q;
	vec7 qdnow = qact.qd;
	a3 = 2.0 * (qnow - qf) / pow(T,3) + (qfdot + qdnow) / pow(T,2);
	a2 = 3.0 * (qf - qnow) / pow(T,2) - (qfdot + qdnow) / T;
	b3 = 2.0 * (qf - q_rest_des) / pow(T,3) + (qfdot) / pow(T,2);
	b2 = 3.0 * (q_rest_des - qf) / pow(T,2) - (qfdot) / T;

	static const double dt = 0.002;
	int N_hit = (int)T/dt;
	vec times_hit = linspace<vec>(dt,T,N_hit);
	int N_return = (int)time2return/dt;
	vec times_ret = linspace<vec>(dt,time2return,N_return);

	mat Q_hit, Qd_hit, Qdd_hit, Q_ret, Qd_ret, Qdd_ret;
	Q_hit = Qd_hit = Qdd_hit = zeros<mat>(NDOF,N_hit);
	Q_ret = Qd_ret = Qdd_ret = zeros<mat>(NDOF,N_return);

	gen_3rd_poly(times_hit,a3,a2,qdnow,qnow,Q_hit,Qd_hit,Qdd_hit);
	gen_3rd_poly(times_ret,a3,a2,qfdot,qf,Q_ret,Qd_ret,Qdd_ret);
	Q = join_horiz(Q_hit,Q_ret);
	Qd = join_horiz(Qd_hit,Qd_ret);
	Qdd = join_horiz(Qdd_hit,Qdd_ret);
}

/*
 * Function that calculates a racket strategy : positions, velocities and racket normal
 * for each point on the predicted ball trajectory (ballMat)
 * to return it a desired point (ballLand) at a desired time (landTime)
 *
 */
void Player::calc_racket_strategy(const mat & balls_predicted) {

	int N = balls_predicted.n_cols;
	mat balls_out_vel = zeros<mat>(3,N);
	mat racket_des_normal = zeros<mat>(3,N);
	mat racket_des_vel = zeros<mat>(3,N);
	calc_des_ball_out_vel(balls_predicted,balls_out_vel);
	calc_des_racket_normal(balls_predicted.rows(DX,DZ),balls_out_vel,racket_des_normal);
	calc_des_racket_vel(balls_predicted.rows(DX,DZ),balls_out_vel,racket_des_normal,racket_des_vel);

	racket_params.pos = balls_predicted;
	racket_params.normal = racket_des_normal;
	racket_params.vel = racket_des_vel;
}

/*
 * Calculate desired racket normal assuming mirror law
 */
void Player::calc_des_racket_normal(const mat & v_in, const mat & v_out, mat & normal) const {

	normal = v_out - v_in;
	// normalize
	normalise(normal);
}

/*
 *
 *  Computes the desired outgoing velocity of the ball after contact
 *	to hit the goal on a desired landing position on the
 *	opponents court
 *
 *
 */
void Player::calc_des_ball_out_vel(const mat & balls_predicted, mat & balls_out_vel) const {

	static double z_table = floor_level - table_height + ball_radius;

	// elementwise division
	balls_out_vel.row(X) = (ball_land_des(X) - balls_predicted.row(X)) / time_land_des;
	balls_out_vel.row(Y) = (ball_land_des(Y) - balls_predicted.row(Y)) / time_land_des;
	balls_out_vel.row(Z) = (z_table - balls_predicted.row(Z) +
			                0.5 * gravity * pow(time_land_des,2)) / time_land_des;

	//TODO: consider air drag, hack for now
	balls_out_vel.row(X) *= 1.1;
	balls_out_vel.row(Y) *= 1.1;
	balls_out_vel.row(Z) *= 1.2;
}

/*
 * Initialize an Extended Kalman Filter
 * useful for passing to Player constructor
 *
 */
EKF init_filter() {

	double std = 0.001;
	mat C = eye<mat>(3,6);
	mat66 Q = zeros<mat>(6,6);
	mat33 R = std * eye<mat>(3,3);
	return EKF(calc_next_ball,C,Q,R);
}

/*
 * Generate matrix of joint angles, velocities and accelerations
 */
void gen_3rd_poly(const vec & times, const vec7 & a3, const vec7 & a2, const vec7 & a1, const vec7 & a0,
		     mat & Q, mat & Qd, mat & Qdd) {

	// IN MATLAB:
	//	qStrike(m,:) = a(1)*t.^3 + a(2)*t.^2 + a(3)*t + a(4);
	//	qdStrike(m,:) = 3*a(1)*t.^2 + 2*a(2)*t + a(3);
	//	qddStrike(m,:) = 6*a(1)*t + 2*a(2);

	for(int i = 0; i < NDOF; i++) {
		Q.row(i) = a3(i) * pow(times,3) + a2(i) * pow(times,2) + a1(i) * times + a0;
		Qd.row(i) = 3*a3(i) * pow(times,2) + 2*a2(i) * times + a1(i);
		Qdd.row(i) = 6*a3(i) * times + 2*a2(i);
	}
}

/*
 * Calculate desired racket velocity given ball incoming and
 * outgoing velocities
 * Assuming a mirror law
 * Assumes no desired spin, i.e. racket velocity along the racket will be set to zero
 *
 * Output is the last parameter: racketVel
 *
 */
void calc_des_racket_vel(const mat & vel_ball_in, const mat & vel_ball_out,
		                 const mat & racket_normal, mat & racket_vel) {

	int N = vel_ball_in.n_elem;
	for (int i = 0; i < N; i++) {
		racket_vel.col(i) = dot((vel_ball_out.col(i) + CRR * vel_ball_in.col(i) / (1 + CRR)),
								racket_normal.col(i)) * racket_normal.col(i);
	}
}

/*
 * Check if there is a bounce among the predicted ball states
 * If it is exactly one ball will be legal!
 *
 * TODO: no need to check after ball passes table
 *
 */
bool check_legal_ball(const mat & balls_predicted) {

	int num_bounces = 0;
	int N = balls_predicted.n_elem;

	// if sign of z-velocity changes then the ball bounces
	for (int i = 0; i < N-1; i++) {
		if (balls_predicted(Z,i) < 0 && balls_predicted(Z,i+1) > 0) {
			num_bounces++;
		}
	}

	// multiple bounces are predicted
	if (num_bounces > 1) {
		cout << "Multiple bounces predicted. Not moving" << endl;
		return false;
	}
	// no bounce is predicted
	if (num_bounces == 0) {
		cout << "No bounce predicted. Not moving\n" << endl;
		return false;
	}

	return true;
}

/*
 * Checks to see if the observation is new (updated)
 *
 * The blobs need to be at least tol = 1e-3 apart from each other in distance
 *
 */
bool check_new_obs(const vec3 & obs) {

	static vec3 last_obs = zeros<vec>(3);
	static const double tol = 1e-3;

	if (norm(obs - last_obs) > tol) {
		last_obs = obs;
		return true;
	}
	return false;
}
