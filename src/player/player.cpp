/*! \mainpage Optimal Control based Table Tennis Trajectory Generation
 *
 * \section intro_sec Introduction
 *
 * Documentation for the 'polyoptim' repository starts here.
 *
 *
 * \section install_sec Installation
 *
 * After pulling run 'make install'.
 * This will allow us to run the unit tests
 * where we can validate the results found in the paper.
 *
 *
 */


/**
 * @file player.cpp
 *
 * @brief Table Tennis player class and the functions it uses are stored here.
 *
 * Player class launches 3 different optimization algorithms
 * for striking trajectory generation.
 *
 *  Created on: Feb 8, 2017
 *  Author: okoc
 */

#include <armadillo>
#include <thread>
#include <string>
#include "stdlib.h"
#include "player.hpp"

#include "kinematics.h"
#include "kinematics.hpp"
#include "constants.h"
#include "kalman.h"
#include "tabletennis.h"
#include "utils.h"
#include "optim.h"

using namespace arma;

/**
 * @brief Initialize Table Tennis Player.
 *
 * Table Tennis Player class that can run 3 different trajectory
 * generation algorithms.
 * VHP and FP try to return the ball to the centre of the opponents court,
 * LP tries to just return the ball to the opponents court.
 *
 * TODO: remove callocs and replace with new! and add unique/shared pointers

 * @param q0 Initial joint positions.
 * @param filter_ Reference to an input filter (must be initialized in a separate line).
 * @param alg_ The algorithm for running trajectory optimization: VHP, FP, LP
 * @param mpc_ Flag for turning on/off model predictive control (i.e. corrections).
 * @param verbose_ Flag for changing verbosity level (0 = OFF, 1 = PLAYER, 2 = PLAYER + OPTIM).
 */
Player::Player(const vec7 & q0, EKF & filter_, algo alg_, bool mpc_, int verbose_, mode_operate mode_)
               : filter(filter_), alg(alg_), mpc(mpc_), verbose(verbose_), mode(mode_) {

	time_land_des = 0.8;
	time2return = 1.0;
	num_obs = 0;
	valid_obs = true;
	t_cum = 0.0;
	game_state = AWAITING;

	ball_land_des(X) = 0.0;
	ball_land_des(Y) = dist_to_table - 3*table_length/4;
	q_rest_des = q0;

	int N = 500;
	double** pos = my_matrix(0,NCART,0,N);
	double** vel = my_matrix(0,NCART,0,N);
	double** normal = my_matrix(0,NCART,0,N);
	racket_params = {pos, vel, normal, 0.002, N};

	double* qzerodot = (double*)calloc(NDOF,sizeof(double));
	double* qzerodot2 = (double*)calloc(NDOF,sizeof(double));
	double* qzero = (double*)calloc(NDOF, sizeof(double));
	double* qinit = (double*)calloc(NDOF, sizeof(double));
	double* qrest = (double*)calloc(NDOF, sizeof(double));
	double *lb = (double*)calloc(OPTIM_DIM,sizeof(double));
	double *ub = (double*)calloc(OPTIM_DIM,sizeof(double));
	double SLACK = 0.02;
	double Tmax = 1.0;
	set_bounds(lb,ub,SLACK,Tmax);

	for (int i = 0; i < NDOF; i++) {
		qinit[i] = qrest[i] = qzero[i] = q0(i);
	}

	optim_params = {qzero, qzerodot, 0.5, false, false};
	coparams = {qinit, qzerodot2, qrest, lb, ub, time2return, false, verbose > 1};

}

/*
 *
 * Deconstructor for the Player class.
 * Frees the memory using free() as in C-style since calloc() was called.
 *
 */
Player::~Player() {

	free(optim_params.qf);
	free(optim_params.qfdot);
	free(coparams.lb);
	free(coparams.ub);
	free(coparams.q0dot);
	free(coparams.qrest);
	free(coparams.q0);
	my_free_matrix(racket_params.normal,0,NCART,0,racket_params.Nmax);
	my_free_matrix(racket_params.pos,0,NCART,0,racket_params.Nmax);
	my_free_matrix(racket_params.vel,0,NCART,0,racket_params.Nmax);
}

/*
 * Filter the blob information with a Kalman Filter.
 * (Extended) KF is used both in simulation mode and for real robot.
 *
 * Checking for new ball that is at least 1 mm away from last observation
 * Checking for also outliers.
 * Resets if the ball suddenly appears on the opponent's court.
 *
 * Returns the valid ball flag.
 *
 * TODO: we're assuming that time elasped dt = DT = 0.002 seconds every time!
 *
 */
void Player::estimate_ball_state(const vec3 & obs) {

	// observation matrix
	static const int min_obs = 5;
	static mat OBS = zeros<mat>(3,min_obs);
	static vec TIMES = zeros<vec>(min_obs);
	bool newball = check_new_obs(obs,1e-3);
	valid_obs = false;

	if (check_reset_filter(newball,verbose,filter)) {
		num_obs = 0;
		game_state = AWAITING;
		t_cum = 0.0; // t_cumulative
	}

	if (num_obs < min_obs) {
		if (newball) {
			TIMES(num_obs) = t_cum;
			OBS.col(num_obs) = obs;
			num_obs++;
			if (num_obs == min_obs) {
				if (verbose >= 1)
					cout << "Estimating initial ball state\n";
				estimate_prior(OBS,TIMES,filter);
				//cout << OBS << TIMES << filter.get_mean() << endl;
			}
		}
	}
	else { // comes here if there are enough balls to start filter
		filter.predict(DT,true);
		if (newball) {
			valid_obs = true;
			if (mode == REAL_ROBOT)
				valid_obs = !filter.check_outlier(obs,verbose > 0);
		}
		if (valid_obs) {
			filter.update(obs);
			//cout << "Updating...\n"
			//     << "OBS\t" << obs.t() << "FILT\t" << filter.get_mean().t();
		}

	}
	t_cum += DT;
}

/**
 *
 * @brief Public interface for estimating ball state.
 *
 * This interface allows us to test/debug ball state estimation
 * (which is private).
 *
 * @param obs Ball position observation as a 3-vector.
 * @return Ball state as a 6-vector, if filter is not initialized,
 * returns the observation as positions and zeroes as velocities.
 */
vec6 Player::filt_ball_state(const vec3 & obs) {

	estimate_ball_state(obs);
	try {
		return filter.get_mean();
	}
	catch (const std::exception & exception) {
		return join_vert(obs,zeros<vec>(3));
	}
}

/**
 * @brief Play Table Tennis.
 *
 * Main function for playing Table Tennis. Calls one of three different
 * trajectory generation algorithms (depending on initialization) and
 * updates the desired joint states when the optimization threads have finished.
 *
 * @param qact Actual joint positions, velocities, accelerations.
 * @param ball_obs Ball observations (positions as 3-vector).
 * @param qdes Desired joint positions, velocities, accelerations.
 */
void Player::play(const joint & qact,
		           const vec3 & ball_obs,
				   joint & qdes) {

	estimate_ball_state(ball_obs);

	// initialize optimization and get the hitting parameters
	switch (alg) {
		case FIXED:
			optim_fixedp_param(qact);
			break;
		case VHP:
			optim_vhp_param(qact);
			break;
		case LAZY:
			optim_lazy_param(qact);
			break;
		default:
			throw ("Algorithm is not recognized!\n");
	}

	// generate movement or calculate next desired step
	calc_next_state(qact, qdes);

}

/**
 * @brief Cheat Table Tennis by getting the exact ball state in simulation.
 *
 * Similar to play() method, this method receives from the simulator the
 * exact ball states, which bypasses then the ball estimation method.
 * Useful for debugging the internal filter used.
 *
 * @param qact Actual joint positions, velocities, accelerations.
 * @param ball_obs Ball observations (positions as 3-vector).
 * @param qdes Desired joint positions, velocities, accelerations.
 */
void Player::cheat(const joint & qact, const vec6 & ballstate, joint & qdes) {

	// resetting legal ball detecting to AWAITING state
	if (ballstate(Y) < (dist_to_table - table_length) && ballstate(DY) > 2.0)
		game_state = AWAITING;
	filter.set_prior(ballstate,0.01*eye<mat>(6,6));

	switch (alg) {
		case FIXED:
			optim_fixedp_param(qact);
			break;
		case VHP:
			optim_vhp_param(qact);
			break;
		case LAZY:
			optim_lazy_param(qact);
			break;
		default:
			throw ("Algorithm is not recognized!\n");
	}

	// generate movement or calculate next desired step
	calc_next_state(qact, qdes);
}

/*
 * Calculate hitting parameters qf, qfdot
 * on the Virtual Hitting Plane (VHP) by running Inverse Kinematics
 *
 * The inverse kinematics routine runs an optimization to minimize
 * the distance to a rest posture
 *
 *
 */
void Player::optim_vhp_param(const joint & qact) {

	double time_pred;
	vec6 balls_pred;

	// if ball is fast enough and robot is not moving consider optimization
	if (check_update(qact)) {
		if (predict_hitting_point(balls_pred,time_pred)) { // ball is legal and reaches VHP
			calc_racket_strategy(balls_pred,ball_land_des,
					time_land_des,racket_params);
			for (int i = 0; i < NDOF; i++) {
				coparams.q0[i] = qact.q(i);
				coparams.q0dot[i] = qact.qd(i);
			}
			optim_params.T = time_pred;
			// run optimization in another thread
			std::thread t(&nlopt_vhp_run,
					&coparams,&racket_params,&optim_params);
			//t.join();
			t.detach();
		}
	}

}

/*
 * FIXED PLAYER
 *
 * Calculate the optimization parameters using an NLOPT nonlinear optimization algorithm
 * in another thread
 *
 * The optimized parameters are: qf, qf_dot, T
 * assuming T_return and q0 are fixed
 *
 */
void Player::optim_fixedp_param(const joint & qact) {

	static mat balls_pred;

	// if ball is fast enough and robot is not moving consider optimization
	if (check_update(qact)) {
		predict_ball(balls_pred);
		if (check_legal_ball(filter.get_mean(),balls_pred,game_state)) { // ball is legal
			calc_racket_strategy(balls_pred,ball_land_des,time_land_des,racket_params);
			for (int i = 0; i < NDOF; i++) {
				coparams.q0[i] = qact.q(i);
				coparams.q0dot[i] = qact.qd(i);
			}
			//cout << filter.get_mean().t() << endl;
			// run optimization in another thread
			std::thread t(&nlopt_optim_fixed_run,
					&coparams,&racket_params,&optim_params);
			//t.join();
			t.detach();
		}
	}
}

/*
 * LAZY PLAYER
 *
 * Calculate the optimization parameters using an NLOPT nonlinear optimization algorithm
 * in another thread
 *
 * The optimized parameters are: qf, qf_dot, T
 *
 * TODO: delete ballpred in destructor?
 *
 */
void Player::optim_lazy_param(const joint & qact) {

	static double** ballpred = my_matrix(0,2*NCART,0,racket_params.Nmax);
	static mat balls_pred;

	// if ball is fast enough and robot is not moving consider optimization
	if (check_update(qact)) {
		predict_ball(balls_pred);
		if (check_legal_ball(filter.get_mean(),balls_pred,game_state)) { // ball is legal
			calc_racket_strategy(balls_pred,ball_land_des,time_land_des,racket_params);
			for (int i = 0; i < NDOF; i++) {
				coparams.q0[i] = qact.q(i);
				coparams.q0dot[i] = qact.qd(i);
			}
			for (int i = 0; i < racket_params.Nmax; i++) {
				for (int j = 0; j < 2*NCART; j++) {
					ballpred[j][i] = balls_pred(j,i);
				}
			}
			// run optimization in another thread
			std::thread t(&nlopt_optim_lazy_run,
					ballpred,&coparams,&racket_params,&optim_params);
			//t.join();
			t.detach();
		}
	}


}

/*
 * Check MPC flag and update if possible
 *
 * IF MPC IS TURNED OFF
 * if ball is incoming and robot is not moving consider optimization
 *
 * IF MPC IS TURNED ON
 * then additionally consider (after running initial optimization)
 * relaunching optimization if ball is valid (new ball and not an outlier)
 * the frequency of updates is respected, and the ball has not passed the y-limit
 *
 */
bool Player::check_update(const joint & qact) const {

	vec6 state_est;
	bool update;
	static int counter;
	racket robot_racket;
	//static int num_updates;
	static const double FREQ_MPC = 10.0;
	static wall_clock timer;
	bool activate, passed_lim = false;

	try {
		state_est = filter.get_mean();
		counter++;
		update = !optim_params.update && !optim_params.running;
		// ball is incoming
		if (mpc && coparams.moving) {
			calc_racket_state(qact,robot_racket);
			activate = (timer.toc() > (1.0/FREQ_MPC));
			//activate = (counter % 20 == 0);
			passed_lim = state_est(Y) > robot_racket.pos(Y);
			update = update && valid_obs && activate && !passed_lim;
		}
		else {
			update = update && !coparams.moving && state_est(DY) > 0.0 && (state_est(Y) > (dist_to_table - table_length/2.0));
		}
	}
	catch (const std::exception & not_init_error) {
		update = false;
	}
	if (update) {
		//cout << num_updates++ << endl;
		timer.tic();
	}
	return update;
}

/*
 * Predict ball with the models fed into the filter
 *
 * Number of prediction steps is given by Nmax in racket
 * parameters
 *
 */
void Player::predict_ball(mat & balls_pred) const {

	//static wall_clock timer;
	//timer.tic();
	int N = racket_params.Nmax;
	//cout << filter.get_mean() << endl;
	balls_pred = filter.predict_path(racket_params.dt,N);
	//cout << timer.toc() << endl;
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

	static unsigned idx = 0;
	static mat Q_des, Qd_des, Qdd_des;

	// this should be only for MPC?
	if (optim_params.update) {
		if (verbose) {
			std::cout << "Launching/updating strike" << std::endl;
		}
		coparams.moving = true;
		optim_params.update = false;
		// call polynomial generation
		generate_strike(optim_params,qact,q_rest_des,time2return,Q_des,Qd_des,Qdd_des);
		idx = 0;
	}

	// make sure we update after optim finished
	if (coparams.moving) {
		qdes.q = Q_des.col(idx);
		qdes.qd = Qd_des.col(idx);
		qdes.qdd = Qdd_des.col(idx);
		idx++;
		if (idx == Q_des.n_cols) {
			// hitting process will finish
			coparams.moving = false;
			qdes.q = q_rest_des;
			qdes.qd = zeros<vec>(NDOF);
			qdes.qdd = zeros<vec>(NDOF);
			idx = 0;
		}
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
bool Player::predict_hitting_point(vec6 & ball_pred, double & time_pred) {

	static const double time_min = 0.05;
	static mat balls_path;
	bool valid_hp = false;
	predict_ball(balls_path);
	uvec vhp_index;
	unsigned idx;

	if (check_legal_ball(filter.get_mean(),balls_path,game_state)) {
		vhp_index = find(balls_path.row(Y) >= VHPY, 1);
		if (vhp_index.n_elem == 1) {
			idx = as_scalar(vhp_index);
			ball_pred = balls_path.col(idx);
			time_pred = racket_params.dt * (idx + 1);
			if (time_pred > time_min)
				valid_hp = true;
		}
	}

	return valid_hp;
}

/**
 * @brief Method useful for testing performance of different players.
 *
 * Using many balls in simulation requires fast resetting
 * Setting a time threshold as a resetting condition won't work in this case.
 *
 */
void Player::reset_filter(double std_model, double std_noise) {

	filter = init_filter(std_model,std_noise);
	num_obs = 0;
	game_state = AWAITING;
	t_cum = 0.0;
}

/*
 * Function that calculates a racket strategy : positions, velocities and racket normal
 * for each point on the predicted ball trajectory (ballMat)
 * to return it a desired point (ballLand) at a desired time (landTime)
 *
 */
racketdes calc_racket_strategy(const mat & balls_predicted,
		                       const vec2 & ball_land_des, const double time_land_des,
							   racketdes & racket_params) {

	//static wall_clock timer;
	//timer.tic();
	static TableTennis tennis = TableTennis(false,false);

	int N = balls_predicted.n_cols;
	mat balls_out_vel = zeros<mat>(3,N);
	mat racket_des_normal = zeros<mat>(3,N);
	mat racket_des_vel = zeros<mat>(3,N);
	tennis.calc_des_ball_out_vel(ball_land_des,time_land_des,balls_predicted,balls_out_vel);
	tennis.calc_des_racket_normal(balls_predicted.rows(DX,DZ),balls_out_vel,racket_des_normal);
	tennis.calc_des_racket_vel(balls_predicted.rows(DX,DZ),balls_out_vel,racket_des_normal,racket_des_vel);

	// place racket centre on the predicted ball

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < NCART; j++) {
			racket_params.pos[j][i] = balls_predicted(j,i);
			racket_params.vel[j][i] = racket_des_vel(j,i);
			racket_params.normal[j][i] = racket_des_normal(j,i);
		}
	}
	//cout << timer.toc() << endl;
	return racket_params;
}

/*
 * Create batch hitting and returning joint state 3rd degree polynomials
 *
 */
void generate_strike(const optim & params, const joint & qact,
		             const vec7 & q_rest_des, const double time2return,
		            mat & Q, mat & Qd, mat & Qdd) {

	// first create hitting polynomials
	vec7 a2, a3;
	vec7 b2, b3; // for returning
	double T = params.T;
	vec7 qf(params.qf);
	vec7 qfdot(params.qfdot);
	vec7 qnow = qact.q;
	vec7 qdnow = qact.qd;
	a3 = 2.0 * (qnow - qf) / pow(T,3) + (qfdot + qdnow) / pow(T,2);
	a2 = 3.0 * (qf - qnow) / pow(T,2) - (qfdot + 2.0*qdnow) / T;
	b3 = 2.0 * (qf - q_rest_des) / pow(time2return,3) + (qfdot) / pow(time2return,2);
	b2 = 3.0 * (q_rest_des - qf) / pow(time2return,2) - (2.0*qfdot) / time2return;

	int N_hit = T/DT;
	rowvec times_hit = linspace<rowvec>(DT,T,N_hit);
	int N_return = time2return/DT;
	rowvec times_ret = linspace<rowvec>(DT,time2return,N_return);

	mat Q_hit, Qd_hit, Qdd_hit, Q_ret, Qd_ret, Qdd_ret;
	Q_hit = Qd_hit = Qdd_hit = zeros<mat>(NDOF,N_hit);
	Q_ret = Qd_ret = Qdd_ret = zeros<mat>(NDOF,N_return);

	gen_3rd_poly(times_hit,a3,a2,qdnow,qnow,Q_hit,Qd_hit,Qdd_hit);
	gen_3rd_poly(times_ret,b3,b2,qfdot,qf,Q_ret,Qd_ret,Qdd_ret);
	Q = join_horiz(Q_hit,Q_ret);
	Qd = join_horiz(Qd_hit,Qd_ret);
	Qdd = join_horiz(Qdd_hit,Qdd_ret);
}

/*
 * Initialize an Extended Kalman Filter
 * useful for passing to Player constructor
 *
 */
EKF init_filter(double std_model, double std_noise) {

	mat C = eye<mat>(3,6);
	mat66 Q = std_model * eye<mat>(6,6);
	mat33 R = std_noise * eye<mat>(3,3);
	return EKF(calc_next_ball,C,Q,R);
}

/*
 * Generate matrix of joint angles, velocities and accelerations
 */
void gen_3rd_poly(const rowvec & times, const vec7 & a3, const vec7 & a2, const vec7 & a1, const vec7 & a0,
		     mat & Q, mat & Qd, mat & Qdd) {

	// IN MATLAB:
	//	qStrike(m,:) = a(1)*t.^3 + a(2)*t.^2 + a(3)*t + a(4);
	//	qdStrike(m,:) = 3*a(1)*t.^2 + 2*a(2)*t + a(3);
	//	qddStrike(m,:) = 6*a(1)*t + 2*a(2);

	for(int i = 0; i < NDOF; i++) {
		Q.row(i) = a3(i) * pow(times,3) + a2(i) * pow(times,2) + a1(i) * times + a0(i);
		Qd.row(i) = 3*a3(i) * pow(times,2) + 2*a2(i) * times + a1(i);
		Qdd.row(i) = 6*a3(i) * times + 2*a2(i);
	}
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
void check_legal_bounce(const vec6 & ball_est, game & game_state) {


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
bool check_legal_ball(const vec6 & ball_est, const mat & balls_predicted, game & game_state) {

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
 * Checks to see if the observation is new (updated)
 *
 * The blobs need to be at least tol apart from each other in distance
 *
 */
bool check_new_obs(const vec3 & obs, double tol) {

	static vec3 last_obs = zeros<vec>(3);

	if (norm(obs - last_obs) > tol) {
		last_obs = obs;
		return true;
	}
	return false;
}

/*
 * Least squares to estimate prior given
 * matrix of observations Y of column length N, each row
 * corresponding to one
 *
 * If observations arrive as column vectors then we take
 * transpose of it.
 *
 * Velocity estimation is biased, we multiply velocities by 1.1
 * since they often underestimate actual velocities.
 *
 *
 */
void estimate_prior(const mat & observations,
		            const vec & times,
					EKF & filter) {

	vec6 x; mat66 P;
	int num_samples = times.n_elem;
	mat M = zeros<mat>(num_samples,3);
	vec3 vel_multiplier = {1.1, 1.1, 1.1};

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
	//P *= 100.0;
	filter.set_prior(x,P);
	filter.update(observations.col(0));

	double dt;
	for (unsigned i = 1; i < times.n_elem; i++) {
		dt = times(i) - times(i-1);
		filter.predict(dt,true);
		filter.update(observations.col(i));
	}
	//x = filter.get_mean();
	//x(span(DX,DZ)) = x(span(DX,DZ)) % vel_multiplier;
	//filter.set_prior(x,P);
}

/*
 * Check to see if we want to reset the filter.
 *
 * Basically if a new ball appears 300 ms later than the last new ball
 * we reset the filter.
 *
 * TODO: check if it works on real ball data!
 *
 */
bool check_reset_filter(const bool newball, const int verbose, EKF & filter) {

	bool reset = false;
	static int reset_cnt = 0;
	static double threshold = 0.3; // 300 miliseconds
	static bool firsttime = true;
	static wall_clock timer;

	if (firsttime) {
		firsttime = false;
		timer.tic();
	}

	if (newball) {
		if (timer.toc() > threshold) {
			reset = true;
			if (verbose > 0) {
				std::cout << "Resetting filter! Count: " << ++reset_cnt << std::endl;
			}
			filter = init_filter();
		}
		timer.tic();
	}

	return reset;
}

/*
 * Set upper and lower bounds on the optimization.
 * First loads the joint limits and then
 */
void set_bounds(double *lb, double *ub, double SLACK, double Tmax) {

	read_joint_limits(lb,ub);
	// lower bounds and upper bounds for qf are the joint limits
	for (int i = 0; i < NDOF; i++) {
		ub[i] -= SLACK;
		lb[i] += SLACK;
		ub[i+NDOF] = MAX_VEL;
		lb[i+NDOF] = -MAX_VEL;
	}
	// constraints on final time
	ub[2*NDOF] = Tmax;
	lb[2*NDOF] = 0.0;
}
