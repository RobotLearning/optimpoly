/*! \mainpage Optimal Control based Table Tennis Trajectory Generation
 *
 * \section intro_sec Introduction
 *
 * Welcome to Table Tennis!
 *
 * Player class is the orchestrator for Table Tennis
 * which can call three different trajectory optimizers for table tennis.
 * These are 3rd order striking and returning polynomials computed
 * differently for each method.
 *
 * We provide three modes for testing/playing table tennis:
 * 1. Unit Tests, here the performances of three different players are compared.
 * 2. SL, here the simulation is real-time so threads are detached.
 * 3. Real-Robot, here filtering and trajectory corrections are more robust,
 *                outlier detection is also considered for instance.
 *
 * For the modes (2) and (3) configuration file "player.cfg" sets the modes
 * for the player class, these can be changed online (assuming robot has no task).
 *
 *
 * \section install_sec Installation
 *
 * After pulling run 'make install'.
 * This will allow us to run the unit tests
 * where we can validate the results found in the paper.
 * Make sure to type 'make test' and run ./unit_tests.o
 *
 * \section test_sec Unit Tests
 *
 * The unit tests, use boost testing framework and do not only
 * consider 'unit' tests for each method, but also general scripts for various
 * test scenarios. For instance
 *
 * 1. Does the ball land on the other side?
 * 2. Which method (1 = VHP, 2 = FP, 3 = LAZY) can land more balls?
 * 3. Is correction of trajectories useful when you have observation noise
 *    or  ball prediction error?
 * 4. Are the filters stable? Can we estimate ball state better as we get more data?
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
 * @param q0 Initial joint positions.
 * @param filter_ Reference to an input filter (must be initialized in a separate line).
 * @param alg_ The algorithm for running trajectory optimization: VHP, FP, LP
 * @param mpc_ Flag for turning on/off model predictive control (i.e. corrections).
 * @param verbose_ Flag for changing verbosity level (0 = OFF, 1 = PLAYER, 2 = PLAYER + OPTIM).
 * @param mode_ Mode of running player (0 = SIM FOR UNIT TESTS, 1 = SL SIM, 2 = REAL ROBOT!)
 */
Player::Player(const vec7 & q0, EKF & filter_, algo alg_, bool mpc_, int verbose_, mode_operate mode_)
               : filter(filter_), alg(alg_), mpc(mpc_), verbose(verbose_), mode(mode_) {

	time_land_des = 0.8;
	time2return = 1.0;
	num_obs = 0;
	valid_obs = true;
	moving = false;
	t_cum = 0.0;
	observations = zeros<mat>(3,min_obs);
	times = zeros<vec>(min_obs);
	game_state = AWAITING;

	ball_land_des(X) = 0.0;
	ball_land_des(Y) = dist_to_table - 3*table_length/4;
	q_rest_des = q0;
	opt_params.Nmax = 500;

	double lb[2*NDOF+1];
	double ub[2*NDOF+1];
	double qrest[NDOF];
	double SLACK = 0.02;
	double Tmax = 1.0;
	set_bounds(lb,ub,SLACK,Tmax);

	for (int i = 0; i < NDOF; i++) {
		qrest[i] = q0(i);
	}
	switch (alg) {
		case FIXED:
			opt = new FocusedOptim(qrest,lb,ub);
			break;
		case VHP:
			opt = new HittingPlane(qrest,lb,ub);
			break;
		case LAZY:
			opt = new LazyOptim(qrest,lb,ub);
			break;
		default:
			throw ("Algorithm is not recognized!\n");
	}
	opt->set_verbose(verbose > 1);
	opt->set_detach(mode != TEST_SIM);
}

/*
 *
 * Deconstructor for the Player class.
 * Frees the memory using free() as in C-style since calloc() was called.
 *
 */
Player::~Player() {

	delete opt;
}

/*
 * Filter the blob information with a Kalman Filter.
 * (Extended) KF is used both in simulation mode and for real robot.
 *
 * Checking for new ball that is at least 1 mm away from last observation
 * Checking for also outliers.
 * Resets if the ball suddenly appears on the opponent's court.
 *
 * Ball is valid if ball is a new ball and (in real robot mode)
 * it is not an outlier!
 *
 * Note: we're assuming that time elasped dt = DT = 0.002 seconds every time!
 *
 */
void Player::estimate_ball_state(const vec3 & obs) {

	bool newball = check_new_obs(obs,1e-3);
	valid_obs = false;

	if (check_reset_filter(newball,verbose)) {
		filter = init_filter(0.03,0.001,mode == REAL_ROBOT);
		num_obs = 0;
		game_state = AWAITING;
		t_cum = 0.0; // t_cumulative
	}

	if (num_obs < min_obs) {
		if (newball) {
			times(num_obs) = t_cum;
			observations.col(num_obs) = obs;
			num_obs++;
			if (num_obs == min_obs) {
				if (verbose >= 1)
					cout << "Estimating initial ball state\n";
				estimate_prior(observations,times,filter,mode == REAL_ROBOT);
				//cout << OBS << TIMES << filter.get_mean() << endl;
			}
		}
	}
	else { // comes here if there are enough balls to start filter
		filter.predict(DT,true); //true);
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
void Player::play(const joint & qact,const vec3 & ball_obs, joint & qdes) {

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
	vec6 ball_pred;
	double q0[NDOF], q0dot[NDOF];

	// if ball is fast enough and robot is not moving consider optimization
	if (check_update(qact)) {
		if (predict_hitting_point(ball_pred,time_pred,filter,game_state)) { // ball is legal and reaches VHP
			calc_racket_strategy(ball_pred,ball_land_des,time_land_des,opt_params);
			for (int i = 0; i < NDOF; i++) {
				q0[i] = qact.q(i);
				q0dot[i] = qact.qd(i);
			}
			opt->set_des_params(&opt_params);
			opt->update_init_state(q0,q0dot,time_pred);
			opt->run();
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

	mat balls_pred;
	double q0[NDOF], q0dot[NDOF];

	// if ball is fast enough and robot is not moving consider optimization
	if (check_update(qact)) {
		predict_ball(1.0,balls_pred,filter);
		if (check_legal_ball(filter.get_mean(),balls_pred,game_state)) { // ball is legal
			calc_racket_strategy(balls_pred,ball_land_des,time_land_des,opt_params);
			for (int i = 0; i < NDOF; i++) {
				q0[i] = qact.q(i);
				q0dot[i] = qact.qd(i);
			}
			opt->set_des_params(&opt_params);
			opt->update_init_state(q0,q0dot,0.5);
			opt->run();
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
 *
 */
void Player::optim_lazy_param(const joint & qact) {

	mat balls_pred;
	double q0[NDOF], q0dot[NDOF];

	// if ball is fast enough and robot is not moving consider optimization
	if (check_update(qact)) {
		predict_ball(1.0,balls_pred,filter);
		if (check_legal_ball(filter.get_mean(),balls_pred,game_state)) { // ball is legal
			opt_params.ball_pos = balls_pred.rows(X,Z);
			opt_params.ball_vel = balls_pred.rows(DX,DZ);
			for (int i = 0; i < NDOF; i++) {
				q0[i] = qact.q(i);
				q0dot[i] = qact.qd(i);
			}
			opt->set_des_params(&opt_params);
			opt->update_init_state(q0,q0dot,0.5);
			opt->run();
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

	static vec6 state_last = -10 * ones<vec>(6);
	static wall_clock timer;
	const double FREQ_MPC = (mode == REAL_ROBOT) ? 10.0 : 40.0;
	vec6 state_est;
	bool update;
	static int counter;
	racket robot_racket;
	bool activate, passed_lim, incoming = false;

	try {
		state_est = filter.get_mean();
		counter++;
		update = !opt->check_update() && !opt->check_running();
		// ball is incoming
		if (mpc && moving) {
			calc_racket_state(qact,robot_racket);
			activate = (mode == TEST_SIM) ? counter % 5 == 0 :
					                        timer.toc() > (1.0/FREQ_MPC);
			passed_lim = state_est(Y) > robot_racket.pos(Y);
			incoming = state_est(Y) > state_last(Y);
			update = update && valid_obs && activate
					&& !passed_lim && incoming;
		}
		else {
			update = update && !moving
					&& (state_est(Y) > (dist_to_table - table_length/2.0))
					&& (state_est(DY) > 2.0);
		}
	}
	catch (const std::exception & not_init_error) {
		update = false;
	}
	if (update) {
		//cout << num_updates++ << endl;
		timer.tic();
		state_last = state_est;
	}
	return update;
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

	// this should be only for MPC?
	if (opt->get_params(qact,poly)) {
		if (verbose) {
			std::cout << "Launching/updating strike" << std::endl;
		}
		moving = true;
		opt->set_moving(true);
	}

	// make sure we update after optim finished
	if (moving) {
		update_next_state(poly,q_rest_des,time2return,qdes);
	}

}

/**
 * @brief Method useful for testing performance of different players.
 *
 * Using many balls in simulation requires fast resetting
 * Setting a time threshold as a resetting condition won't work in this case.
 *
 */
void Player::reset_filter(double std_model, double std_noise) {

	filter = init_filter(std_model,std_noise,mode == REAL_ROBOT);
	num_obs = 0;
	game_state = AWAITING;
	t_cum = 0.0;
}

/*
 * Predict hitting point on the Virtual Hitting Plane (VHP)
 * if the ball is legal (legal detected bounce or legal predicted bounce)
 * and there is enough time (50 ms threshold)
 *
 * The location of the VHP is defined as a constant (constants.h)
 *
 */
bool predict_hitting_point(vec6 & ball_pred, double & time_pred,
		                   EKF & filter, game & game_state) {

	const double time_min = 0.05;
	mat balls_path;
	bool valid_hp = false;
	predict_ball(1.0,balls_path,filter);
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
void predict_ball(const double & time_pred, mat & balls_pred, EKF & filter) {

	int N = (int)(time_pred/DT);
	balls_pred = filter.predict_path(DT,N);
}

/*
 * Function that calculates a racket strategy : positions, velocities and racket normal
 * for each point on the predicted ball trajectory (ballMat)
 * to return it a desired point (ballLand) at a desired time (landTime)
 *
 *
 */
optim_des calc_racket_strategy(const mat & balls_predicted,
		                       const vec2 & ball_land_des, const double time_land_des,
							   optim_des & racket_params) {

	//static wall_clock timer;
	//timer.tic();
	TableTennis tennis = TableTennis(false,false);

	int N = balls_predicted.n_cols;
	mat balls_out_vel = zeros<mat>(3,N);
	mat racket_des_normal = zeros<mat>(3,N);
	mat racket_des_vel = zeros<mat>(3,N);
	tennis.calc_des_ball_out_vel(ball_land_des,time_land_des,balls_predicted,balls_out_vel);
	tennis.calc_des_racket_normal(balls_predicted.rows(DX,DZ),balls_out_vel,racket_des_normal);
	tennis.calc_des_racket_vel(balls_predicted.rows(DX,DZ),balls_out_vel,racket_des_normal,racket_des_vel);

	// place racket centre on the predicted ball
	racket_params.racket_pos = balls_predicted.rows(X,Z);
	racket_params.racket_vel = racket_des_vel;
	racket_params.racket_normal = racket_des_normal;

	//cout << timer.toc() << endl;
	return racket_params;
}

/*
 * Update state using hitting and returning joint state 3rd degree polynomials
 *
 */
void update_next_state(const spline_params & poly,
		           const vec7 & q_rest_des,
				   const double time2return, joint & qdes) {

	static double t = DT;
	mat a,b;
	double tbar;

	if (t <= poly.time2hit) {
		a = poly.a;
		qdes.q = a.col(0)*t*t*t + a.col(1)*t*t + a.col(2)*t + a.col(3);
		qdes.qd = 3*a.col(0)*t*t + 2*a.col(1)*t + a.col(2);
		qdes.qdd = 6*a.col(0)*t + 2*a.col(1);
	}
	else if (t <= poly.time2hit + time2return) {
		b = poly.b;
		tbar = t - poly.time2hit;
		qdes.q = b.col(0)*tbar*tbar*tbar + b.col(1)*tbar*tbar + b.col(2)*tbar + b.col(3);
		qdes.qd = 3*b.col(0)*tbar*tbar + 2*b.col(1)*tbar + b.col(2);
		qdes.qdd = 6*b.col(0)*tbar + 2*b.col(1);
	}
	else {
		t = 0.0; // hitting finished
		qdes.q = q_rest_des;
		qdes.qd = zeros<vec>(NDOF);
		qdes.qdd = zeros<vec>(NDOF);
	}
	t += DT;
}

/*
 * Initialize an Extended Kalman Filter
 * useful for passing to Player constructor
 *
 */
EKF init_filter(double std_model, double std_noise, bool spin) {

	mat C = eye<mat>(3,6);
	mat66 Q = std_model * eye<mat>(6,6);
	mat33 R = std_noise * eye<mat>(3,3);
	if (spin)
		return EKF(calc_spin_ball,C,Q,R);
	else
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
 * Velocity estimation is biased, we multiply velocities by 0.8
 * since LSE often overestimates the model with spin.
 *
 *
 */
void estimate_prior(const mat & observations,
		            const vec & times,
					EKF & filter,
					bool real_robot) {

	vec6 x; mat66 P;
	int num_samples = times.n_elem;
	mat M = zeros<mat>(num_samples,3);
	vec times_z = times - times(0); // times zeroed

	// and create the data matrix
	for (int i = 0; i < num_samples; i++) {
		M(i,0) = 1.0;
		M(i,1) = times_z(i);
		M(i,2) = times_z(i) * times_z(i);
	}
	// solving for the parameters
	//cout << "Data matrix:" << endl << M << endl;
	mat Beta = solve(M,observations.t());
	//cout << "Parameters:" << endl << Beta << endl;
	x = join_horiz(Beta.row(0),Beta.row(1)).t(); //vectorise(Beta.rows(0,1));
	P.eye(6,6);
	if (real_robot)	{
		x(span(DX,DZ)) *= 0.8;
		P *= 1e3;
	}
	//cout << "Data:\n" << observations << endl;
	//cout << "Initial est:\n" << x << endl;
	filter.set_prior(x,P);
	filter.update(observations.col(0));

	double dt;
	for (unsigned i = 1; i < times.n_elem; i++) {
		dt = times_z(i) - times_z(i-1);
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
 */
bool check_reset_filter(const bool newball, const int verbose) {

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
