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
 * 2. Which method (1 = VHP, 2 = FP, 3 = DP) can land more balls?
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
#include "lookup.h"

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
 * @param flags Flags/options for player class, initialized with c++11 (see player.hpp)
 *              or through player.cfg file (see sl_interface)
 */
Player::Player(const vec7 & q0, EKF & filter_, player_flags & flags)
               : filter(filter_), pflags(flags) {

	ball_land_des(X) += pflags.ball_land_des_offset[X];
	ball_land_des(Y) = dist_to_table - 3*table_length/4 + pflags.ball_land_des_offset[Y];
	q_rest_des = q0;
	observations = zeros<mat>(3,pflags.min_obs); // for initializing filter
	times = zeros<vec>(pflags.min_obs); // for initializing filter
	//load_lookup_table(lookup_table);

	double lb[2*NDOF+1];
	double ub[2*NDOF+1];
	double SLACK = 0.02;
	double Tmax = 1.0;
	set_bounds(lb,ub,SLACK,Tmax);

	switch (pflags.alg) {
		case FOCUS: {
			opt = new FocusedOptim(q0,lb,ub);
			pred_params.Nmax = 1000;
			break; }
		case VHP: {
			opt = new HittingPlane(q0,lb,ub);
			break; }
		case DP: {
			opt = new DefensiveOptim(q0,lb,ub,true,true); // use lookup
			DefensiveOptim *dp = static_cast<DefensiveOptim*>(opt);
			dp->set_weights(pflags.weights);
			dp->set_velocity_multipliers(pflags.mult_vel);
			dp->set_penalty_loc(pflags.penalty_loc);
			pred_params.Nmax = 1000;
			break; }
		default:
			throw ("Algorithm is not recognized!\n");
	}
	opt->set_return_time(pflags.time2return);
	opt->set_verbose(pflags.verbosity > 1);
	opt->set_detach(pflags.detach);
}

/*
 *
 * Deconstructor for the Player class.
 * Frees the pointer to Optimization classes
 *
 */
Player::~Player() {

	delete opt;
}

/**
 * If filter is initialized returns true
 * @return
 */
bool Player::filter_is_initialized() const {
	return init_ball_state;
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

	using std::thread;
	using std::ref;
	int verb = pflags.verbosity;
	bool newball = check_new_obs(obs,1e-3);
	valid_obs = false;

	if (check_reset_filter(newball,verb,pflags.t_reset_thresh)) {
		filter = init_filter(pflags.var_model,pflags.var_noise,pflags.spin,pflags.out_reject_mult);
		num_obs = 0;
		init_ball_state = false;
		game_state = AWAITING;
		//cout << obs << endl;
		t_obs = 0.0; // t_cumulative
	}

	if (num_obs < pflags.min_obs && newball) {
		times(num_obs) = t_obs;
		observations.col(num_obs) = obs;
		num_obs++;
		if (num_obs == pflags.min_obs) {
			if (verb >= 1)
				cout << "Estimating initial ball state\n";
			thread t = thread(estimate_prior,ref(observations),ref(times),
					          ref(pflags.verbosity),ref(init_ball_state),ref(filter));
			if (pflags.detach)
				t.detach();
			else
				t.join();
			//estimate_prior(observations,times,pflags.verbosity > 2,filter);
			//cout << OBS << TIMES << filter.get_mean() << endl;
		}

	}
	else if (init_ball_state) { // comes here if there are enough balls to start filter
		filter.predict(DT,true); //true);
		if (newball) {
			valid_obs = true;
			if (pflags.outlier_detection)
				valid_obs = !filter.check_outlier(obs,verb > 2);
		}
		if (valid_obs) {
			filter.update(obs);
			//vec x = filter.get_mean();
			//mat P = (filter.get_covar());
			//cout << "OBS:" << obs.t() << "STATE:" << x.t() << "VAR:" << P.diag().t() << endl;
		}

	}
	t_obs += DT;
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

	switch (pflags.alg) {
		case FOCUS:
			optim_fp_param(qact);
			break;
		case VHP:
			optim_vhp_param(qact);
			break;
		case DP:
			optim_dp_param(qact);
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
 * @param ballstate Ball state (positions AND velocities).
 * @param qdes Desired joint positions, velocities, accelerations.
 */
void Player::cheat(const joint & qact, const vec6 & ballstate, joint & qdes) {

	// resetting legal ball detecting to AWAITING state
	if (ballstate(Y) < (dist_to_table - table_length) && ballstate(DY) > 2.0)
		game_state = AWAITING;
	filter.set_prior(ballstate,0.01*eye<mat>(6,6));

	switch (pflags.alg) {
		case FOCUS:
			optim_fp_param(qact);
			break;
		case VHP:
			optim_vhp_param(qact);
			break;
		case DP:
			optim_dp_param(qact);
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

	vec6 ball_pred;
	double time_pred;

	// if ball is fast enough and robot is not moving consider optimization
	if (check_update(qact)) {
		if (predict_hitting_point(pflags.VHPY,pflags.check_bounce,ball_pred,time_pred,filter,game_state)) { // ball is legal and reaches VHP
			calc_racket_strategy(ball_pred,ball_land_des,pflags.time_land_des,pred_params);
			HittingPlane *vhp = static_cast<HittingPlane*>(opt);
			vhp->set_des_params(&pred_params);
			vhp->fix_hitting_time(time_pred);
			vhp->update_init_state(qact);
			vhp->run();
		}
	}
}

/*
 * FOCUSED PLAYER
 *
 * Calculate the optimization parameters using an NLOPT nonlinear optimization algorithm
 * in another thread
 *
 * The optimized parameters are: qf, qf_dot, T
 * assuming T_return and q0 are fixed
 *
 */
void Player::optim_fp_param(const joint & qact) {

	mat balls_pred;

	// if ball is fast enough and robot is not moving consider optimization
	if (check_update(qact)) {
		predict_ball(2.0,balls_pred,filter);
		if (!pflags.check_bounce || check_legal_ball(filter.get_mean(),balls_pred,game_state)) { // ball is legal
			//lookup_soln(filter.get_mean(),1,qact);
			calc_racket_strategy(balls_pred,ball_land_des,pflags.time_land_des,pred_params);
			FocusedOptim *fp = static_cast<FocusedOptim*>(opt);
			fp->set_des_params(&pred_params);
			fp->update_init_state(qact);
			fp->run();
		}
		else {
			//cout << "Ball is not legal!\n";
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
void Player::optim_dp_param(const joint & qact) {

	mat balls_pred;

	// if ball is fast enough and robot is not moving consider optimization
	if (check_update(qact)) {
		predict_ball(2.0,balls_pred,filter);
		if (!pflags.check_bounce || check_legal_ball(filter.get_mean(),balls_pred,game_state)) { // ball is legal
			//lookup_soln(filter.get_mean(),1,qact);
			//calc_racket_strategy(balls_pred,ball_land_des,time_land_des,pred_params);
			pred_params.ball_pos = balls_pred.rows(X,Z);
			pred_params.ball_vel = balls_pred.rows(DX,DZ);
			pred_params.Nmax = balls_pred.n_cols;
			DefensiveOptim *dp = static_cast<DefensiveOptim*>(opt);
			dp->set_des_params(&pred_params);
			dp->update_init_state(qact);
			dp->run();
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

	static int firsttime = true;
	static int counter;
	static vec6 state_last = zeros<vec>(6);
	static wall_clock timer;
	vec6 state_est;
	bool update = false;
	racket robot_racket;
	bool activate, passed_lim, incoming, feasible = false;

	if (firsttime) {
		timer.tic();
		firsttime = false;
	}

	try {
		state_est = filter.get_mean();
		counter++;
		feasible = (state_est(DY) > 0.5) &&
				   (state_est(Y) > (dist_to_table - table_length + pflags.optim_offset));
		update = !opt->check_update() && !opt->check_running();
		// ball is incoming
		if (pflags.mpc) {// && t_poly > 0.0) {
			calc_racket_state(qact,robot_racket);
			activate = (!pflags.detach) ? counter % 5 == 0 :
					                        timer.toc() > (1.0/pflags.freq_mpc);
			passed_lim = state_est(Y) > robot_racket.pos(Y);
			incoming = state_est(Y) > state_last(Y);
			update = update && valid_obs && activate && feasible && !passed_lim && incoming;
		}
		else {
			update = update && (t_poly == 0.0) && feasible; // only once
		}
		state_last = state_est;
		if (update) {
			//cout << "Consider reacting to ball: " << state_est.t() << endl;
			//cout << num_updates++ << endl;
			timer.tic();
		}
	}
	catch (const std::exception & not_init_error) {
		update = false;
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

	using std::thread;
	using std::ref;
	// this should be only for MPC?
	if (opt->get_params(qact,poly)) {
		if (pflags.verbosity) {
			std::cout << "Launching/updating strike" << std::endl;
		}
		t_poly = DT;
		opt->set_moving(true);
	}

	// make sure we update after optim finished
	if (t_poly > 0.0) {
		if (!update_next_state(poly,q_rest_des,pflags.time2return,t_poly,qdes)) {
			opt->set_moving(false);
			// optimize to find a better resting state close to predicted balls
			if (pflags.optim_rest_posture)
				opt->run_qrest_optim(q_rest_des);
		}
	}

}

/**
 * @brief Method useful for testing performance of different players.
 *
 * Using many balls in simulation requires fast resetting
 * Setting a time threshold as a resetting condition won't work in this case.
 *
 */
void Player::reset_filter(double var_model, double var_noise) {

	filter = init_filter(var_model,var_noise,pflags.spin,pflags.out_reject_mult);
	init_ball_state = false;
	num_obs = 0;
	game_state = AWAITING;
	t_obs = 0.0;
}

/**
 * @brief Start moving pre-optim based on lookup if lookup flag is turned ON.
 *
 */
void Player::lookup_soln(const vec6 & ball_state, const int k, const joint & qact) {

	double time2return = pflags.time2return;
	if (t_poly == 0.0) {
		if (pflags.verbosity) {
			cout << "Starting movement based on lookup, k = 5\n"; // kNN parameter k = 5
		}
		vec::fixed<15> robot_params;
		vec6 ball_est = ball_state;
		//cout << "Init ball est:" << ball_params << endl;
		predict_till_net(ball_est);
		//cout << "Net ball est:" << ball_params << endl;
		knn(lookup_table,ball_est,k,robot_params);
		vec7 qf, qfdot;
		for (int i = 0; i < NDOF; i++) {
			qf(i) = robot_params(i);
			qfdot(i) = robot_params(i+NDOF);
		}
		double T = robot_params(2*NDOF);
		vec7 qnow = qact.q;
		vec7 qdnow = qact.qd;
		poly.a.col(0) = 2.0 * (qnow - qf) / pow(T,3) + (qfdot + qdnow) / pow(T,2);
		poly.a.col(1) = 3.0 * (qf - qnow) / pow(T,2) - (qfdot + 2.0*qdnow) / T;
		poly.a.col(2) = qdnow;
		poly.a.col(3) = qnow;
		//cout << "A = \n" << p.a << endl;
		poly.b.col(0) = 2.0 * (qf - q_rest_des) / pow(time2return,3) + (qfdot) / pow(time2return,2);
		poly.b.col(1) = 3.0 * (q_rest_des - qf) / pow(time2return,2) - (2.0*qfdot) / time2return;
		poly.b.col(2) = qfdot;
		poly.b.col(3) = qf;
		poly.time2hit = T;
		t_poly = DT;
	}
}

/**
 * @brief Predict hitting point on the Virtual Hitting Plane (VHP)
 *
 * If the ball is legal (legal detected bounce or legal predicted bounce)
 * and there is enough time (50 ms threshold) predict loc. on VHP.
 *
 * The location of the VHP is given as an argument (vhp-y-location vhpy)
 *
 */
bool predict_hitting_point(const double & vhpy, const bool & check_bounce,
		                   vec6 & ball_pred, double & time_pred,
		                   EKF & filter, game & game_state) {

	const double time_min = 0.05;
	mat balls_path;
	bool valid_hp = false;
	predict_ball(2.0,balls_path,filter);
	uvec vhp_index;
	unsigned idx;

	if (!check_bounce || check_legal_ball(filter.get_mean(),balls_path,game_state)) { // ball is legal
		vhp_index = find(balls_path.row(Y) >= vhpy, 1);
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

/**
 * @brief Predict ball with the models fed into the filter
 *
 * Number of prediction steps is given by Nmax in racket
 * parameters
 *
 */
void predict_ball(const double & time_pred, mat & balls_pred, EKF & filter) {

	//static wall_clock timer;
	//timer.tic();
	int N = (int)(time_pred/DT);
	balls_pred = filter.predict_path(DT,N);
	//cout << "Pred. ball time: " << 1000 * timer.toc() << " ms." << endl;
}

/**
 * @brief Compute desired racket pos,vel,normals and/or ball positions, vels.
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
	racket_params.ball_pos = balls_predicted.rows(X,Z);
	racket_params.ball_vel = balls_predicted.rows(DX,DZ);
	racket_params.racket_pos = balls_predicted.rows(X,Z);
	racket_params.racket_vel = racket_des_vel;
	racket_params.racket_normal = racket_des_normal;

	//cout << "Pred. racket time: " << 1000 * timer.toc() << " ms." << endl;
	return racket_params;
}

/**
 * @brief Compute spin-based desired racket pos,vel,normals and/or ball positions, vels.
 *
 * Outgoing ball velocity is computed using a boundary value problem
 * for the spin model, which is initialized using ballistic model.
 *
 * Function that calculates a racket strategy : positions, velocities and racket normal
 * for each point on the predicted ball trajectory (ballMat)
 * to return it a desired point (ballLand) at a desired time (landTime)
 *
 *
 */
optim_des calc_spin_racket_strategy(const mat & balls_predicted,
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
	racket_params.ball_pos = balls_predicted.rows(X,Z);
	racket_params.ball_vel = balls_predicted.rows(DX,DZ);
	racket_params.racket_pos = balls_predicted.rows(X,Z);
	racket_params.racket_vel = racket_des_vel;
	racket_params.racket_normal = racket_des_normal;

	//cout << "Pred. racket time: " << 1000 * timer.toc() << " ms." << endl;
	return racket_params;
}


/**
 *
 * @brief Checks for legal ball bounce
 * If an incoming ball has bounced before
 * it is declared ILLEGAL (legal_bounce as DATA MEMBER of player class)
 *
 * Bounce variable is static variable of estimate_ball_state method of player class
 * which is reset each time an incoming ball from ball gun is detected.
 *
 * TODO: also consider detecting HIT by robot racket
 * We can turn this off in configuration file
 *
 */
void check_legal_bounce(const vec6 & ball_est, game & game_state) {

	static double last_y_pos = 0.0;
	static double last_z_vel = 0.0;
	bool incoming = (ball_est(DY) > 0.0);
	bool on_opp_court = (ball_est(Y) < (dist_to_table - (table_length/2.0)));
	bool bounce = (last_z_vel < 0.0 && ball_est(DZ) > 0.0)
			       && (fabs(ball_est(Y) - last_y_pos) < 0.1);

	if (bounce && incoming) {
		// incoming ball has bounced
		if (game_state == LEGAL) {
			cout << "Ball bounced twice\n";
			game_state = ILLEGAL;
		}
		else if (game_state == AWAITING && on_opp_court) {
			cout << "Ball bounced on opponents court!\n";
			game_state = ILLEGAL;
		}
		else {
			cout << "Legal bounce occurred!" << endl;
			game_state = LEGAL;
		}
	}
	last_y_pos = ball_est(Y);
	last_z_vel = ball_est(DZ);
}

/**
 * @brief Check if the table tennis trial is LEGAL (hence motion planning can be started).
 *
 * If it is exactly one predicted bounce when in awaiting mode
 * (before an actual legal bounce was detected) trial will be legal!
 *
 * If ball has bounced legally bounce, then there should be no more bounces.
 *
 * TODO: no need to check after ball passes table
 * We can turn this check off in the configuration file
 *
 */
bool check_legal_ball(const vec6 & ball_est, const mat & balls_predicted, game & game_state) {

	int num_bounces = 0;
	int N = balls_predicted.n_cols;

	check_legal_bounce(ball_est, game_state);

	// if sign of z-velocity changes then the ball bounces
	for (int i = 0; i < N-1; i++) {
		if (balls_predicted(DZ,i) < 0.0 && balls_predicted(DZ,i+1) > 0.0) {
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

/**
 * Get players strategy (if exists)
 * @param pos_land_des
 * @param time_land_des
 */
void Player::get_strategy(vec2 & pos_land_des, double & time_land_des) {

	switch (pflags.alg) {
		case FOCUS:
			pos_land_des = this->ball_land_des;
			time_land_des = this->pflags.time_land_des;
			break;
		case VHP:
			pos_land_des = this->ball_land_des;
			time_land_des = this->pflags.time_land_des;
			break;
		case DP:
			throw ("DP does not have a fixed return strategy!\n");
			break;
		default:
			throw ("Algorithm is not recognized!\n");
	}
}
