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
using namespace optim;

namespace player {

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

Player::~Player() {

	delete opt;
}

bool Player::filter_is_initialized() const {
	return init_ball_state;
}

void Player::estimate_ball_state(const vec3 & obs, const double & dt) {

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
		filter.predict(dt,true); //true);
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
	t_obs += dt;
}

vec6 Player::filt_ball_state(const vec3 & obs) {

	estimate_ball_state(obs);
	try {
		return filter.get_mean();
	}
	catch (const std::exception & exception) {
		return join_vert(obs,zeros<vec>(3));
	}
}

void Player::play(const joint & qact,
                  const vec3 & ball_obs,
                  joint & qdes) {

	estimate_ball_state(ball_obs,const_tt::DT);

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

void Player::cheat(const joint & qact,
                    const vec6 & ballstate,
                    joint & qdes) {

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

void Player::calc_next_state(const joint & qact, joint & qdes) {

	using std::thread;
	using std::ref;
	// this should be only for MPC?
	if (opt->get_params(qact,poly)) {
		if (pflags.verbosity) {
			std::cout << "Launching/updating strike" << std::endl;
		}
		t_poly = const_tt::DT;
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

void Player::reset_filter(double var_model, double var_noise) {

	filter = init_filter(var_model,var_noise,pflags.spin,pflags.out_reject_mult);
	init_ball_state = false;
	num_obs = 0;
	game_state = AWAITING;
	t_obs = 0.0;
}

void Player::lookup_soln(const vec6 & ball_state,
                         const int k,
                         const joint & qact) {

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
		t_poly = const_tt::DT;
	}
}

bool predict_hitting_point(const double & vhpy,
                           const bool & check_bounce,
		                   vec6 & ball_pred,
		                   double & time_pred,
		                   EKF & filter,
		                   game & game_state) {

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
			time_pred = const_tt::DT * (idx + 1);
			if (time_pred > time_min)
				valid_hp = true;
		}
	}

	return valid_hp;
}

void predict_ball(const double & time_pred,
                    mat & balls_pred,
                    EKF & filter) {

	//static wall_clock timer;
	//timer.tic();
	int N = (int)(time_pred/const_tt::DT);
	balls_pred = filter.predict_path(const_tt::DT,N);
	//cout << "Pred. ball time: " << 1000 * timer.toc() << " ms." << endl;
}

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

bool check_legal_ball(const vec6 & ball_est,
                        const mat & balls_predicted,
                        game & game_state) {

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
			std::cout << "DP does not have a fixed return strategy!\n";
			pos_land_des = zeros<vec>(2);
			time_land_des = 0.0;
			break;
		default:
			throw ("Algorithm is not recognized!\n");
	}
}

}
