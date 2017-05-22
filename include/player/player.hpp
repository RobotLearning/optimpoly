/**
 * @file player.hpp
 *
 * @brief Main class for playing Table Tennis declared here.
 *
 *  Created on: Feb 9, 2017
 *      Author: okoc
 */

#ifndef PLAYER_HPP_
#define PLAYER_HPP_

#include "kalman.h"
#include "optim.h"

using namespace arma;

enum algo {
	VHP,
	FIXED,
	LAZY,
};

enum game { //trial state
	AWAITING,
	ILLEGAL,
	LEGAL,
	HIT,
};

enum mode_operate { // mode of operations
	TEST_SIM, // optim wont be detached
	SL_SIM, // sim but optim will be detached
	REAL_ROBOT, // outlier detection
};

/**
 * @brief Desired/actual joint positions, velocities, accelerations.
 *
 * output of main Player function play()
 */
typedef struct {
	vec7 q;
	vec7 qd;
	vec7 qdd;
} joint;

/**
 *
 * @brief Table Tennis Player class for playing Table Tennis.
 *
 * The methods play() or cheat() must be called every DT milliseconds.
 */
class Player {

private:

	// data fields
	EKF & filter; // filter for the ball estimation
	vec2 ball_land_des; // desired landing position
	double time_land_des; // desired landing time
	double time2return; // fixed return time for robot
	racketdes racket_params;
	optim optim_params;
	coptim coparams;
	vec7 q_rest_des; // desired resting joint state
	double t_cum; // counting time stamps for resetting filter
	mat observations; // for initializing filter
	mat times; // for initializing filter

	// flags and related fields
	algo alg; // algorithm (fixed player, vhp, etc.)
	game game_state; // ball awaiting, detected bouncing legally/illegally, or was hit
	bool mpc; // apply corrections
	bool valid_obs; // ball observed is valid (new ball and not an outlier)
	int verbose; // level of verbosity (printing, OFF = 0, LOW = 1, HIGH = 2)
	int num_obs; // number of observations received
	mode_operate mode; // sim vs. real robot

	// ball estimation
	void estimate_ball_state(const vec3 & obs);

	// optimization for different players
	void optim_fixedp_param(const joint & qact); // run optimizer for FIXED player
	void optim_lazy_param(const joint & qact);
	void optim_vhp_param(const joint & qact); // run VHP player

	bool check_update(const joint & qact) const; // flag for (re)running optimization
	bool predict_hitting_point(vec6 & ball_pred, double & time_pred);
	void predict_ball(mat & balls_pred) const;
	void calc_next_state(const joint & qact, joint & qdes);

public:

	Player(const vec7 & q0, EKF & filter,
			algo alg = FIXED, bool mpc = false,
			int verbose = 0, mode_operate mode = TEST_SIM);
	~Player();

	// auxiliary function, public interface for filter test performance
	vec6 filt_ball_state(const vec3 & obs);
	void reset_filter(double std_model, double std_noise);

	// main function
	void play(const joint & qact, const vec3 & ball_obs, joint & qdes);

	// cheat function for simulation (with exact knowledge of ball state)
	void cheat(const joint & qact, const vec6 & ballstate, joint & qdes);

};

// ball estimation and filter constructor/state initialization
EKF init_filter(double std_model = 0.001, double std_noise = 0.001, bool spin = false);
void estimate_prior(const mat & observations,
		            const vec & times,
					EKF & filter,
					bool mode);
bool check_new_obs(const vec3 & obs, double tol);
bool check_reset_filter(const bool newball, const int verbose);

// movement generation
void generate_strike(const optim & params, const joint & qact,
		             const vec7 & q_rest_des, const double time2return,
		            mat & Q, mat & Qd, mat & Qdd);
void gen_3rd_poly(const rowvec & times, const vec7 & a3, const vec7 & a2, const vec7 & a1, const vec7 & a0,
		     mat & Q, mat & Qd, mat & Qdd);
void set_bounds(double *lb, double *ub, double SLACK, double Tmax);

// racket calculations
racketdes calc_racket_strategy(const mat & balls_predicted,
		                       const vec2 & ball_land_des, const double time_land_des,
							   racketdes & racket_params);
bool check_legal_ball(const vec6 & ball_est, const mat & balls_predicted, game & game_state);
void check_legal_bounce(const vec6 & ball_est, game & game_state);

#endif /* PLAYER_HPP_ */
