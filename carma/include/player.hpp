/*
 * player.hpp
 *
 *  Created on: Feb 9, 2017
 *      Author: okoc
 */

#ifndef PLAYER_HPP_
#define PLAYER_HPP_

#include "kalman.h"
#include "optimpoly.h"

using namespace arma;

typedef struct {
	vec7 q;
	vec7 qd;
	vec7 qdd;
} joint; // output of main Player function play()

class Player { // Table Tennis Player

private:

	EKF filter; // filter for the ball estimation
	vec2 ball_land_des; // desired landing position
	double time_land_des; // desired landing time
	double time2return; // fixed return time for robot
	racketdes racket_params;
	optim optim_params;
	coptim coparams;
	vec7 q_rest_des; // desired resting joint state
	bool moving; // robot is moving or not

	void estimate_ball_state(const vec3 & obs);
	void estimate_prior(const mat & observations, const vec & times);
	void calc_optim_param(const joint & qact); // run optimizer
	void predict_ball(mat & balls_pred);
	racketdes calc_racket_strategy(const mat & balls_predicted);
	void calc_next_state(const joint & qact, joint & qdes);
	void generate_strike(const joint & qact, mat & Q, mat & Qd, mat & Qdd) const;

public:

	Player();
	Player(const vec7 & q0, const EKF & filter);

	// auxiliary function, public interface for filter test performance
	vec6 filt_ball_state(const vec3 & obs);

	// friend function to test racket strategy
	friend racketdes send_racket_strategy(Player & robot);

	// main function
	void play(const joint & qact, const vec3 & ball_obs, joint & qdes);
};

EKF init_filter(); // init filter for ball estimation
void gen_3rd_poly(const rowvec & times, const vec7 & a3, const vec7 & a2, const vec7 & a1, const vec7 & a0,
		     mat & Q, mat & Qd, mat & Qdd);
void calc_des_ball_out_vel(const vec2 & ball_land_des,
						   const double time_land_des,
						   const mat & balls_predicted, mat & balls_out_vel);
void calc_des_racket_vel(const mat & vel_ball_in, const mat & vel_ball_out,
		                 const mat & racket_normal, mat & racket_vel);
void calc_des_racket_normal(const mat & v_in, const mat & v_out, mat & normal);
bool check_legal_ball(const mat & balls_predicted);
bool check_new_obs(const vec3 & obs);
racketdes send_racket_strategy(const vec7 & qinit, const vec6 & ball_state,
		                    const double T);
void set_bounds(double *lb, double *ub, double SLACK, double Tmax);
double** my_matrix(int nrl, int nrh, int ncl, int nch);


#endif /* PLAYER_HPP_ */
