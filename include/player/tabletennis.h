/**
 * @file tabletennis.h
 *
 * @brief Header for table tennis ball prediction functions.
 *
 * Mostly taken from table_tennis_common.h file.
 *
 *  Created on: Feb 1, 2017
 *      Author: okoc
 */

#ifndef TABLETENNIS_H_
#define TABLETENNIS_H_

#include "constants.h"
#include "table.h"

using namespace arma;

/**
 * @brief Racket positions, velocities and normal
 * used to predict racket contact.
 */
struct racket {
	vec3 pos;
	vec3 vel;
	vec3 normal;
};

/**
 * @brief Information about the status of the game
 *
 * This information is useful in SIMULATION. For instance
 * when we're running tests for performance of optimizers, these
 * are useful to collect statistics.
 */
struct status {
	bool hit = false; //!< robot hit the incoming ball
	bool has_bounced = false; //!< the ball has bounced on robot court (without checking for legal bounce)
	bool legal_bounce = false; //!< the ball has bounced only once on the robot court
	bool has_landed = false; //!< the ball has landed on the other side (without checking for legal land)
	bool legal_land = false; //!< ball has landed legally on the other side (bounced once only on robot court)
	bool touched_ground = false; //!< the ball touched the ground level (vel. zeroed)
};

/**
 * @brief Ball parameters used to predict future ball path
 * and to calculate desired racket parameters.
 */
struct ball_params {

	/* Contact Coefficients */
	double CRT = 0.88; //!< coefficient of restitution for the table (i.e. rebound z-velocity multiplier)
	double CFTY = 0.72; //!< coefficient of table contact model on Y-direction
	double CFTX = 0.68; //!< coefficient of table contact model on X-direction
	double CRR = 0.78; //!< coefficent of restitution for racket

	double Cdrag = 0.1414; 	//!< Air drag coefficient
	double gravity = -9.802; //!< gravity
	double Clift = 0.001; //!< coefficient of lift for the magnus force
	double mu = 0.10; //!< dynamic coefficient of friction
	double init_topspin = -50.0; //!< initial topspin amount
};

// flags for table tennis players
static const bool CHECK_CONTACTS = true; // turn off for simplified debugging

/**
 * @brief Table Tennis ball prediction methods
 *
 * Spin and spinless ball flight models can be used the predict the
 * next ball state. Checking for contact, moreover, enables to predict the
 * future of the table tennis trial, given racket and the table.
 */
class TableTennis {

private:

	bool SPIN_MODE; // turn on prediction with a spin model
	bool VERBOSE;
	status stats;

	ball_params params; // ball prediction parameters
	vec3 ball_pos;
	vec3 ball_vel;
	vec3 ball_spin; // ball angular velocity = 0 if spin mode is turned OFF

	// init topspin function useful in spin initialization
	void init_topspin(const double val = -50);

	vec3 flight_model() const;
	vec3 drag_flight_model() const;
	vec3 table_contact_model(const vec3 & ball_vel_in) const;
	void symplectic_euler(const double dt, const vec3 & ball_acc,
			              vec3 & ball_cand_pos, vec3 & ball_cand_vel) const;

	// contact models
	void check_contact(const racket & robot_racket,
			           vec3 & ball_cand_pos, vec3 & ball_cand_vel); // calls the contact functions below

	void check_legal_bounce(const vec3 & ball_cand_pos, const vec3 & ball_cand_vel);
	void check_legal_land(const vec3 & ball_cand_pos, const vec3 & ball_cand_vel);
	// Check contact to table
	void check_ball_table_contact(const vec3 & ball_cand_pos, vec3 & ball_cand_vel);
	// Check contact with net
	void check_ball_net_contact(vec3 & ball_cand_pos, vec3 & ball_cand_vel) const;
	// Check contact with racket
	void check_ball_racket_contact(const racket & robot,
			       const vec3 & ball_cand_pos, vec3 & ball_cand_vel);
	// Check if it hits the ground...
	void check_ball_ground_contact(vec3 & ball_cand_pos, vec3 & ball_cand_vel);

public:

	// initialization
	TableTennis(bool spin = false, bool verbose = false);
	TableTennis(const vec6 & ball_state, bool spin = false, bool verbose = false);

	bool has_legally_landed() const;
	bool has_legally_bounced() const;
	void reset_stats();
	void load_params(const std::string & file_name_relative);

	vec3 get_ball_position() const;
	vec3 get_ball_velocity() const;
	vec6 get_ball_state() const;
	void set_ball_state(const vec6 & ball_state);
	void set_topspin(const double val);
	// set reasonable ball positions and velocities for table tennis
	void set_ball_gun(double std, int ballgun_side = 1);

	// ball prediction functions
	void integrate_ball_state(const double dt);
	void integrate_ball_state(const racket & robot, const double dt);

	// calculating desired racket
	void calc_des_racket_normal(const mat & v_in, const mat & v_out, mat & normal) const;
	void calc_des_ball_out_vel(const vec2 & ball_land_des,
							   const double time_land_des,
							   const mat & balls_predicted, mat & balls_out_vel) const;
	void calc_des_racket_vel(const mat & vel_ball_in, const mat & vel_ball_out,
			                 const mat & racket_normal, mat & racket_vel) const;
};

// Functions set as pointers for Extended Kalman Filtering
vec calc_next_ball(const vec & xnow, const double dt, const void *fp);
vec calc_next_ball(const racket & robot, const vec & xnow, const double dt);
vec calc_spin_ball(const vec & xnow, const double dt, const void *fp);

// Used for lookup table
void predict_till_net(vec6 & ball_est);

#endif /* end of TABLETENNIS_H_ */
