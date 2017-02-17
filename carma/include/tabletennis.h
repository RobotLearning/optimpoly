/*
 * tabletennis.h
 *
 * Header for table tennis ball prediction functions.
 * Mostly taken from table_tennis_common.h file residing in SL/barrett.
 *
 *  Created on: Feb 1, 2017
 *      Author: okoc
 */

#ifndef TABLETENNIS_H_
#define TABLETENNIS_H_

#include "table.h"
#include "constants.h"

using namespace arma;

// flags for table tennis players
static const bool SPIN_MODE = false; // turn on prediction with a spin model
static const bool CHECK_CONTACTS = true; // turn off for simplified debugging
static const bool VERBOSE = true;

// functions outside of Table Tennis class
mat33 quat2mat(const vec4 & q);
void table_contact_model(vec3 & ball_spin, vec3 & ball_vel);
void racket_contact_model(const vec3 & racket_vel, const vec3 & racket_normal, vec3 ball_vel);
vec calc_next_ball(const vec & xnow, double dt);

class TableTennis { // Table Tennis ball prediction methods

private:

	vec3 ball_pos;
	vec3 ball_vel;
	vec3 ball_spin; // ball angular velocity = 0 if spin mode is turned OFF

	vec3 racket_pos;
	vec3 racket_vel;
	vec4 racket_orient;

	// init topspin function useful in spin initialization
	void init_topspin();

	vec3 flight_model() const;
	vec3 drag_flight_model() const;
	void symplectic_euler(const vec3 & ball_acc,
			              vec3 & ball_cand_pos, vec3 & ball_cand_vel, double dt) const;

	// contact models
	void check_contact(vec3 & ball_cand_pos, vec3 & ball_cand_vel); // calls the contact functions below
	// Check contact to table
	void check_ball_table_contact(const vec3 & ball_cand_pos, vec3 & ball_cand_vel);
	// Check contact with net
	void check_ball_net_contact(vec3 & ball_cand_pos, vec3 & ball_cand_vel) const;
	// Check contact with racket
	void check_ball_racket_contact(const vec3 & ball_cand_pos, vec3 & ball_cand_vel) const;
	// Check if it hits the ground...
	void check_ball_ground_contact(vec3 & ball_cand_pos, vec3 & ball_cand_vel) const;

public:

	// initialization
	TableTennis();
	TableTennis(vec6 ball_state);

	vec3 get_ball_position() const;
	vec3 get_ball_velocity() const;

	// set reasonable ball positions and velocities for table tennis
	void set_ball_state(double std);

	// ball prediction functions
	void integrate_ball_state(double dt);

	// used as an interface for the ball prediction
	friend vec calc_next_ball(const vec & xnow, double dt);
};

#endif /* end of TABLETENNIS_H_ */
