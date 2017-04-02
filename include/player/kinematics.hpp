/*
 * kinematics.hpp
 *
 *  Created on: Feb 12, 2017
 *      Author: okoc
 */

#ifndef _KINEMATICS_HPP_
#define _KINEMATICS_HPP_

#include "player.hpp"
#include "tabletennis.h"

/*! defines for the preference and config files */
#define CONFIG   "config/"
#define PREFS    "prefs/"

// robot constants
static const double ZSFE = 0.346; // z height of SAA axis above ground
static const double ZHR = 0.505;  // length of upper arm until 4.5cm before elbow link
static const double YEB = 0.045;  // elbow y offset
static const double ZEB = 0.045;  // elbow z offset
static const double YWR = -0.045; // elbow y offset (back to forewarm)
static const double ZWR = 0.045;  // elbow z offset (back to forearm)
static const double ZWFE = 0.255; // forearm length (minus 4.5cm)

static const int PALM = 5; // for kinematics

using namespace arma;

typedef struct {
	vec3 x; // pos
	vec3 o; // orientation
} eff; // endeffector

typedef struct {
	vec3 x; // positions
	vec4 q; // orientations
} pose; // used for base pose of the robot


void calc_racket_state(const joint & robot_joint,
		               racket & robot_racket);
void calc_racket_orient(vec4 & quat);
void mult_two_quats(const vec4 & q1, const vec4 & q2, vec4 & q3);
void rotate_to_quat(const mat33 & R, vec4 & quat);
void revolute_jac_col(const vec3 & p, const vec3 & pi, const vec3 & zi, vec6 & col);
void kinematics(const vec7 & q, mat & Xlink, mat & Xorigin, mat & Xaxis, cube & Amats);
void jacobian(const mat & lp, const mat & jop, const mat & jap, mat & jac);

#endif /* _KINEMATICS_HPP_ */
