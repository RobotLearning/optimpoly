/**
 * @file kinematics.hpp
 *
 * @brief Kinematics C++ functions to calculate racket state
 *
 * Useful for calling from Player class.
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

//! robot constants
#ifndef KINEMATICS_H_
static const double ZSFE = 0.346; //! z height of SAA axis above ground
static const double ZHR = 0.505;  //! length of upper arm until 4.5cm before elbow link
static const double YEB = 0.045;  //! elbow y offset
static const double ZEB = 0.045;  //! elbow z offset
static const double YWR = -0.045; //! elbow y offset (back to forewarm)
static const double ZWR = 0.045;  //! elbow z offset (back to forearm)
static const double ZWFE = 0.255; //! forearm length (minus 4.5cm)
#endif

static const int PALM = 5; //! for kinematics

using namespace arma;

/**
 * @brief Endeffector positions and normal.
 */
struct eff {
	vec3 x; //! pos
	vec3 o; //! orientation
};
/**
 * @brief Pose used for base pose of the robot
 */
struct pose {
	vec3 x; //! positions
	vec4 q; //! orientations
};

void calc_racket_state(const joint & robot_joint,
		               racket & robot_racket);
void calc_racket_orient(vec4 & quat);
vec3 get_jacobian(const vec7 & q, mat::fixed<6,7> & jac);

#endif /* _KINEMATICS_HPP_ */
