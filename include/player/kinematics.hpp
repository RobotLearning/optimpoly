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

using namespace arma;
using namespace optim;

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

/**
 * @brief Calculates cartesian racket pos, vel and normal
 * given joint positions and velocities.
 *
 * C++ version of the same C-code using ARMADILLO library.
 * In the optimization we do not call this function, but stick to C version.
 * Can be used to calculate racket positions in player class.
 *
 * @param robot_joint Robot joint positions, velocities and accelerations.
 * @param robot_racket Robot racket positions, velocities and normal.
 */
void calc_racket_state(const optim::joint & robot_joint,
		               racket & robot_racket);

/**
 * @brief Rotate racket by 90 degrees to get
 * racket orientation from endeffector orientation.
 *
 * @param quat Endeffector orientation as a quaternion.
 */
void calc_racket_orient(vec4 & quat);

/**
 * @brief Return the Jacobian (linear) matrix at joint values q
 *
 * Returns the jacobian and the cartesian coordinates of endeffector
 * @param q joint values q
 * @param jac Jacobian matrix to be updated
 */
vec3 get_jacobian(const vec7 & q, mat::fixed<6,7> & jac);

#endif /* _KINEMATICS_HPP_ */
