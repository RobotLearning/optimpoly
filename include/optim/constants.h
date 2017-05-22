/**
 * @file constants.h
 *
 * @brief Constants are located here
 *
 * Constants are indices, number of degrees of freedom,
 * table tennis parameters (gravity, coefficients, ...).
 *
 *  Created on: Feb 2, 2017
 *      Author: okoc
 */

#ifndef CONSTANTS_H_
#define CONSTANTS_H_

// indices for simplifying code comprehension
#define X 0
#define Y 1
#define Z 2
#define W 3 // for quaternion
#define DX 3
#define DY 4
#define DZ 5
#define NCART 3
#define NDOF 7

// TODO: make some of these non-const so they can be loaded from player.cfg?

const double VHPY = -0.3; // VHP in front of the robot
const double DT = 0.002; // 500 Hz robot operation

/* Table Tennis Ball Variables */
const double ball_radius  = 0.02;

/* Table Tennis Racket Radius */
const double racket_radius = 0.076; // shorter axis about 15.2 cm, longer approx 15.5 - 15.6

/* Minimum number of observations to start filter */
const int min_obs = 5;

#endif /* INCLUDE_CONSTANTS_H_ */
