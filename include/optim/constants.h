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

static const double VHPY = -0.3; // VHP in front of the robot
static const double DT = 0.002; // 500 Hz robot operation

/************************************* constants  ***********************/

/* Contact Coefficients */
/* coefficient of restitution for the table (i.e. rebound z-velocity multiplier) */
static const double CRT = 0.88;
/* coefficient of table contact model on Y-direction */
static const double CFTY = 0.72;
/* coefficient of table contact model on X-direction */
static const double CFTX = 0.68;
/* coefficent of restitution for racket */
static const double CRR = 0.78;

/* Air drag coefficient */
static const double Cdrag = 0.1414;

/* for simulating different gravities */
static const double gravity = -9.802;
/* coefficient of lift for the magnus force */
static const double Clift = 0.001;

/* Table Tennis Ball Variables */
static const double ball_radius  = 0.02;

/* Table Tennis Racket Radius */
static const double racket_radius = 0.076; // shorter axis about 15.2 cm, longer approx 15.5 - 15.6



#endif /* INCLUDE_CONSTANTS_H_ */
