/**
 * @file constants.h
 *
 * @brief Constants are located here
 *
 * Constants are indices, number of degrees of freedom,
 * table tennis parameters (gravity, coefficients, radii, ...).
 *
 *  Created on: Feb 2, 2017
 *      Author: okoc
 */

#ifndef CONSTANTS_H_
#define CONSTANTS_H_

namespace const_tt {

// indices for simplifying code comprehension
const int X = 0;
const int Y = 1;
const int Z = 2;
const int W = 3; //!< for quaternion, not used
const int DX = 3;
const int DY = 4;
const int DZ = 5;
const int NBLOBS = 2; //!< In the table tennis setup, we have 4 cameras = 2 blobs.
const int NCART = 3;
const int NDOF = 7;

const double DT = 0.002; //!< 500 Hz robot operation

/* Table Tennis Ball Variables */
const double ball_radius  = 0.02; //!< table tennis standard ball radius
const double racket_radius = 0.076; //!< shorter axis about 15.2 cm, longer approx 15.5 - 15.6

}

#endif /* INCLUDE_CONSTANTS_H_ */
