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

// indices for simplifying code comprehension
#define X 0
#define Y 1
#define Z 2
#define W 3 //!< for quaternion, not used
#define DX 3
#define DY 4
#define DZ 5
#define NBLOBS 2 //!< In the table tennis setup, we have 4 cameras = 2 blobs.
#define NCART 3
#define NDOF 7

const double DT = 0.002; //!< 500 Hz robot operation

/* Table Tennis Ball Variables */
const double ball_radius  = 0.02; //!< table tennis standard ball radius
const double racket_radius = 0.076; //!< shorter axis about 15.2 cm, longer approx 15.5 - 15.6



#endif /* INCLUDE_CONSTANTS_H_ */
