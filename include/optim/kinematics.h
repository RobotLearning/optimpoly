/**
 * @file kinematics.h
 *
 * @brief Here we include the kinematics related functions taken from SL
 *
 *  Created on: Jun 22, 2016
 *      Author: okoc
 */

#ifndef KINEMATICS_H_
#define KINEMATICS_H_

#include "math.h"

#define NDOF 7
#define NCART 3
#define NQUAT 4
#define NLINK 6

/*! defines that are used to parse the config and prefs files */
#define MIN_THETA    1
#define MAX_THETA    2
#define THETA_OFFSET 3

// dimensions of the robot
#define ZSFE 0.346              //!< z height of SAA axis above ground
#define ZHR  0.505              //!< length of upper arm until 4.5cm before elbow link
#define YEB  0.045              //!< elbow y offset
#define ZEB  0.045              //!< elbow z offset
#define YWR -0.045              //!< elbow y offset (back to forewarm)
#define ZWR  0.045              //!< elbow z offset (back to forearm)
#define ZWFE 0.255              //!< forearm length (minus 4.5cm)

// FROM MDEFS.H FILE
#define Power(x, y)	(pow((double)(x), (double)(y)))
#define Sqrt(x)		(sqrt((double)(x)))

#define Abs(x)		(fabs((double)(x)))

#define Exp(x)		(exp((double)(x)))
#define Log(x)		(log((double)(x)))

#define Sin(x)		(sin((double)(x)))
#define Cos(x)		(cos((double)(x)))
#define Tan(x)		(tan((double)(x)))

#define ArcSin(x)       (asin((double)(x)))
#define ArcCos(x)       (acos((double)(x)))
#define ArcTan(x)       (atan((double)(x)))

#define Sinh(x)          (sinh((double)(x)))
#define Cosh(x)          (cosh((double)(x)))
#define Tanh(x)          (tanh((double)(x)))


#ifndef E
#define E		2.71828182845904523536029
#endif
#ifndef Pi
#define Pi		3.14159265358979323846264
#endif
#define Degree		0.01745329251994329576924

/**
 * @brief Calculates cartesian racket pos, vel and normal
 * given joint positions and velocities.
 *
 * Calls the kinematics and the (geometric) jacobian functions
 * to compute the racket quantities.
 *
 * @param q Joint positions.
 * @param qdot Joint velocities.
 * @param pos Racket positions.
 * @param vel Racket velocities.
 * @param normal Racket normals.
 */
void calc_racket_state(const double q[NDOF],
		               const double qdot[NDOF],
					   double pos[NCART],
					   double vel[NCART],
					   double normal[NCART]);

/**
 * @brief Return racket positions, normal and the jacobian
 *
 * Makes the input pos vector equal to racket positions for given joints q
 * Makes the input matrix equal to the jacobian at input q
 *
 * Useful to test derivatives of kinematics
 */
void calc_racket_state(const double q[NDOF],
                        double pos[NCART],
                        double normal[NCART],
                        double jacobi[NCART][NDOF]);

/** @brief Returns the cartesian racket positions */
void get_position(double q[NDOF]);

#endif /* KINEMATICS_H_ */
