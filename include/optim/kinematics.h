/**
 * @file kinematics.h
 *
 * @brief Here we include the kinematics related functions taken from SL
 *
 *  Created on: Jun 22, 2016
 *      Author: okoc
 */

#pragma once
#include "constants.h"
#include "math.h"

using namespace const_tt;

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
void calc_racket_state(const double q[NDOF], const double qdot[NDOF],
                       double pos[NCART], double vel[NCART],
                       double normal[NCART]);

/**
 * @brief Return racket positions, normal and the jacobian
 *
 * Makes the input pos vector equal to racket positions for given joints q
 * Makes the input matrix equal to the jacobian at input q
 *
 * Useful to test derivatives of kinematics
 */
void calc_racket_state(const double q[NDOF], double pos[NCART],
                       double normal[NCART], double jacobi[2 * NCART][NDOF]);

/** @brief Returns the cartesian racket positions */
void get_position(double q[NDOF]);
