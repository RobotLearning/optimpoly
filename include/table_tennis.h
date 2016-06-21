/*
 * table_tennis.h
 *
 *  Created on: Jun 21, 2016
 *      Author: okoc
 */

#ifndef TABLE_TENNIS_H_
#define TABLE_TENNIS_H_

#include "SL.h"
#include "SL_user.h"

// defines
#define CART 3
#define DOF 7
#define OPTIM_DIM 2*DOF+1
#define CONSTR_DIM 3*CART
#define MAX_VEL 200
#define MAX_ACC 200
#define TSTEP 0.01
#define TPRED 1.0

/* racket computations */
void calc_racket_strategy(double *ballLand, double landTime);
void calc_ball_vel_out(SL_Cstate hitPoint, Vector landPoint, double time2land, Vector velOut);
void calc_racket_normal(Vector bin, Vector bout, Vector normal);
void calc_racket_vel(Vector velBallIn, Vector velBallOut, Vector normalRacket, Vector velRacket);

#endif /* TABLE_TENNIS_H_ */
