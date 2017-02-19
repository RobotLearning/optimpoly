/*
 * table_tennis.h
 *
 *  Created on: Jun 21, 2016
 *      Author: okoc
 */

#ifndef TABLE_TENNIS_H_
#define TABLE_TENNIS_H_

#include "math.h"
#include "string.h"
#include "SL.h"
#include "constants.h"
#include "table.h"

// defines
#define CART 3
#define TSTEP 0.01
#define TPRED 1.0
#define DOF 7

extern Matrix ballMat; // predicted ball pos and vel values for T_pred time seconds
extern Matrix racketMat; // racket strategy: desired positions, velocities and normal saved in matrix

/* racket computations */
void calc_racket_strategy(double *ballLand, double landTime);
void calc_ball_vel_out(SL_Cstate hitPoint, Vector landPoint, double time2land, Vector velOut);
void calc_racket_normal(Vector bin, Vector bout, Vector normal);
void calc_racket_vel(Vector velBallIn, Vector velBallOut, Vector normalRacket, Vector velRacket);

/* ball related methods */
void set_des_land_param(double *ballLand, double *landTime);
void predict_ball_state(double *b0, double *v0);
void integrate_ball_state(SL_Cstate ballState, SL_Cstate *ballPred, double deltat, int *bounce); //ball prediction
int check_ball_table_contact(SL_Cstate state);

/* first order hold to get racket parameters at time T */
void first_order_hold(double *ballPred, double *racketVel, double *racketNormal, double T);

/* lookup table, i.e. kNN with k = 1 */
int lookup(const Matrix lookupTable, const double* b0, const double* v0, double *x);

#endif /* TABLE_TENNIS_H_ */
