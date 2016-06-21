/*
 * table_tennis.c
 *
 *  Created on: Jun 21, 2016
 *      Author: okoc
 *
 *  Includes the common functions that we know from table_tennis_common.c in SL
 *
 */

// SL variables and kinematics

#include "table_tennis.h"

/*
 * Function that calculates a racket strategy : positions, velocities and orientations
 * for each point on the predicted ball trajectory (ballMat)
 * to return it a desired point (ballLand) at a desired time (landTime)
 *
 * TODO: shall we switch to a structure for racketMat? Should it be global?
 */
void calc_racket_strategy(Vector ballLand, double landTime) {

	int N = TPRED/TSTEP;
	racketMat = my_matrix(1, N, 1, 3*CART);
	static SL_Cstate ballIncomingPos;
	Vector ballOutVel = my_vector(1,CART);
	Vector ballInVel = my_vector(1,CART);
	Vector racketVel = my_vector(1,CART);
	Vector racketNormal = my_vector(1,CART);
	int i,j;

	for (i = 1; i <= N; i++) {
		for (j = 1; j <= CART; j++) {
			ballIncomingPos.x[j] = ballMat[i][j];
			ballInVel[j] = ballMat[i][j+CART];
		}

		// determine the desired outgoing velocity of the ball at contact
		calc_ball_vel_out(ballIncomingPos, ballLand, landTime, ballOutVel);

		//print_vec("ball out vel: ", ballOutVel);
		calc_racket_normal(ballInVel, ballOutVel, racketNormal);
		calc_racket_vel(ballInVel, ballOutVel, racketNormal, racketVel);

		//print_vec("racket vel = ",racketVel);
		//print_vec("racket normal = ",racketNormal);

		for (j = 1; j <= CART; j++) {
			racketMat[i][j] = ballIncomingPos.x[j];
			racketMat[i][j+CART] = racketVel[j];
			racketMat[i][j+2*CART] = racketNormal[j];
		}
	}

	//print_mat("Racket matrix:", racketMat);
}

/*
 * Calculate desired racket velocity given ball incoming and
 * outgoing velocities
 * Assuming a mirror law
 * Assumes no desired spin, i.e. racket velocity along the racket will be set to zero
 *
 * Output is the last parameter: racketVel
 *
 * TODO: is the CRR value accurate?
 *
 */
void calc_racket_vel(Vector velBallIn, Vector velBallOut, Vector normalRacket, Vector velRacket) {

	double velBallInAlongNormal;
	double velBallOutAlongNormal;
	double eps = CRR;

	velBallInAlongNormal = vec_mult_inner(velBallIn, normalRacket);
	velBallOutAlongNormal = vec_mult_inner(velBallOut, normalRacket);
	velBallInAlongNormal = eps * velBallInAlongNormal;

	velBallOutAlongNormal = (velBallOutAlongNormal + velBallInAlongNormal) / (1+eps);

	vec_mult_scalar(normalRacket, velBallOutAlongNormal, velRacket);
}

/*
 * Calculate desired racket normal using the mirror law
 */
void calc_racket_normal(Vector bin, Vector bout, Vector normal) {

	vec_sub(bout, bin, normal);
	// normalize
	vec_mult_scalar(normal, 1./sqrt(vec_mult_inner(normal,normal)), normal);
}

/*
 *
 *  Computes the desired outgoing velocity of the ball after contact
 *		            to hit the goal on the opponents court
 *  Input:
 *
 *		SL_CState	hitPoint:   hitting point at contact
 *		Vector		landPoint: 	landing point on the opponents court
 *
 *	Output:
 *
 *		Vector		velOut:     desired velocity of the ball
 *
 */
void calc_ball_vel_out(SL_Cstate hitPoint, Vector landPoint, double time2land, Vector velOut) {

	double x, y, z, time2Land;
	double ynet, znet, zoc;
	double alpha;
	int sign, check, i;

	static int firsttime = TRUE;
	static double zTable;

	if (firsttime) {
		firsttime = FALSE;
		zTable = floor_level - table_height + ball_radius;
	}

	velOut[_X_] = (landPoint[1] - hitPoint.x[_X_]) / time2land;
	velOut[_Y_] = (landPoint[2] - hitPoint.x[_Y_]) / time2land;
	velOut[_Z_] = (zTable - hitPoint.x[_Z_] + 0.5 * intern_gravity * sqr(time2land)) / time2land;

	//TODO: consider the air drag case
	// hack for now
	velOut[_X_] = 1.1 * velOut[_X_];
	velOut[_Y_] = 1.1 * velOut[_Y_];
	velOut[_Z_] = 1.2 * velOut[_Z_];

}

/*
 * Integrate the ball state dt seconds later. Used for prediction
 * Using symplectic Euler.
 *
 */
void integrate_ball_state(SL_Cstate ballState, SL_Cstate *ballPred, double deltat, int *bounce) {

	int i;
	double velBall;
	static double slack = 0.0001;

	// does it touch the floor?
	if (ballState.x[_Z_] >= floor_level) {

		/*****************************************/
		// Symplectic Euler for No-Contact-Situation (Flight model)

		velBall = sqrt(sqr(ballState.xd[_X_]) + sqr(ballState.xd[_Y_]) + sqr(ballState.xd[_Z_]));

		if (fabs(velBall) > slack) {

			ballPred->xdd[_X_] = -ballState.xd[_X_] * Cdrag * velBall;
			ballPred->xdd[_Y_] = -ballState.xd[_Y_] * Cdrag * velBall;
			ballPred->xdd[_Z_] = -intern_gravity - ballState.xd[_Z_] * Cdrag * velBall;
		}
		else {
			ballPred->xdd[_X_] = 0.;
			ballPred->xdd[_Y_] = 0.;
			ballPred->xdd[_Z_] = -intern_gravity;
		}

		ballPred->xd[_X_] = ballState.xd[_X_] + ballPred->xdd[_X_] * deltat;
		ballPred->xd[_Y_] = ballState.xd[_Y_] + ballPred->xdd[_Y_] * deltat;
		ballPred->xd[_Z_] = ballState.xd[_Z_] + ballPred->xdd[_Z_] * deltat;
		ballPred->x[_X_] =  ballState.x[_X_] + ballPred->xd[_X_] * deltat;
		ballPred->x[_Y_] =  ballState.x[_Y_] + ballPred->xd[_Y_] * deltat;
		ballPred->x[_Z_] =  ballState.x[_Z_] + ballPred->xd[_Z_] * deltat;

		if (check_ball_table_contact(*ballPred)) {

			//printf("Expecting a bounce!\n");
			ballPred->x[_Z_]  = floor_level - table_height  + ball_radius;
			ballPred->xd[_Z_] = -CRT * ballPred->xd[_Z_];
			ballPred->xd[_Y_] = CFTY * ballPred->xd[_Y_];
			ballPred->xd[_X_] = CFTX * ballPred->xd[_X_];
			*bounce = TRUE;
		}
	}
	else { // it touches the floor
		for (i = 1; i <= N_CART; i++) {
			ballPred->xd[i] = 0.;
			ballPred->x[i]  = ballState.x[i];
		}
	}

}

/*
 * Condition to determine if ball hits the table
 * Useful for prediction including a rebound model
 * Useful also in KF/EKF filtering.
 *
 * TODO: can we do this also for the racket?
 */
int check_ball_table_contact(SL_Cstate state) {

	int sign;
	if (dist_to_table > 0)
		sign = 1;
	else
		sign = -1;

	double center_table = dist_to_table + sign * 0.5 * table_length;
	double dist_to_center_y = fabs(center_table - state.x[_Y_]);
	double dist_to_center_x = fabs(table_center - state.x[_X_]);
	double dist_to_table_plane = state.x[_Z_] - (floor_level - table_height + ball_radius);

	if (dist_to_center_y < table_length/2.
		&& dist_to_center_x <= table_width/2.
		&& dist_to_table_plane <= 0 && state.xd[_Z_] <= 0)
		return TRUE;
	else
		return FALSE;
}

/*
 * Predict the ball for Tpred seconds
 * Using ballflight + bounce model
 *
 * TODO: do we need ballPred structure?
 */
void predict_ball_state() {

	int N = TPRED/TSTEP;
	ballMat = my_matrix(1, N, 1, 2*CART);
	int i,j;

	for (j = 1; j <= CART; j++) {
		ballMat[0][j] = ballPred.x[j];
		ballMat[0][j+CART] = ballPred.xd[j];
	}

	// predict Tpred seconds into the future
	int bounce = FALSE;

	for (i = 1; i <= N; i++) {
		integrate_ball_state(ballPred,&ballPred,TSTEP,&bounce);
		for (j = 1; j <= CART; j++) {
			ballMat[i][j] = ballPred.x[j];
			ballMat[i][j+CART] = ballPred.xd[j];
		}
	}

	//print_mat("Ball pred matrix: ", ballMat);
	/*for (i = 1; i <= N; i++) {
		for (j = 1; j <= 2*CART; j++)
			printf("%.4f  ", ballMat[i][j]);
		printf("\n");
	}*/

}

/*
 * Set desired landing position to the centre of the table and time to a
 * reasonable value
 *
 * TODO: shall we rename this? Or make it more extendible by including different strategies?
 */
void set_des_land_param(Vector ballLand, double *landTime) {

	*landTime = 0.8;

	ballLand[_X_] = 0.0;
	ballLand[_Y_] = dist_to_table - 3*table_length/4; // centre of opponents court
	ballLand[_Z_] = floor_level - table_height;// + ball_radius;
}
