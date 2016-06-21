/*
 * table_tennis.c
 *
 *  Created on: Jun 21, 2016
 *      Author: okoc
 */

// SL variables and kinematics

#include "table_tennis.h"

extern Matrix ballMat; // predicted ball pos and vel values for T_pred time seconds
extern Matrix racketMat; // racket strategy

/*
 * Function that calculates a racket strategy : positions, velocities and orientations
 * for each point on the predicted ball trajectory (ballMat)
 * to return it a desired point (ballLand) at a desired time (landTime)
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
