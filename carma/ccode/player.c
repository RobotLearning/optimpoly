/*
 *
 * Table Tennis Player C code
 *
 * player.c
 *
 *  Created on: Feb 6, 2017
 *      Author: okoc
 */

#include "SL0.h"
#include "optims.h"
#include "player.h"
#include "constants.h"
#include "table.h"
#include "stdio.h"
#include "math.h"
#include "stdlib.h"
#include "utils.h"

/*
 * Static global variable
 */
static const int MPC_MODE = FALSE; // turn on mpc

// flags for table tennis players
static const int PRED_SPIN_MODE = FALSE; // turn on prediction with a spin model

/*
 * Table Tennis Player
 * plays every 2ms with this function
 *
 * Estimation, planning and execution are done here.
 *
 * Updates the global joint_des_state structure
 *
 */
void play() {

	int i,j;
	static int firsttime = TRUE;
	// time for moving from desired hitting state to desired resting state
	static double returnToRestTime = 1.0;
	// reset kalman filter with this variable
	static int initKF;
	// we have this new variable in optim task
	static double landTime;
	static Vector ballLand;
	/* thread optimization variables */
	static nlopt_thread_data *optim_data;
	/* bounce variable is used in legal ball detection */
	static int bounce = FALSE;
	/* number of observations so far needed for MPC to be more cautious*/
	static int numObs = 0;
	/* whether the optimization was launched or not */
	static int launch_optim = FALSE;

	if (firsttime) {
		firsttime = FALSE;
		ballLand = my_vector(1,N_CART);
		set_robot_land_param(ballLand, &landTime);
		initKF = TRUE;
		optim_data = (nlopt_thread_data *)malloc(sizeof(nlopt_thread_data));
		for (i = 1; i <= N_DOFS; i++) {
			optim_data->target[i].th = init_joint_state[i].th;
			optim_data->target[i].thd = 0.0;
		}
		optim_data->hitTime = 1.0;
	}

	// estimate ball state with a filter
	calc_ball_state(&initKF, &bounce, &numObs);

	// initialize optimization and get the hitting parameters
	calc_param(ballLand, landTime, optim_data, &launch_optim, bounce, numObs);

	// generate movement or calculate next desired step
	calc_next_state(*optim_data, &launch_optim, returnToRestTime);

}

/*
 * Estimate the balls position and velocity with a filter
 *
 * If the SIM_FILTER_MODE is OFF then we directly copy the sim_ball_state
 * in simulation mode. This is useful for debugging filters.
 *
 */
void calc_ball_state(int *init_filter, int *bounce, int *num_obs) {

	int i;
	//running_window_est(&ballPred, &timeOfFlight);
	//*bounce = rkf(&ballPred, racket_state, racket_orient);
	//ekf_carma(&ballPred, blobs[BLOB1].blob, racket_state, racket_orient, init_filter);
	//*bounce = kf(&ballPred, blobs[BLOB1].blob, racket_state, racket_orient, init_filter);
	*num_obs += 1;
	for(i = 1; i <= N_CART; i++) {
		ballPred.x[i] = sim_ball_state.x[i];
		ballPred.xd[i] = sim_ball_state.xd[i];
	}

}

/*
 * If MPC mode is turned on
 * then we use MPC to correct for new ball observation
 *
 * Returns optimization launched flag
 */
void calc_param(Vector ballLand, double landTime,
		nlopt_thread_data *optim_data, int *moving, int bounce, int numobs) {

	if (MPC_MODE) {
		calc_mpc_param(ballLand, landTime, optim_data, bounce, numobs);
	}
	else {
		// initialize optimization and get the hitting parameters
		calc_optim_param(ballLand, landTime, moving, optim_data);
	}
}


/*
 * Calculate the optimization parameters using an NLOPT nonlinear optimization algorithm
 * in another thread
 *
 * The optimized parameters are: qf, qf_dot, T
 * assuming T_return and q0 are fixed
 *
 * Returns optimization launched flag
 *
 */
static void calc_optim_param(Vector ballLand, double landTime, int *moving,
		                     nlopt_thread_data *optim_data) {

	int i;
	static pthread_t nlopt_thread;
	static int rc;

	if (!(*moving) && ballPred.x[_Y_] > (dist_to_table - table_length/2)
			&& ballPred.xd[_Y_] > 1.0) {

		if (predict_ball_path(&ballPred, FALSE)) { // ball is legal
			*moving = TRUE;
			//print_filter_accuracy();
			//double initTime = getTime();
			calc_racket_strategy(ballLand, landTime);
			//printf("Racket strategy took %f ms\n", (getTime() - initTime)/1e3);
			// run optimization in another thread
			if ((rc = pthread_create(&nlopt_thread, NULL, mthread_nlopt_optim, optim_data))) {
				fprintf(stderr, "Error: pthread_create, rc: %d\n", rc);
				//freeze();
			}
		}
		else {
			//printf("Ball is illegal! Not moving!\n");
		}
	}
}

/*
 * Calculate the optimization parameters using an NLOPT nonlinear optimization algorithm
 * in another thread
 *
 * The optimized parameters are: qf, qf_dot, T
 * assuming T_return and q0 are fixed. First using a lookup to initialize a reasonable estimate.
 *
 * We activate MPC only if
 * 1. new observations are available
 * 2. only with a certain frequency [we can increase to 40 on sim]
 * 3. only if optimization thread is not active
 * 4. only if new observations (e.g. filtered state) indicate the ball is infront of the robot
 * 5. only if the ball is still legal
 * 6. only if number of ball observations are higher than a limit
 *
 *
 */
static void calc_mpc_param(Vector ballLand, double landTime, nlopt_thread_data *optim_data,
		                   int bounce, int numobs) {

	int i;
	static pthread_t nlopt_thread;
	static int rc;
	static int firsttime = TRUE;
	static double lastTime = 0.0;
	static int FREQ_MPC = 40;
	static int MIN_OBS = 20;
	static double MPC_LIMIT_Y;
	double currentTime = get_time();
	int activateUpdate = (currentTime - lastTime) > (1e6) * (1.0/FREQ_MPC);
	int passedRobot = ballPred.x[_Y_] > cart_state[1].x[_Y_];

	if (firsttime) {
		MPC_LIMIT_Y = dist_to_table/2;
		firsttime = FALSE;
	}

	if (check_new_ball_obs(ballPred) && activateUpdate &&
			(!passedRobot) && (ballPred.xd[_Y_] > 1.0) &&
			(ballPred.x[_Y_] < MPC_LIMIT_Y) && (numobs >= MIN_OBS)) {

		if (predict_ball_path_mpc(&ballPred,bounce)) {
			calc_racket_strategy(ballLand,landTime);

			lastTime = currentTime;
			// for good initialization, and to immediately start hitting
			//copy_lookup_param(optim_data,target_joint_state,hitTime);

			for (i = 1; i <= N_DOFS; i++) {
				current_joint_state[i].th = joint_state[i].th;
				current_joint_state[i].thd = joint_state[i].thd;
			}
			// run optimization in another thread
			if (!busy) {
				printf("Predicted when ball_y : %f\n", ballPred.x[_Y_]);
				if ((rc = pthread_create(&nlopt_thread, NULL, mthread_nlopt_optim, optim_data))) {
					fprintf(stderr, "Error: pthread_create, rc: %d\n", rc);
					//freeze();
				}
			}
		}
	}
}

/*
 * Unfold the next desired state of the 3rd order polynomials in joint space
 * If movement finishes then the desired state velocities and accelerations are zeroed.
 *
 * Multithreading : if after initial lookup, the computation in the other thread terminates, then we
 * synchonize the values between the threads and compute the next desired states given these new polynomial
 * parameters (qf, qf_dot and T_hit)
 *
 */
void calc_next_state(nlopt_thread_data optim_data, int* moving, double returnToRestTime) {

	int i;
	static int firsttime = TRUE;
	// the time it takes to reach from initial state (start) to hitting point
	// = predicted time it takes the ball from current (est.) position to the hitting point
	static double hitTime = 1.0;
	// desired arm configuration, can be used for hitting state
	static SL_DJstate target_joint_state[N_DOFS+1];

	if (firsttime) {
		firsttime = FALSE;
		bzero((char *)&(target_joint_state[1]),
				N_DOFS * sizeof(target_joint_state[1]));
		for (i = 1; i <= N_DOFS; i++) {
			target_joint_state[i].th = joint_state[i].th;
			target_joint_state[i].thd = joint_state[i].thd;
		}
	}

	// this should be only for MPC?
	if (mthread_sync_optim(optim_data, target_joint_state, &hitTime)) {
		*moving = TRUE;
		trj_time = DT;
	}

	// make sure we update after optim finished
	if (*moving) {
		if (trj_time <= DT + DT/2) {
			for (i = 1; i <= N_DOFS; i++) {
				joint_des_state[i].th = joint_state[i].th;
				joint_des_state[i].thd = joint_state[i].thd;
			}
		}
		if (trj_time <= hitTime + returnToRestTime + DT/2) {
			// update the desired state
			calcNextStateForPoly3rd(joint_des_state, hitTime, target_joint_state,
					hitTime + returnToRestTime, init_joint_state, trj_time);
			trj_time = trj_time + DT;
		}
		else {
			trj_time = 0.0; // hitting process will finish
			printf("Hitting process finished\n");
			*moving = FALSE;
			for (i = 1; i <= N_DOFS; i++) {
				joint_des_state[i].thd = 0.0;
				joint_des_state[i].thdd = 0.0;
			}
		}
	}
}

/*
 * Function that calculates a racket strategy : positions, velocities and racket normal
 * for each point on the predicted ball trajectory (ballMat)
 * to return it a desired point (ballLand) at a desired time (landTime)
 *
 */
void calc_racket_strategy(Vector ballLand, double landTime) {

	int N = TPRED/TSTEP;
	//racketDes = malloc((N+1)* sizeof(Racket));
	racketDes = (Racket *)calloc(N+1, sizeof(Racket));
	Vector ballOutVel = my_vector(1,N_CART);
	Vector ballInVel = my_vector(1,N_CART);
	Vector racketVel = my_vector(1,N_CART);
	Vector racketNormal = my_vector(1,N_CART);
	int i,j;

	for (i = 1; i <= N; i++) {
		for (j = 1; j <= N_CART; j++) {
			ballInVel[j] = ballPath[i].xd[j];
		}
		// determine the desired outgoing velocity of the ball at contact
		calc_des_ball_out_vel(ballPath[i], ballLand, landTime, ballOutVel);
		// calculate the desired racket normal at the ball position
		calc_des_racket_normal(ballInVel, ballOutVel, racketNormal);
		// calculate the desired racket velocity at the ball position
		calc_des_racket_vel(ballInVel, ballOutVel, racketNormal, racketVel);

		//print_vec("ball out vel = ", ballOutVel);
		//print_vec("racket vel = ",racketVel);
		//print_vec("racket normal = ",racketNormal);

		for (j = 1; j <= N_CART; j++) {
			racketDes[i].x[j] = ballPath[i].x[j];
			racketDes[i].xd[j] = racketVel[j];
			racketDes[i].n[j] = racketNormal[j];
		}
	}

	//printf("Land time = %f\n", landTime);
}

/*
 * Predict the ball for Tpred seconds
 * Using ballflight + bounce model
 *
 * RETURNS TRUE if ball is valid (i.e. bounces only once on the robot side within TPRED seconds)
 *
 */
static int predict_ball_path(SL_Cstate *ballPred, int verbose) {

	int N = TPRED/TSTEP;
	int ball_is_legal = TRUE;
	ballPath = (SL_Cstate*)calloc((N+1), sizeof(SL_Cstate));
	int i,j;

	for (j = 1; j <= N_CART; j++) {
		ballPath[0].x[j] = ballPred->x[j];
		ballPath[0].xd[j] = ballPred->xd[j];
	}

	// predict Tpred seconds into the future
	int bounce = FALSE;
	int numBounces = 0;

	for (i = 1; i <= N; i++) {
		integrate_ball_state(ballPred,TSTEP,&bounce);
		for (j = 1; j <= N_CART; j++) {
			ballPath[i].x[j] = ballPred->x[j];
			ballPath[i].xd[j] = ballPred->xd[j];
		}
		if (bounce) {
			numBounces++;
			bounce = FALSE;
		}
		// multiple bounces are predicted
		if (numBounces > 1) {
			if (verbose)
				printf("Multiple bounces predicted. Not moving\n");
			ball_is_legal = FALSE;
		}
		// no bounce is predicted
		if (ballPath[i].x[_Y_] > dist_to_table && numBounces == 0) {
			if (verbose)
				printf("No bounce predicted. Not moving\n");
			ball_is_legal = FALSE;
		}
	}

	/* revert back to initial ballpred */
	for (j = 1; j <= N_CART; j++) {
		ballPred->x[j] = ballPath[0].x[j];
		ballPred->xd[j] = ballPath[0].xd[j];
	}

	//printf("numBounces = %d\n", numBounces);
	if (ballPath[N].x[_Y_] < dist_to_table && numBounces == 0) {
		if (verbose) {
			printf("Prediction time is not enough to ascertain bounce.\n");
			printf("Returning TRUE (legal ball)\n");
		}
		ball_is_legal = TRUE;
	}

	return ball_is_legal;

}

/*
 * Check for new ball observation. If they differ then we return TRUE.
 *
 * As opposed to the version in filtering (check_new_blob) here we check the
 * differences in ballPred structure [updated through filtering]
 */
static int check_new_ball_obs(SL_Cstate ballNewState) {

	int i,j;
	static int firsttime = TRUE;
	static double last_blob[N_CART+1];

	for (i = 1; i <= N_CART; i++) {
		if (ballNewState.x[i] != last_blob[i]) {
			for (j = 1; j <= N_CART; j++) {
				last_blob[j] = ballNewState.x[j]; // copy to last_blob
			}
			return TRUE;
		}
	}
	return FALSE;

}

/*
 * Predict the ball for Tpred seconds
 * Using ballflight + bounce model
 *
 * RETURNS TRUE if ball is valid (i.e. bounces only once on the robot side within TPRED seconds)
 *
 * The bounce variable is acquired from the Kalman Filter.
 * If it is true (1) then ball is valid unless it is incremented.
 *
 * TODO: do not need to predict for TPRED seconds all the time!
 */
static int predict_ball_path_mpc(SL_Cstate *ballPred, int bounce) {

	int N = TPRED/TSTEP;
	ballPath = (SL_Cstate*)calloc((N+1), sizeof(SL_Cstate));
	int i,j;
	for (j = 1; j <= N_CART; j++) {
		ballPath[0].x[j] = ballPred->x[j];
		ballPath[0].xd[j] = ballPred->xd[j];
	}
	// predict Tpred seconds into the future
	int bouncePred = FALSE;
	int numBounces = bounce;

	for (i = 1; i <= N; i++) {
		integrate_ball_state(ballPred,TSTEP,&bouncePred);
		for (j = 1; j <= N_CART; j++) {
			ballPath[i].x[j] = ballPred->x[j];
			ballPath[i].xd[j] = ballPred->xd[j];
		}
		if (bouncePred) {
			numBounces++;
			bouncePred = FALSE;
		}
		// multiple or no bounces are predicted
		if (ballPath[i].x[_Y_] > dist_to_table && numBounces != 1) {
			printf("No bounce predicted. Not updating\n");
			/* revert back to initial ballpred */
			for (j = 1; j <= N_CART; j++) {
				ballPred->x[j] = ballPath[0].x[j];
				ballPred->xd[j] = ballPath[0].xd[j];
			}
			return FALSE;
		}
	}
	/* revert back to initial ballpred */
	for (j = 1; j <= N_CART; j++) {
		ballPred->x[j] = ballPath[0].x[j];
		ballPred->xd[j] = ballPath[0].xd[j];
	}

	//printf("numBounces = %d\n", numBounces);
	if (ballPath[N].x[_Y_] < dist_to_table && numBounces == 0) {
		printf("Prediction time is not enough to ascertain bounce.\n");
		printf("Returning TRUE (legal ball)\n");
		return TRUE;
	}
	return TRUE;
}

/*
 * Planning based on 3rd order polynomials.
 * Considers the hitting via state when computing the (table tennis) trajectories
 * Computes the polynomial once and then updates the next state.
 * Date: March 2016
 *
 */
int calcNextStateForPoly3rd(SL_DJstate desState[N_DOFS+1], double timeToVia, SL_DJstate viaState[N_DOFS+1],
		            double totalTime, SL_DJstate goalState[N_DOFS+1], double t) {

	int i;
	static int order = 3;
	static double coefInitToVia[4*N_DOFS+1], coefViaToGoal[4*N_DOFS+1];
	static double coef[4*N_DOFS+1]; // 3rd order polynomial coeffs
	static double timeStep = 1./SAMP_FREQ;

	// will be run only once
	// robot is at hitting stage
	if (t <= timeStep + timeStep/2) {
		poly3rd(desState, viaState, timeToVia, coefInitToVia);
		for (i = 1; i <= (order + 1) * N_DOFS; i++)
			coef[i] = coefInitToVia[i];
		//poly3rd(viaState, goalState, totalTime - timeToVia, coefViaToGoal);
	}

	// robot is at return stage
	if (t >= timeToVia + timeStep/2 && t <= timeToVia + 3*timeStep/2) {
		poly3rd(viaState, goalState, totalTime - timeToVia, coefViaToGoal);
		for (i = 1; i <= (order + 1) * N_DOFS; i++)
			coef[i] = coefViaToGoal[i];
	}

	if (t <= totalTime + timeStep/2) {
		if (t > timeToVia + timeStep/2)
			t = t - timeToVia;
		for (i = 1; i <= N_DOFS; i++) {
			desState[i].th = coef[(i-1)*4+1]*t*t*t + coef[(i-1)*4+2]*t*t + coef[(i-1)*4+3]*t + coef[(i-1)*4+4];
			desState[i].thd = 3.0*coef[(i-1)*4+1]*t*t + 2.0*coef[(i-1)*4+2]*t + coef[(i-1)*4+3];
			desState[i].thdd = 6.0*coef[(i-1)*4+1]*t + 2.0*coef[(i-1)*4+2];
		}
	}

	return TRUE;
}

/*
 * Compute the coefficients of a 3rd order polynomial that takes
 * current (mostly initial) state to a goal state in tau time.
 */
int poly3rd(SL_DJstate state[N_DOFS+1], SL_DJstate goal[N_DOFS+1], double T, double *coef) {

	double a1, a2, a3, a4;
	int    i,j;

	coef[0] = 0;

	for (j = 1; j <= N_DOFS; ++j) {

		a1 = 2*(state[j].th - goal[j].th)/(T*T*T) + (goal[j].thd + state[j].thd)/(T*T);
		a2 = 3*(goal[j].th - state[j].th)/(T*T) - (goal[j].thd + 2*state[j].thd)/T;
		a3 = state[j].thd;
	    a4 = state[j].th;

		coef[(j-1)*4+1] = a1;
		coef[(j-1)*4+2] = a2;
		coef[(j-1)*4+3] = a3;
		coef[(j-1)*4+4] = a4;

	}

	return TRUE;
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
void calc_des_racket_vel(Vector velBallIn, Vector velBallOut, Vector normalRacket, Vector velRacket) {

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
void calc_des_racket_normal(Vector bin, Vector bout, Vector normal) {

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
void calc_des_ball_out_vel(SL_Cstate hitPoint, Vector landPoint, double time2land, Vector velOut) {


	static int firsttime = TRUE;
	static double zTable;

	if (firsttime) {
		firsttime = FALSE;
		zTable = floor_level - table_height + ball_radius;
	}

	velOut[_X_] = (landPoint[1] - hitPoint.x[_X_]) / time2land;
	velOut[_Y_] = (landPoint[2] - hitPoint.x[_Y_]) / time2land;
	velOut[_Z_] = (zTable - hitPoint.x[_Z_] + 0.5 * gravity * sqr(time2land)) / time2land;

	//TODO: consider the air drag case
	// hack for now
	velOut[_X_] = 1.1 * velOut[_X_];
	velOut[_Y_] = 1.1 * velOut[_Y_];
	velOut[_Z_] = 1.2 * velOut[_Z_];

}

/*
 * Integrate the ball state dt seconds later.
 * Used for prediction, not for simulation!
 * Using symplectic Euler.
 *
 */
void integrate_ball_state(SL_Cstate *ballState, double dt, int *bounce) {

	int i;
	static Vector spin;
	static int firsttime = TRUE;
	static int resetSim = TRUE; // for ball ground checking to not print anything
	Vector ballvec = my_vector(1,2*N_CART); // temporary ball positions and velocities

	if (firsttime) {
		firsttime = FALSE;
		spin = my_vector(1,N_CART); // constant spin
		init_topspin(spin, PRED_SPIN_MODE);
	}

	// Symplectic Euler for No-Contact-Situation (Flight model)
	//nonlinear_flight_model(ballState);
	spinning_flight_model(ballState, spin);
	symplectic_euler(*ballState,ballvec,dt);
	if (check_ball_table_contact(ballvec,spin,PRED_SPIN_MODE)) {
		*bounce = TRUE;
	}
	check_ball_ground_contact(ballState,ballvec,&resetSim);
	copy_vec_state(ballState,ballvec);
}

/*
 * Utility function used to copy ball state to vector x
 */
void copy_state_vec(SL_Cstate state, Vector x) {

	int i;
	for(i = 1; i <= N_CART; i++) {
		x[i] = state.x[i];
		x[i + N_CART] = state.xd[i];
	}
	x[2*N_CART + 1] = 1;
}

/*
 * Utility function used to copy vector x to ball state
 */
void copy_vec_state(SL_Cstate *state, Vector x) {

	int i;
	for (i = 1; i <= N_CART; i++) {
		state->x[i] = x[i];
		state->xd[i] = x[i + N_CART];
	}
}

/*
 * Set desired landing position to the centre of the table and time to a
 * reasonable value
 *
 */
void set_robot_land_param(Vector ballLand, double *landTime) {

	*landTime = 0.8;

	ballLand[_X_] = 0.0;
	ballLand[_Y_] = dist_to_table - 3*table_length/4; // centre of opponents court
	ballLand[_Z_] = floor_level - table_height;// + ball_radius;

	//print_vec("ball Land = ", ballLand);
}


/*
 * Initialize constant angular velocity (a.k.a. spin) for the spinning ball
 */
void init_topspin(Vector spin, int flag) {

	if (flag) {
		spin[1] = -50*2*PI; // constant 3000 rpm topsin
		// others are zero
	}
	else {
		// do nothing, they are all zero - spinless model
	}
}

/*
 * Spinning ball flight model [including airdrag and Magnus force]
 *
 * Assuming constant spin, not part of the state for convenience
 */
void spinning_flight_model(SL_Cstate *ballState, Vector spin) {

	int i;
	Vector magnus = my_vector(1,N_CART); // acceleration due to magnus force

	nonlinear_flight_model(ballState);
	// add Magnus force
	cross_prod(spin,ballState->xd,magnus);
	vec_mult_scalar(magnus,Clift,magnus);

	for (i = 1; i <= N_CART; i++) {
		ballState->xdd[i] += magnus[i];
	}

}

/*
 * Spin free nonlinear flight model including airdrag (Cdrag as parameter)
 */
void nonlinear_flight_model(SL_Cstate *ballState) {

	static double slack = 0.0001;
	double velBall = sqrt(sqr(ballState->xd[_X_]) + sqr(ballState->xd[_Y_]) + sqr(ballState->xd[_Z_]));

	if (fabs(velBall) > slack) {

		ballState->xdd[_X_] = -ballState->xd[_X_] * Cdrag * velBall;
		ballState->xdd[_Y_] = -ballState->xd[_Y_] * Cdrag * velBall;
		ballState->xdd[_Z_] = gravity - ballState->xd[_Z_] * Cdrag * velBall;
	}
	else {
		ballState->xdd[_X_] = 0.;
		ballState->xdd[_Y_] = 0.;
		ballState->xdd[_Z_] = gravity;
	}
}

/*
 * Crossproduct needed to calculate the magnus force
 */
void cross_prod(Vector spin, double vel[], double out[]) {

	out[1] = spin[2] * vel[3] - spin[3] * vel[2];
	out[2] = spin[3] * vel[1] - spin[1] * vel[3];
	out[3] = spin[1] * vel[2] - spin[2] * vel[1];
}

/*
 * First integrating the accelerations to velocities by dt
 * Then integrating the velocities to positions by dt
 * These (pos and vel) are kept in the ball Vector
 */
void symplectic_euler(SL_Cstate ballState, Vector ball, double dt) {

	ball[4] = ballState.xd[_X_] + ballState.xdd[_X_] * dt;
	ball[5] = ballState.xd[_Y_] + ballState.xdd[_Y_] * dt;
	ball[6] = ballState.xd[_Z_] + ballState.xdd[_Z_] * dt;
	ball[1] = ballState.x[_X_] + ball[4] * dt;
	ball[2] = ballState.x[_Y_] + ball[5] * dt;
	ball[3] = ballState.x[_Z_] + ball[6] * dt;
}

/*
 * Checks contact with racket.
 * If there is contact, then the predicted state will be transformed according to a racket-ball contact model.
 * If racket is going to hit the ball, static hit variable is set to TRUE so that
 * we do not hit it the next time.
 *
 */
int check_ball_racket_contact(SL_Cstate racketState, SL_quat racketOrient, Vector ballvec) {

	int i;
	static Vector racket2Ball; // vector between ball and racket, used for contact modelling
	double velBallAlongRacketNormal;
	double velRacket;
	double normalDistBall2Racket;
	double distBall2Racket;
	double parallelDistBall2Racket;
	static Vector Zdir; // z-direction (0,0,1)
	static Vector racketNormal; // normal of the racket
	static Vector racketVel; // velocity of the racket
	static Vector ballVel; // incoming and outgoing velocities of the ball (in case of contact)

	// flags
	static int hit = FALSE;
	static int firsttime = TRUE;

	if (firsttime) {
		firsttime = FALSE;
		racket2Ball = my_vector(1,N_CART);
		ballVel = my_vector(1,N_CART);
		racketNormal = my_vector(1,N_CART);
		racketVel = my_vector(1,N_CART);
		Zdir = my_vector(1,N_CART);
		Zdir[1] = 0; Zdir[2] = 0; Zdir[3] = 1;
	}

	for (i = 1; i <= N_CART; i++) {
		racketVel[i] = racketState.xd[i];
		racket2Ball[i] = racketState.x[i] - ballvec[i];
		ballVel[i] = ballvec[i + N_CART]; // outgoing,incoming ball vels are initialized to predicted ball velocity
	}

	// calculate the normal of the racket
	quat_vec_mult(racketOrient, Zdir, racketNormal);
	normalDistBall2Racket = vec_mult_inner(racketNormal, racket2Ball);
	distBall2Racket = sqrt(sqr(racket2Ball[1]) + sqr(racket2Ball[2]) + sqr(racket2Ball[3]));
	parallelDistBall2Racket = sqrt(fabs(sqr(distBall2Racket) - sqr(normalDistBall2Racket)));

	// rotate the difference (pos and vel) vectors from global to local coordinates
	// by using inverse quaternion (since we have to take transpose of rotation matrix R)
	for(i = 2; i <= 4; i++) {
		racketOrient.q[i] *= -1;
	}
	rotate_quat(racketOrient, racket2Ball, racket2Ball);

	// when ball passes to the human side reset the hitting condition
	if (ballvec[2] < dist_to_table - 0.5 * table_length)
		hit = FALSE;

	// check for contact with racket
	if ((parallelDistBall2Racket < racket_radius) && (fabs(normalDistBall2Racket) < ball_radius) && (!hit)) {
		hit = TRUE;
		printf("Robot hitting at: x = %f, y = %f, z = %f.\n", ballvec[1], ballvec[2], ballvec[3]);
		// Reflect to back
		racket_contact_model(racketVel, racketNormal, ballVel);
		ballvec[4] = ballVel[_X_];
		ballvec[5] = ballVel[_Y_];
		ballvec[6] = ballVel[_Z_];
		return TRUE;
	}

	return FALSE;
}

/*
 * Update the incoming ball velocity with outgoing ball velocity using MIRROR LAW
 *
 * The racket contact model in vector form is O = I + (1 + eps_R)*N*N'*(V - I)
 * where I is the incoming ball velocity
 *       N is the racket normal
 *       V is the racket velocity
 *       eps_R is the coefficient of restitution of the racket
 *
 * TODO: add a spinning contact model
 *
 */
void racket_contact_model(Vector racketVel, Vector racketNormal, Vector ballVel) {

	static double speed = 0.0;
	static int firsttime = TRUE;
	static Vector diffVel;
	static Vector normalMultSpeed;

	if (firsttime) {
		firsttime = FALSE;
		diffVel = my_vector(1,N_CART);
		normalMultSpeed = my_vector(1,N_CART);
	}

	vec_sub(racketVel, ballVel, diffVel);
	speed = (1 + CRR) * vec_mult_inner(racketNormal, diffVel);
	vec_mult_scalar(racketNormal,speed,normalMultSpeed);
	vec_add(normalMultSpeed,ballVel,ballVel);

}

/*
 * Simple contact model that uses spin if spin mode is turned on.
 * Assuming roll contact instead of slide contact in rebound calculations for simplicity.
 * Spin vector is also modified.
 * Coeff of restitution and friction used.
 */
void table_contact_model(Vector ballvec, Vector spin, int flag) {

	if (flag) { // if spin mode is on ballvec is not a null pointer
		// compute effect on spin
		spin[1] -= (3*(1-CFTX)/(2*ball_radius))*ballvec[5] + (3*(1-CFTX)/2)*spin[1];
		spin[2] += (3*(1-CFTY)/(2*ball_radius))*ballvec[4] - (3*(1-CFTY)/2)*spin[2];
		// in z-direction spin is preserved
		// compute effect on velocity
		ballvec[6] = -CRT * ballvec[6];
		ballvec[5] = CFTY * ballvec[5] - (1-CFTY) * ball_radius * spin[1];
		ballvec[4] = CFTX * ballvec[4] + (1-CFTX) * ball_radius * spin[2];
	}
	else {
		// reflect ball velocity
		ballvec[6] = -CRT * ballvec[6];
		ballvec[5] = CFTY * ballvec[5];
		ballvec[4] = CFTX * ballvec[4];
	}
}


/*
 * Checks contact with ground and zeros the accelerations and velocities
 * Reset sim is also set to TRUE
 */
void check_ball_ground_contact(SL_Cstate *ballState, Vector ballvec, int *resetSim) {

	int i;
	if (ballState->x[_Z_] <= floor_level) {
		for (i = 1; i <= 3; i++)
			ballState->xdd[i] = 0;
		// zero the velocities
		ballvec[4] = 0.0; ballvec[5] = 0.0; ballvec[6] = 0.0;
		// positions remain the same
		ballvec[1] = ballState->x[_X_];
		ballvec[2] = ballState->x[_Y_];
		ballvec[3] = ballState->x[_Z_];
		if (!(*resetSim))
			printf("Contact with floor! \n");
		*resetSim = TRUE;
	}
}

/*
 * Checks contact with net.
 * Curious way to check it: if the net's distance to integrated y-state and distance to current y-state
 * signs do not match, it means that the ball is in contact with the net
 */
void check_ball_net_contact(SL_Cstate ballState, Vector ballvec, int *resetSim) {

	// velocity indexes
	static int _DX_ = 4;
	static int _DY_ = 5;
	static int _DZ_ = 6;
	static int firsttime = TRUE;
	static double contact_table_level;
	static double eps;
	static double eps2;

	if (firsttime) {
		firsttime = FALSE;
		contact_table_level = floor_level - table_height + ball_radius;
	}

	// Check contact with net
	if ((ballvec[_Z_] <= contact_table_level + net_height)
			&& (fabs(ballvec[_X_]) <= table_width/2. + net_overhang)) {
		eps = (dist_to_table -0.5 * table_length) - ballState.x[_Y_];
		eps2 = (dist_to_table -0.5 * table_length) - ballvec[_Y_];

		// If on the other side of the net after integration
		// apply super simplistic model for contact with net
		if ((sign(eps) != sign(eps2))) {

			if (eps >= 0) {
				// Reflect to Front
				ballvec[_DY_] = -(net_restitution) * ballState.xd[_Y_];
				ballvec[_Y_]  = dist_to_table - 0.5 * table_length + (0.5 * net_thickness + ball_radius);
			}
			else {
				// Reflect to Back
				ballvec[_DY_] = -(net_restitution) * ballState.xd[_Y_];
				ballvec[_Y_]  = dist_to_table - 0.5 * table_length - (0.5 * net_thickness + ball_radius);
			}
			printf("Contact with net!\n");
			*resetSim = TRUE;
		}
	}
}

/*
 * Condition to determine if ball hits the table
 * Useful for prediction including a rebound model
 * Useful also in KF/EKF filtering.
 */
int check_ball_table_contact(Vector ballvec, Vector spin, int spin_flag) {

	static double contact_table_level;
	static int firsttime = TRUE;
	static int land = FALSE;

	if (firsttime) {
		contact_table_level = floor_level - table_height + ball_radius;
		firsttime = FALSE;
	}

	// when ball is coming towards the human from robot court reset the landing condition
	if (ballvec[2] > dist_to_table - 0.5 * table_length && ballvec[5] < 0.0)
		land = FALSE;

	// check to see if ball is over the table
	if ((ballvec[2] > dist_to_table - table_length) && (ballvec[2] < dist_to_table) &&
			(fabs(ballvec[1] - table_center) <= table_width/2.0)) {
		// check if the ball hits the table coming from above
		if ((ballvec[3] <= contact_table_level) && ballvec[6] < 0.0) { //revised by Yanlong Huang, Aug 12, 2014.
			// give landing info when ball lands on human side of table coming from robot side
			if ((ballvec[2] < dist_to_table - table_length/2.0) && (ballvec[5] < 0.0) && (!land)) {
				printf("Success! Landing at: x = %f, y = %f.\n",ballvec[1], ballvec[2]);
				land = TRUE;
			}
			table_contact_model(ballvec, spin, spin_flag);
			return TRUE;
		}
	}

	return FALSE;
}

/*
 * Checking if the predicted ball within predTime seconds
 * will bounce once on the robot's court
 *
 * RETURNS FALSE if ball is predicted to bounce less or more than once
 * TODO: check also intersection with robot workspace
 */
int check_legal_ball(SL_Cstate state, double predTime) {

	int i;
	double dt = 0.01;
	int bounce = FALSE;
	int numBounces = 0;
	int numSteps = (int)(predTime / dt);

	for (i = 1; i <= numSteps; i++) {
		integrate_ball_state(&state, dt, &bounce);
		if (bounce) {
			numBounces++;
			bounce = FALSE;
		}
		// multiple bounces are predicted
		if (numBounces > 1) {
			printf("Multiple bounces predicted. Not moving\n");
			return FALSE;
		}
		// no bounce is predicted
		if ((state.x[_Y_] > dist_to_table) && numBounces == 0) {
			printf("No bounce predicted. Not moving\n");
			return FALSE;
		}
	}

	//printf("numBounces = %d\n", numBounces);

	if ((state.x[_Y_] < dist_to_table) && numBounces == 0) {
		printf("Prediction time is not enough to ascertain bounce.\n");
		printf("Returning TRUE (legal ball)\n");
	}

	return TRUE;
}

/*
 * Rotate vector v by a quaternion q to get v_new
 *
 * Note: checked for correctness with the quaternion calculation code and quat_vec_mult!
 * Also checks out with the Rodrigues formula
 *
 */
void rotate_quat(SL_quat q, Vector v, Vector v_new) {

	double a, b, c, d, t2, t3, t4, t5, t6, t7, t8, t9, t10;

	a = q.q[1];
	b = q.q[2];
	c = q.q[3];
	d = q.q[4];
	t2 =   a*b;
	t3 =   a*c;
	t4 =   a*d;
	t5 =  -b*b;
	t6 =   b*c;
	t7 =   b*d;
	t8 =  -c*c;
	t9 =   c*d;
	t10 = -d*d;

	v_new[1] = 2 * ((t8 + t10)*v[1] + (t6 - t4)*v[2] + (t3 + t7)*v[3]) + v[1];
	v_new[2] = 2 * ((t4 + t6)*v[1] + (t5 + t10)*v[2] + (t9 - t2)*v[3]) + v[2];
	v_new[3] = 2 * ((t7 - t3)*v[1] + (t2 + t9)*v[2] + (t5 + t8)*v[3]) + v[3];
}

/*
 * Multiplies the quaternion with a vector to get a vector output (3rd parameter).
 * First forms the rotation matrix that corresponds to the quaternion and then performs matrix multiplication.
 * Note: checked for correctness
 *
 * The rotation matrix version of rotate_quat. rotate_quat should be faster!
 *
 *
 */
void quat_vec_mult(SL_quat q, Vector in, Vector out) {

	static Matrix Q;
	static int firsttime = TRUE;

	double q0, q1, q2, q3;

	if(firsttime) {
		firsttime = FALSE;
		Q = my_matrix(1, N_CART, 1, N_CART);
	}

	q0 = q.q[_Q0_];
	q1 = q.q[_Q1_];
	q2 = q.q[_Q2_];
	q3 = q.q[_Q3_];

	Q[1][1] = 2*sqr(q0) - 1 + 2*sqr(q1);
	Q[1][2] = 2*q1*q2 - 2*q0*q3;
	Q[1][3] = 2*q1*q3 + 2*q0*q2;
	Q[2][1] = 2*q1*q2 + 2*q0*q3;
	Q[2][2] = 2*sqr(q0) - 1 + 2*sqr(q2);
	Q[2][3] = 2*q2*q3 - 2*q0*q1;
	Q[3][1] = 2*q1*q3 - 2*q0*q2;
	Q[3][2] = 2*q2*q3 + 2*q0*q1;
	Q[3][3] = 2*sqr(q0) - 1 + 2*sqr(q3);

	mat_vec_mult(Q,in,out);
}

