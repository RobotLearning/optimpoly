/*
 * Table tennis common includes all the
 * functions used for modelling table tennis interactions between
 * ball, table, racket (and robots, including possible opponents).
 *
 * Most of them are used in simulation mode, such as
 * the simulate_ball function.
 *
 * We support in simulation opponent modelling, simple spin modelling, etc.
 *
 * TODO: simplify inverse kinematics functions for table tennis, maybe push
 * to inverse_kinematics.c
 *
 *
 */

#include "table_tennis_common.h"

/*
 * Global variables to be shared across tasks
 *
 */
SL_Cstate  racket_state; /* Racket state and orientation */
SL_quat    racket_orient;
SL_Cstate  sim_ball_state; /* Ball state in simulation */
SL_Cstate ballPred; // ball status variable used for prediction
double trj_time; // variable for trajectory generation and planning
SL_DJstate init_joint_state[N_DOFS+1]; // initial state
int simulation = TRUE;
int moving = FALSE;  //indicate moving/waiting
int bounce = FALSE;

int SIM_SPIN_MODE = FALSE; // turn on spin based ball in integration in simulation
int SIM_FILTER_MODE = FALSE; // turn on ball filtering in simulation mode
int PRED_SPIN_MODE = FALSE; // turn on prediction with a spin model

// static variables that we dont make global
static int OPPONENT_MODE = FALSE; // turn on opponent player mode

/*
 * Draws the simulated ball and the racket every 2 ms for simulation
 */
void display_sim_ball(void) {

	double pos[14];

	pos[4] = racket_state.x[_X_];
	pos[5] = racket_state.x[_Y_];
	pos[6] = racket_state.x[_Z_];
	pos[7] = racket_orient.q[_Q0_];
	pos[8] = racket_orient.q[_Q1_];
	pos[9] = racket_orient.q[_Q2_];
	pos[10] = racket_orient.q[_Q3_];

	//ball simulated
	pos[11] = sim_ball_state.x[_X_];
	pos[12] = sim_ball_state.x[_Y_];
	pos[13] = sim_ball_state.x[_Z_];

	sendUserGraphics("table_tennis", &(pos[0]), 14*sizeof(double));
}

/*
 * Draws the predicted ball and the racket every 2 ms for simulation
 */
void display_pred_ball(void) {

	double pos[14];

	pos[4] = racket_state.x[_X_];
	pos[5] = racket_state.x[_Y_];
	pos[6] = racket_state.x[_Z_];
	pos[7] = racket_orient.q[_Q0_];
	pos[8] = racket_orient.q[_Q1_];
	pos[9] = racket_orient.q[_Q2_];
	pos[10] = racket_orient.q[_Q3_];
 
	//ball predicted
	pos[11] = ballPred.x[_X_];
	pos[12] = ballPred.x[_Y_];
	pos[13] = ballPred.x[_Z_];

	sendUserGraphics("table_tennis", &(pos[0]), 14*sizeof(double));
}

/*
 *
 * Reset the simulated ball state.
 * Adds noise to the initial ball launch velocity.
 * Unless you set the seed to a reasonable value - Okan. *
 *
 */
void reset_sim_ball() {

	int i;
	static int firsttime = TRUE;
	static SL_Cstate ball_cannon;
	static int seed = 2;
	static double std = 0.1;
	//srand(time(0));

	if (firsttime) {
		srand(seed);
		firsttime = FALSE;
		bzero((char *)&(ball_cannon), sizeof(ball_cannon));
		init_ball_cannon(&ball_cannon);
	}

	sim_ball_state.x[_X_]  = ball_cannon.x[_X_] + gaussian(0,std);
	sim_ball_state.x[_Y_]  = ball_cannon.x[_Y_] + gaussian(0,std);
	sim_ball_state.x[_Z_]  = ball_cannon.x[_Z_] + gaussian(0,std);

	//double rand_val = ((double)rand())/RAND_MAX - 0.5;
	sim_ball_state.xd[_X_] = -1.08 + gaussian(0,std);
	sim_ball_state.xd[_Y_] = 4.80 + gaussian(0,std);
	sim_ball_state.xd[_Z_] = 3.84 + gaussian(0,std);

	/*printf("Init pos:\n");
	printf("x = [%.2f,%.2f,%.2f] \t", sim_ball_state.x[_X_],sim_ball_state.x[_Y_],
			  	  	  	  	  	  	  sim_ball_state.x[_Z_]);
	printf("xd = [%.2f,%.2f,%.2f]. \n", sim_ball_state.xd[_X_],sim_ball_state.xd[_Y_],
									  sim_ball_state.xd[_Z_]);*/

}

/*
 * Initialize ball cannon somewhere behind the table
 */
void init_ball_cannon(SL_Cstate* ball_cannon) {

	ball_cannon->x[_X_] = table_center + 0.4;
	ball_cannon->x[_Y_] = dist_to_table - table_length - 0.2; // since dist_to_table is negative
	ball_cannon->x[_Z_] = floor_level - table_height + 0.15;
}

/*
 * Date	: August 2007
 * Modified: July-August 2016
 *
 * Simulates a ball and integrates for the next ball states.
 * Checking contacts with environment, i.e. racket, net, table, ground.
 * TODO: implement RK4
 *
 * Takes around 1mu sec to run
 *
 */
int simulate_ball(SL_Cstate *bstate, SL_Cstate *rstate, SL_quat *rorient, int *resetSim) {

	static Vector spin;
	static int firsttime = TRUE;
	Vector ballvec = my_vector(1,2*N_CART); // temporary ball positions and velocities

	if (firsttime) {
		firsttime = FALSE;
		spin = my_vector(1,N_CART); // constant spin
		init_topspin(spin,SIM_SPIN_MODE);
	}

	// Symplectic Euler for No-Contact-Situation (Flight model)
	// nonlinear_flight_model(bstate);
	spinning_flight_model(bstate, spin);
	symplectic_euler(*bstate, ballvec, DT);

	// Check contact to table
	check_ball_table_contact(ballvec,spin,SIM_SPIN_MODE);
	// Check contact with net
	check_ball_net_contact(*bstate, ballvec, resetSim);
	// Check contact with racket
	check_ball_racket_contact(*rstate,*rorient,ballvec);
	// Check contact with opponent
	check_ball_opponent_contact(ballvec);
	// Check if it hits the ground...
	check_ball_ground_contact(bstate,ballvec,resetSim);

	// Pass the computed ball variables to bstate
	copy_vec_state(bstate,ballvec);
	//update_blobs(bstate); // draw new ball state

	return TRUE;
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
 * Spin free nonlinear flight model including airdrag (Cdrag as parameter)
 */
void nonlinear_flight_model(SL_Cstate *ballState) {

	static double slack = 0.0001;
	double velBall = sqrt(sqr(ballState->xd[_X_]) + sqr(ballState->xd[_Y_]) + sqr(ballState->xd[_Z_]));

	if (fabs(velBall) > slack) {

		ballState->xdd[_X_] = -ballState->xd[_X_] * Cdrag * velBall;
		ballState->xdd[_Y_] = -ballState->xd[_Y_] * Cdrag * velBall;
		ballState->xdd[_Z_] = -intern_gravity - ballState->xd[_Z_] * Cdrag * velBall;
	}
	else {
		ballState->xdd[_X_] = 0.;
		ballState->xdd[_Y_] = 0.;
		ballState->xdd[_Z_] = -intern_gravity;
	}
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
 * Check contact with opponent.
 *
 * If the ball is crossing the 2nd half of opponent's court
 * Then the opponent picks a random point on the incoming ball
 * When the ball goes close to that point
 * It will be reflected with a velocity calculated to land the point
 * at a desired point at a desired landing time
 *
 */
int check_ball_opponent_contact(Vector ballvec) {

	int i;
	static int firsttime = TRUE;
	static double landTime;
	static Vector ballLand;
	static double std = 0.1; // standard deviation for random ball velocity
	static double opp_hit_begin;
	static double opp_hit_end;
	static int opp_decided = FALSE;
	Vector velOut = my_vector(1,N_CART);
	static SL_Cstate hitPoint;
	double opp_hit_point;
	double rand_val;

	if (firsttime) {
		firsttime = FALSE;
		opp_hit_begin = dist_to_table - (3.6*table_length)/4.0;
		opp_hit_end = dist_to_table - 5*table_length/4.0;
		ballLand = my_vector(1,2);
		bzero((char *)&(hitPoint), sizeof(hitPoint));
	}

	if (OPPONENT_MODE) {
		// TODO: there should be a check for legal ball here
		if ((ballvec[2] < opp_hit_begin) && (ballvec[2] > opp_hit_end) &&
		    (ballvec[5] < -0.0) && (!opp_decided)) {

			opp_decided = TRUE;
			// pick a random hitting point
			rand_val = ((double)rand())/RAND_MAX;
			opp_hit_point = opp_hit_begin + rand_val *
					(opp_hit_end - opp_hit_begin);
			// pick a random landing point
			set_opp_land_param(ballLand, &landTime, TRUE);
		}

		if (opp_decided && (ballvec[2] - opp_hit_point) > 0.0) {
			opp_decided = FALSE;
			hitPoint.x[1] = ballvec[1];
			hitPoint.x[2] = ballvec[2];
			hitPoint.x[3] = ballvec[3];
			calc_des_ball_out_vel(hitPoint, ballLand, landTime, velOut);
			// reflect the ball with a random velocity
			ballvec[4] = velOut[1];
			ballvec[5] = velOut[2];
			ballvec[6] = velOut[3];
		}

	}
	return TRUE;
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
	static Vector Z; // z-direction (0,0,1)
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
		Z = my_vector(1,N_CART);
		Z[1] = 0; Z[2] = 0; Z[3] = 1;
	}

	for (i = 1; i <= N_CART; i++) {
		racketVel[i] = racketState.xd[i];
		racket2Ball[i] = racketState.x[i] - ballvec[i];
		ballVel[i] = ballvec[i + N_CART]; // outgoing,incoming ball vels are initialized to predicted ball velocity
	}

	// calculate the normal of the racket
	quat_vec_mult(racketOrient, Z, racketNormal);
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
		if ((macro_sign(eps) != macro_sign(eps2))) {

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
 * Send to graphics the integrated state as the raw blobs in simulation
 * In the WAM computer blobs structure is updated directly
 */
void update_blobs(SL_Cstate *ballstate) {

	static int firsttime = TRUE;
	static char hostname[1024];
	static int last_frame_counter = -999;

	if (firsttime) {
		firsttime = FALSE;
		hostname[1023] = '\0';
		gethostname(hostname,1023);
	}

	/*
	 * For the robot computer we have to update blobs directly
	 */
	if (strcmp(hostname,"wam") == 0) {
		blobs[BLOB1].status = TRUE;
		blobs[BLOB1].blob.x[_X_] = ballstate->x[_X_];
		blobs[BLOB1].blob.x[_Y_] = ballstate->x[_Y_];
		blobs[BLOB1].blob.x[_Z_] = ballstate->x[_Z_];
	}
	else {
		raw_blobs[BLOB1].status = TRUE;
		raw_blobs[BLOB1].x[_X_] = ballstate->x[_X_];
		raw_blobs[BLOB1].x[_Y_] = ballstate->x[_Y_];
		raw_blobs[BLOB1].x[_Z_] = ballstate->x[_Z_];
	}

	// send raw_blobs
	if (last_frame_counter != frame_counter) {
		send_raw_blobs();
		last_frame_counter = frame_counter;
	}

}

/*
 * Set desired landing position to a random position on the robot court
 *
 * Turn on randomization only if randomize is TRUE
 */
void set_opp_land_param(Vector ballLand, double *landTime, int randomize) {

	static int firsttime = TRUE;
	static double robot_court_land_x, robot_court_land_y;
	static double std = 0.1;
	double rand_time = 0.0;
	Vector rand_land = my_vector(1,2);

	if (firsttime) {
		firsttime = FALSE;
		robot_court_land_x = -0.3;
		robot_court_land_y = dist_to_table - 0.4;
	}

	if (randomize) {
		rand_time = 0.1 * ((double)rand())/RAND_MAX;
		rand_land[1] = gaussian(0,std);
		rand_land[2] = gaussian(0,std);
	}

	*landTime = 0.6 + rand_time;
	ballLand[1] = robot_court_land_x + rand_land[1];
	ballLand[2] = robot_court_land_y + rand_land[2];

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
 * Calculate racket cartesian states and orientations
 * using cartesian state and orientation of endeffector
 * INPUTS:
 * 			cart_state and cart_orient
 * OUTPUTS:
 * 			rstate and rorient
 */
void calc_racket(SL_Cstate *racket_state, SL_quat *racket_orient,
		SL_Cstate cart_state, SL_quat cart_orient) {
	int i;
	static int firsttime = TRUE;
	static SL_quat cup_rot;

	if (firsttime) {
		cup_rot.q[_Q0_] = cos(PI/4);
		cup_rot.q[_Q1_] = -sin(PI/4);//0;
		cup_rot.q[_Q2_] = 0;//sin(PI/4);
		cup_rot.q[_Q3_] = 0;

		firsttime = FALSE;
	}

	// copy cartesian state
	for (i = 1; i <= N_CART; i++) {
		racket_state->x[i] = cart_state.x[i];
		racket_state->xd[i] = cart_state.xd[i];
		racket_state->xdd[i] = cart_state.xdd[i];
	}

	// multiply orientations
	mult_quat(&cart_orient, &cup_rot, racket_orient);

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
	velOut[_Z_] = (zTable - hitPoint.x[_Z_] + 0.5 * intern_gravity * sqr(time2land)) / time2land;

	//TODO: consider the air drag case
	// hack for now
	velOut[_X_] = 1.1 * velOut[_X_];
	velOut[_Y_] = 1.1 * velOut[_Y_];
	velOut[_Z_] = 1.2 * velOut[_Z_];

}

/*
 * Used to go to initial posture
 *
 */
int goto_posture(double* q) {

	int i;
	/* go to start posture */
	static SL_DJstate target[N_DOFS+1];
	bzero((char *)&(target[1]), n_dofs * sizeof(target[1]));
	for (i = 1; i <= n_dofs; ++i) {
		target[i] = joint_default_state[i];
		target[i].th = q[i];
	}

	if (!go_target_wait_ID(target)) {
		return FALSE;
	}
	else {
		if (!go_target_wait(target))
			return FALSE;
	}
	return TRUE;
}

/*****
 *
 * Compute RMS error. 
 * Adds also the velocities error.
 *
 * Author: Okan Koc
 */
int rms_traj_error(Matrix traj_pos_act, Matrix traj_vel_act, Matrix traj_pos, Matrix traj_vel) {

	// extract errors

	int i,j;
	static Vector Q;
	static int firsttime = TRUE;
	int trj_len = traj_pos_act[0][0];
	//printf("Trajectory length: %d\n",trj_len);
	int len = trj_len - 1;
	static Vector error;
	double sqnorm = 0.0;
	double finCost = 0.0;

	if (firsttime) {
		firsttime = FALSE;
		error = my_vector(1, len*2*N_DOFS);
		Q = my_vector(1,2*N_DOFS);
		for (i = 1; i <= N_DOFS; i++) {
			Q[i] = 1.0;
		}
		//Q[7] = 0.0;
		for (i = N_DOFS+1; i <= 2*N_DOFS; i++) {
			Q[i] = 1.0;
		}
		//Q[14] = 0.0;
	}

	for (i = 1; i <= len; i++) {
		for (j = 1; j <= N_DOFS; j++) {
			error[2*N_DOFS*(i-1) + j] = traj_pos_act[i+1][j] - traj_pos[i+1][j];
			error[2*N_DOFS*(i-1) + N_DOFS + j] = traj_vel_act[i+1][j] - traj_vel[i+1][j];
		}
	}

	for (i = 1; i <= len; i++) {
		for (j = 1; j <= 2*N_DOFS; j++) {
			sqnorm = sqnorm + error[2*N_DOFS*(i-1) + j]*Q[j]*error[2*N_DOFS*(i-1) + j];
		}
	}

	for (j = 1; j <= 2*N_DOFS; j++) {
		finCost = finCost + error[2*N_DOFS*(len-1) + j]*Q[j]*error[2*N_DOFS*(len-1) + j];
	}

	//dump_vec(error,2*N_DOFS*(traj_len-1),"error_vec");
	//printf("SSE error on trajectory: %f\n", sqnorm);
	printf("RMS error on trajectory: %f\n", sqrt(sqnorm/len));
	printf("Squareroot of final cost: %f\n", sqrt(finCost));

	return 1;
}

/*
 * Return time of day as micro seconds
 */
long get_time() {
	struct timeval tv;
	if (gettimeofday(&tv, (struct timezone *)0) == 0)
		return (tv.tv_sec*1000*1000 + tv.tv_usec);  //us

	return 0.;
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
 * Save ball sim data, pred data, and joint data to a new file.
 */
int save_sim_data() {

	static double timeLabel;
	static int num, j;
	static double timeStart;
	static char dirName[] = "saveData";
	static char fileName[100];
	static int firsttime = TRUE;
	static double timeElapsed = 0.0;
	static time_t date;
	static char dateString[100];
	static struct tm tm;

	if (firsttime == TRUE) {

		firsttime = FALSE;
		timeStart = get_time();
		date = time(NULL);
		tm = *localtime(&date);
		sprintf(dateString, "_%d-%d-%d_%d:%d:%d", tm.tm_mday, tm.tm_mon + 1, tm.tm_year - 100,
				tm.tm_hour, tm.tm_min, tm.tm_sec);
		sprintf(fileName, "%s//sim_data%s.txt", dirName, dateString);
		// do not overwrite each time sl is started
		//if(access(fileName, 0) != -1)
		//	remove(fileName);
	}
	timeElapsed = (get_time() - timeStart) / 1e3;

	FILE * fid;
	fid = fopen(fileName, "at+");
	fprintf(fid, "%.4f  ",timeElapsed);
	fprintf(fid, "%.4f  %.4f  %.4f  ", sim_ball_state.x[1], sim_ball_state.x[2], sim_ball_state.x[3]);
	fprintf(fid, "%.4f  %.4f  %.4f  ", sim_ball_state.xd[1], sim_ball_state.xd[2], sim_ball_state.xd[3]);
	fprintf(fid, "%.4f  %.4f  %.4f  ",ballPred.x[1], ballPred.x[2],ballPred.x[3]);
	fprintf(fid, "%.4f  %.4f  %.4f  ",ballPred.xd[1], ballPred.xd[2],ballPred.xd[3]);
	for (j = 1; j <= 7; j++)
		fprintf(fid, "%.4f  ",joint_des_state[j].th);
	for (j = 1; j <= 7; j++)
		fprintf(fid, "%.4f  ",joint_des_state[j].thd);
	for (j = 1; j <= 7; j++)
		fprintf(fid, "%.4f  ",joint_state[j].th);
	for (j = 1; j <= 7; j++)
		fprintf(fid, "%.4f  ",joint_state[j].thd);
	fprintf(fid,"\n");
	fclose(fid);

	return TRUE;
}

/*
 * Check for limits of a desired joint q_des,qd_des,qdd_des generated by a motion planner
 * Second argument is the slack
 */
int check_des_joint_limits(SL_DJstate *des, double th) {

	int i;
	int flag = TRUE;
	double acc_limit = MAX_ACC;
	double vel_limit = MAX_VEL;

	for (i = 1; i <= N_DOFS; ++i) {
		if (des[i].th + th > joint_range[i][MAX_THETA]) {
			printf("qdes[%d] = %f > max qdes[%d] = %f.\n", i, des[i].th, i, joint_range[i][MAX_THETA] - th);
			flag = FALSE;
		}
		if (des[i].th - th < joint_range[i][MIN_THETA]) {
			printf("qdes[%d] = %f < min qdes[%d] = %f.\n", i, des[i].th, i, joint_range[i][MIN_THETA] + th);
			flag = FALSE;
		}
		if (des[i].thd > vel_limit) {
			printf("qd_des[%d] = %f > max qd_des[%d] = %f.\n", i, des[i].thd, i, vel_limit);
			flag = FALSE;
		}
		if (des[i].thd < -vel_limit) {
			printf("qd_des[%d] = %f < min qd_des[%d] = %f.\n", i, des[i].thd, i, -vel_limit);
			flag = FALSE;
		}
		if (des[i].thdd > acc_limit) {
			printf("qdd_des[%d] = %f > max qdd_des[%d] = %f.\n", i, des[i].thdd, i, acc_limit);
			flag = FALSE;
		}
		if (des[i].thdd < -acc_limit) {
			printf("qdd_des[%d] = %f < min qdd_des[%d] = %f.\n", i, des[i].thdd, i, -acc_limit);
			flag = FALSE;
		}
	}
	return flag;
}

/*
 * Checks for the limits of current states!
 * (as opposed to checkDesJointLimits which only checks for the desired state.)
 */
int check_joint_limits(SL_Jstate *cur, double th) {

	int i;
	int flag = TRUE;
	double acc_limit = MAX_ACC;
	double vel_limit = MAX_VEL;

	for (i = 1; i <= N_DOFS; ++i) {
		if (cur[i].th + th > joint_range[i][MAX_THETA]) {
			printf("q[%d] = %f > max q[%d] = %f.\n", i, cur[i].th, i, joint_range[i][MAX_THETA] - th);
			flag = FALSE;
		}
		if (cur[i].th - th < joint_range[i][MIN_THETA]) {
			printf("q[%d] = %f < min q[%d] = %f.\n", i, cur[i].th, i, joint_range[i][MIN_THETA] + th);
			flag = FALSE;
		}
		/*if (cur[i].thd > vel_limit) {
			printf("qd[%d] = %f > max qd[%d] = %f.\n", i, cur[i].thd, i, vel_limit);
			flag = FALSE;
		}
		if (cur[i].thd < -vel_limit) {
			printf("qd[%d] = %f < min qd[%d] = %f.\n", i, cur[i].thd, i, -vel_limit);
			flag = FALSE;
		}
		if (cur[i].thdd > acc_limit) {
			printf("qdd[%d] = %f > max qdd[%d] = %f.\n", i, cur[i].thdd, i, acc_limit);
			flag = FALSE;
		}
		if (cur[i].thdd < -acc_limit) {
			printf("qdd[%d] = %f < min qdd[%d] = %f.\n", i, cur[i].thdd, i, -acc_limit);
			flag = FALSE;
		}*/
	}
	return flag;
}

/*
 * Detect collision of racket with table
 */
int collision_detection(SL_Cstate racket_state) {
	double y,z;

	y = dist_to_table - racket_state.x[_Y_];
	z = (floor_level - table_height) - racket_state.x[_Z_];

	if (sqrt(sqr(y) + sqr(z)) <= racket_radius) {
		return TRUE;
	}

	return FALSE;
}
