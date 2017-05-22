/*
 * player.h
 *
 *  Created on: Feb 6, 2017
 *      Author: okoc
 */

#ifndef PLAYER_H_
#define PLAYER_H_

#define PI 3.14

//#ifdef __cplusplus
//extern "C" {
//#endif

extern void play();

//#ifdef __cplusplus
//}
//#endif


#ifdef __cplusplus
void calc_ball_state(int *init_filter, int *bounce, int *num_obs);

void calc_param(Vector ballLand, double landTime,
		nlopt_thread_data *optim_data, int *moving, int bounce, int numobs);
static void calc_optim_param(Vector ballLand, double landTime, int *moving,
		                     nlopt_thread_data *optim_data);
static void calc_mpc_param(Vector ballLand, double landTime, nlopt_thread_data *optim_data,
		                   int bounce, int numobs);
static int predict_ball_path(SL_Cstate *ballPred, int verbose);
void calc_next_state(nlopt_thread_data optim_data, int *moving, double returnToRestTime);
void calc_racket_strategy(Vector ballLand, double landTime);
static int check_new_ball_obs(SL_Cstate ballNewState);
static int predict_ball_path_mpc(SL_Cstate *ballPred, int bounce);

int poly3rd(SL_DJstate state[N_DOFS+1], SL_DJstate goal[N_DOFS+1], double T, double *coef);
int calcNextStateForPoly3rd(SL_DJstate state[N_DOFS+1], double timeToVia, SL_DJstate viaState[N_DOFS+1],
		            double totalTime, SL_DJstate goalState[N_DOFS+1], double t);

void calc_des_racket_vel(Vector bin, Vector bout, Vector normal, Vector vracket);
void calc_des_racket_normal(Vector bin, Vector bout, Vector normal);
void calc_des_ball_out_vel(SL_Cstate origin, Vector goal, double time2land, Vector velOut);
void set_robot_land_param(Vector ballLand, double *landTime);

void integrate_ball_state(SL_Cstate *state, double dt, int *bounce);

/* utility functions */
void copy_state_vec(SL_Cstate state,Vector x);
void copy_vec_state(SL_Cstate *state, Vector x);

/* functions for simulating ball */
void symplectic_euler(SL_Cstate ballState, Vector ball, double dt);
void check_ball_ground_contact(SL_Cstate *ballState, Vector ballvec, int *resetSim);
int check_ball_racket_contact(SL_Cstate racketState, SL_quat racketOrient, Vector ballvec);
int check_legal_ball(SL_Cstate state, double predTime);
int check_ball_table_contact(Vector ballvec, Vector spin, int spin_flag);
void check_ball_net_contact(SL_Cstate ballState, Vector ballvec, int *resetSim);
void racket_contact_model(Vector racketVel, Vector racketNormal, Vector ballVel);
void table_contact_model(Vector ballvec, Vector spin, int spin_flag);

/* spinning ball functions */
void init_topspin(Vector spin, int flag);
void spinning_flight_model(SL_Cstate *ballState, Vector spin);
void nonlinear_flight_model(SL_Cstate *ballState);
void cross_prod(Vector spin, double vel[], double out[]); // needed for magnus force calculation

void rotate_quat(SL_quat q, Vector v, Vector v_new);
void quat_vec_mult(SL_quat q, Vector in, Vector out);

#endif

#endif /* PLAYER_H_ */
