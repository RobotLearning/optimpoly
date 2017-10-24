#ifndef table_tennis_common_H_
#define table_tennis_common_H_

// SL general includes of system headers
#include "SL_system_headers.h"

/* private includes */
#include "SL.h"
#include "SL_common.h"
#include "SL_user.h"
#include "utility.h"
#include "SL_tasks.h"
#include "SL_task_servo.h"
#include "SL_man.h"
#include "SL_dynamics.h"
#include "SL_collect_data.h"
#include "SL_shared_memory.h"
#include "SL_kinematics.h"
#include "SL_filters.h"
#include "time.h"
#include "table.h"

/**************************** defines ***********************/

#define JINV_WEIGHTS_FILE  "JinvWeights"
#define SAMP_FREQ 500 // sampling frequency, frequency of the robot in this case
#define MAX_ACC 200
#define MAX_VEL 200
#define SLACK 0.01 // slack for checking joint limits
#define DT 1.0/SAMP_FREQ // sampling frequency, frequency of the trajectory

typedef struct {
	int status; //!< was ball detected reliably in cameras
	double pos[N_CART]; //!< ball center cartesian positions from cameras 1 and 2(after calibration)
} blob_state;

/************************************** global variables *****************/

/* Racket state and orientation */
extern SL_Cstate  racket_state;
extern SL_quat    racket_orient;

/* Ball state in simulation */
extern SL_Cstate  sim_ball_state;

// ball status variable used for prediction
extern SL_Cstate ballPred;

// variable for trajectory generation and planning
extern double trj_time;

// initial state
extern SL_DJstate init_joint_state[N_DOFS+1];

/* flag for resetting simulation, bounce, etc. */
extern int simulation;
extern int moving;
extern int bounce;

// flags for table tennis players
extern int SIM_SPIN_MODE;
extern int SIM_FILTER_MODE;
extern int PRED_SPIN_MODE;

/************************************* constants for simulation ***********************/

/* Contact Coefficients */
/* coefficient of restitution for the table (i.e. rebound z-velocity multiplier) */
static const double CRT = 0.88;//0.91;cannon//0.8; human //0.88 das hatten wir mal experimentell bestimmt, 0.92 fittet aber besser in die beobachteten Daten
/* coefficient of table contact model on Y-direction */
static const double CFTY = 0.72;//0.89;cannon//0.78; human //0.6 is eigentlich das richtige, aber da immer ein wenig spin drin sein wird, passt das am besten in die beobachteten Daten
/* coefficient of table contact model on X-direction */
static const double CFTX = 0.68; //0.74
/* coefficent of restitution for racket */
static const double CRR = 0.78;//0.9;cannon//0.78;human//0.78;

/* Air drag coefficient */
static const double Cdrag = 0.1414;

/* for simulating different gravities */
static const double intern_gravity = 9.802;

static const double Clift = 0.001; // coefficient of lift for the magnus force

/************************************ functions *******************/

/* functions for simulating ball */
int simulate_ball(SL_Cstate *bstate, SL_Cstate *rstate, SL_quat *rorient, int *contact);
void nonlinear_flight_model(SL_Cstate *ballState);
void symplectic_euler(SL_Cstate ballState, Vector ball, double dt);
int check_ball_opponent_contact(Vector ballvec);
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
void cross_prod(Vector spin, double vel[], double out[]); // needed for magnus force calculation

/* Functions for drawing ball */
void display_sim_ball(void);
void display_pred_ball(void);
void reset_sim_ball();
void init_ball_cannon(SL_Cstate *ballcannon);
void update_blobs(SL_Cstate *ballstate);

/* functions for racket calculations */
void set_robot_land_param(Vector ballLand, double *landTime);
void set_opp_land_param(Vector ballLand, double *landTime, int rand);
void calc_racket(SL_Cstate *rstate, SL_quat *rorient, SL_Cstate cart_state, SL_quat cart_orient);
void calc_des_racket_vel(Vector bin, Vector bout, Vector normal, Vector vracket);
void calc_des_racket_normal(Vector bin, Vector bout, Vector normal);
void calc_des_ball_out_vel(SL_Cstate origin, Vector goal, double time2land, Vector velOut);
void integrate_ball_state(SL_Cstate *state, double dt, int *bounce);

/* other utility functions related to robot dynamics */
long get_time(void);
int rms_traj_error(Matrix traj_pos_act, Matrix traj_vel_act, Matrix traj_pos, Matrix traj_vel);
int check_joint_limits(SL_Jstate *act, double th);
int check_des_joint_limits(SL_DJstate *des, double th);
int collision_detection(SL_Cstate racket_state);
//int goto_posture(double* q);

/* utility functions */
void copy_state_vec(SL_Cstate state,Vector x);
void copy_vec_state(SL_Cstate *state, Vector x);

/* printing information */
int save_sim_data();

#endif /* table_tennis_common_H_ */
