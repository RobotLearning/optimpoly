/*
 * optim.h
 *
 *  Created on: Aug 9, 2016
 *      Author: okan
 */

#ifndef OPTIM_H_
#define OPTIM_H_

#include "math.h"
#include <nlopt.h>
#include <pthread.h>
#include "string.h"
#include "kinematics.h"

/* DEFINES */

#define MAX_ACC 200
#define MAX_VEL 200
#define EQ_CONSTR_DIM 3*N_CART
#define OPTIM_DIM 2*N_DOFS + 1
#define INEQ_CONSTR_DIM 2*N_DOFS + 2*N_DOFS // both strike and returning trajectories, min and max
#define TSTEP 0.002
#define TPRED 2.0
#define TIME2RETURN 1.0
#define SLACK 0.01 // slack for checking joint limits
#define SAMP_FREQ 500 // sampling frequency, frequency of the robot in this case
#define DT 1.0/SAMP_FREQ // sampling frequency, frequency of the trajectory
#define BLOB1 1
#define MAX_BLOBS 3

/************************************** global variables *****************/

extern double        joint_range[N_DOFS+1][3+1];     /* various info on joint limits */
extern SL_Jstate     joint_state[N_DOFS+1];          /* current states */
extern SL_DJstate    joint_des_state[N_DOFS+1];      /* desired states */

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

extern SL_VisionBlob blobs[MAX_BLOBS+1];             /* blob info from vision */

/* new structures defined to be used globally */

typedef struct {
	SL_DJstate target[N_DOFS+1];
	double hitTime;
} nlopt_thread_data;

typedef struct { /* end-effector parameters */
  double   x[N_CART+1];   /* end position of endeffector in global coordinates*/
  double   xd[N_CART+1];  /* velocity of endeffector */
  double   n[N_CART+1];   /* orientation of the tool as a normal vector */
} Racket;

/* Global data that cp_task also uses */

extern SL_Cstate *ballPath; // predicted ball pos and vel values for T_pred time seconds
extern Racket *racketDes; // racket strategy
extern SL_DJstate current_joint_state[N_DOFS+1]; // for correction with MPC
extern int busy; // to check between threads whether optimization is running

/* Multi-threading */
void* mthread_nlopt_optim(void *arg);
int mthread_sync_optim(nlopt_thread_data data, SL_DJstate target[], double *hitTime);
int check_thread_termination(nlopt_thread_data data, SL_DJstate target[], double hitTime);

/* Optimization functions */
void nlopt_optim_fixed_run(SL_DJstate target[], double* hitTime);
double costfunc(unsigned n, const double *x, double *grad, void *my_func_params);
double const_costfunc(unsigned n, const double *x, double *grad, void *my_func_params);
void kinematics_eq_constr(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data);
void joint_lim_ineq_constr(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data);
void calc_strike_poly_coeff(double *a1, double *a2, const double *x);
void calc_return_poly_coeff(double *a1, double *a2, const double *x);
void calc_strike_ext_cand(double *a1, double *a2, double T, double *joint_max_cand, double *joint_min_cand);
void calc_return_ext_cand(double *a1, double *a2, const double *x, double *joint_max_cand, double *joint_min_cand);
void set_bounds(double *lb, double *ub, double slack);
void init_soln_to_rest_posture(double *x);
void init_soln_to_last(double *x, SL_DJstate target[], double *hitTime);
void finalize_soln(double *x, SL_DJstate target[], double *hitTime);

/* Utility functions */
double inner_prod(const double *a1, const double *a2);
void vec_minus(double *a1, const double *a2);
void print_optim_vec(double *x);
void first_order_hold(double *ballPred, double *racketVel,
		double *racketNormal, double* T);

/* Debug functions */
int test_optim(double *x, double *params, int verbose);

long get_time();

#endif /* BARRETT_SRC_OPTIM_H_ */
