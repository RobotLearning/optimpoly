/*
 * optimpoly.h
 *
 *  Created on: Jun 7, 2016
 *      Author: okoc
 */

#ifndef OPTIMPOLY_H_
#define OPTIMPOLY_H_

#include <stdio.h>
#include <stdlib.h>

// optimization and math libraries
#include <math.h>
#include <nlopt.h>

// SL variables and kinematics
#include "SL.h"
#include "SL_user.h"
#include "SL_user_common.h"
#include "table.h"

// defines
#define CART 3
#define DOF 7
#define OPTIM_DIM 2*DOF+1
#define CONSTR_DIM 3*CART
#define MAX_VEL 200
#define MAX_ACC 200
#define TSTEP 0.01
#define TPRED 1.0

// global variables
char joint_names[][20]= {
  {"dummy"},
  {"R_SFE"},
  {"R_SAA"},
  {"R_HR"},
  {"R_EB"},
  {"R_WR"},
  {"R_WFE"},
  {"R_WAA"}
};

// ball status variable used for prediction
SL_Cstate ballPred;

static double q0[DOF];
Matrix ballMat; // predicted ball pos and vel values for T_pred time seconds
Matrix racketMat; // racket strategy

// utility method
long get_time();
void vec_minus(double *a1, const double *a2);
void vec_plus(double *a1, const double *a2);
double inner_prod(const double *a1, const double *a2);
void const_vec(const int n, const double val, double * vec);
void first_order_hold(double *ballPred, double *racketVel, double *racketNormal, double T);
void print_optim_vec(double *x);

// optimization related methods
double costfunc(unsigned n, const double *x, double *grad, void *my_func_data);
void kinematics_constr(unsigned m, double *result, unsigned n, const double *x, double *grad, void *f_data);
void calc_poly_coeff(double *a1, double *a2, const double *q0, const double *x);
void optim_poly_nlopt_run();
void guesstimate_soln(double * x);
void init_joint_state();
void set_bounds(double *lb, double *ub);

// loading joint limits from SL files
void load_joint_limits();
int read_sensor_offsets(char *fname);

// SL kinematics functions copied for convenience
void revoluteGJacColumn(Vector p, Vector pi, Vector zi, Vector c);
void setDefaultEndeffector(void);

// debugging methods
void test_constraint(double *x);
void lookup(double *x);

// ball related methods
void init_ball_state();
void set_land_parameters(double *ballLand, double *landTime);
void predict_ball_state();
// ball related functions taken from table_tennis_common
void integrateBallState(SL_Cstate ballState, SL_Cstate *ballPred, double deltat, int *bounce); //ball prediction
int checkForBallTableContact(SL_Cstate state);

#endif /* OPTIMPOLY_H_ */
