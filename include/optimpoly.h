/*
 * optimpoly.h
 *
 *  Created on: Jun 7, 2016
 *      Author: okoc
 */

#ifndef OPTIMPOLY_H_
#define OPTIMPOLY_H_

// optimization and math libraries
#include <math.h>
#include <nlopt.h>

// SL variables and kinematics
#include "SL.h"
#include "SL_user.h"
#include "SL_user_common.h"

// defines
#define DOF 7
#define OPTIM_DIM 2*DOF+1
#define CONSTR_DIM 3*CART
#define MAX_VEL 200
#define MAX_ACC 200

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

// initialization needs to be done for this mapping - used in SL_kinematics_body.h
int  link2endeffmap[] = {0,PALM};

static double q0[DOF]; // starting position of the robot

SL_Cstate ballPred; // ball status variable used for prediction
Matrix ballMat; // predicted ball pos and vel values for T_pred time seconds
Matrix racketMat; // racket strategy

// optimization related methods
double costfunc(unsigned n, const double *x, double *grad, void *my_func_data);
void kinematics_constr(unsigned m, double *result, unsigned n, const double *x, double *grad, void *f_data);
void calc_poly_coeff(double *a1, double *a2, const double *q0, const double *x);
void optim_poly_nlopt_run();
void init_soln(double * x);
void init_joint_state();
void init_ball_state();
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

#endif /* OPTIMPOLY_H_ */
