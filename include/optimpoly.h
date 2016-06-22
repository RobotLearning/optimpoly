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
#include "string.h" //for bzero

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

Matrix ballMat; // predicted ball pos and vel values for T_pred time seconds
Matrix racketMat; // racket strategy

// optimization related methods
double costfunc(unsigned n, const double *x, double *grad, void *my_func_data);
void kinematics_constr(unsigned m, double *result, unsigned n, const double *x, double *grad, void *f_data);
void calc_poly_coeff(double *a1, double *a2, const double *q0, const double *x);
void optim_poly_nlopt_run(double *params);
void init_soln(double *x, double *q0);
void init_joint_state(double *q0);
void init_ball_state(double *b0, double *v0);
void set_bounds(double *lb, double *ub);

// loading joint limits from SL files
void load_joint_limits();

// debugging methods
void test_constraint(double *x, double *params);
void lookup(double *x);

#endif /* OPTIMPOLY_H_ */
