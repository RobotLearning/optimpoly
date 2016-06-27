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
#define EQ_CONSTR_DIM 3*CART
#define INEQ_CONSTR_DIM 2*DOF + 2*DOF // both strike and returning trajectories, min and max
#define MAX_VEL 200
#define MAX_ACC 200
#define TIME2RETURN 1.0

Matrix ballMat; // predicted ball pos and vel values for T_pred time seconds
Matrix racketMat; // racket strategy

// optimization related methods
void nlopt_optim_poly_run(double *x, double *params);
double costfunc(unsigned n, const double *x, double *grad, void *my_func_data);
void kinematics_eq_constr(unsigned m, double *result, unsigned n, const double *x, double *grad, void *f_data);
void joint_limits_ineq_constr(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data);
void calc_strike_poly_coeff(double *a1, double *a2, const double *q0, const double *x);
void calc_return_poly_coeff(double *a1, double *a2, const double *q0, const double *x);
void calc_strike_extrema_cand(double *a1, double *a2, double T, double *q0, double *joint_max_cand, double *joint_min_cand);
void calc_return_extrema_cand(double *a1, double *a2, const double *x, double *joint_max_cand, double *joint_min_cand);
void init_soln_to_rest_posture(double *x, double *q0);
void init_joint_state(double *q0);
void init_ball_state(double *b0, double *v0);
void set_bounds(double *lb, double *ub);

// loading joint limits from SL files
void load_joint_limits();

// debugging methods
void test_optim(double *x, double *params);
void load_lookup_table();

#endif /* OPTIMPOLY_H_ */
