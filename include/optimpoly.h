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
#include <pthread.h>
#include "string.h" //for bzero

// SL variables and kinematics
#include "SL.h"

// defines
#define DOF 7
#define OPTIM_DIM 2*DOF+1
#define EQ_CONSTR_DIM 3*CART
#define INEQ_CONSTR_DIM 2*DOF + 2*DOF // both strike and returning trajectories, min and max
#define MAX_VEL 200
#define MAX_ACC 200

extern Matrix ballMat; // predicted ball pos and vel values for T_pred time seconds
extern Matrix racketMat; // racket strategy

// optimization related methods
void nlopt_optim_poly_run(double *x);
double costfunc(unsigned n, const double *x, double *grad, void *my_func_data);
double const_costfunc(unsigned n, const double *x, double *grad, void *my_func_params) ;
void kinematics_eq_constr(unsigned m, double *result, unsigned n,
		                  const double *x, double *grad, void *f_data);
void joint_limits_ineq_constr(unsigned m, double *result,
		                      unsigned n, const double *x, double *grad, void *data);

void calc_strike_poly_coeff(const double *q0, const double *q0dot, const double *x,
		                    double *a1, double *a2);
void calc_return_poly_coeff(const double *q0, const double *q0dot,
		                    const double *x, const double time2return,
		                    double *a1, double *a2);
void calc_strike_extrema_cand(const double *a1, const double *a2, const double T,
		                      const double *q0, const double *q0dot,
							  double *joint_max_cand, double *joint_min_cand);
void calc_return_extrema_cand(const double *a1, const double *a2,
		                      const double *x, const double time2return,
							  double *joint_max_cand, double *joint_min_cand);

void init_soln_to_rest_posture(double *x, double *q0);
void init_joint_state(double *q0);
void init_ball_state(double *b0, double *v0);
void set_bounds(double *lb, double *ub, double SLACK);

// debugging methods
void test_optim(double *x);
void load_lookup_table(Matrix lookupTable);

#endif /* OPTIMPOLY_H_ */
