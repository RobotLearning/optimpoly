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

// defines
#define NDOF 7
#define OPTIM_DIM 2*NDOF+1
#define EQ_CONSTR_DIM 3*NCART
#define INEQ_CONSTR_DIM 2*NDOF + 2*NDOF // both strike and returning trajectories, min and max
#define MAX_VEL 200
#define MAX_ACC 200

typedef struct {
	double** pos;
	double** vel;
	double** normal;
	int Nmax; // max column length
} cracket;

typedef struct {
	double* q0;
	double* q0dot;
	double* qrest;
	double* lb;
	double* ub;
	double time2return;
} coptim; // optimization co-parameters

typedef struct {
	double* qf;
	double* qfdot;
	double T;
} optim; // optimization variables

// interface
double nlopt_optim_poly_run(coptim * coparams, cracket * racket, optim * params);

// termination
static double test_optim(double *x, coptim *params, cracket *racket, int info);
static void finalize_soln(const double* x, optim * params);

// optimization related methods
static double costfunc(unsigned n, const double *x, double *grad, void *my_func_data);
static double const_costfunc(unsigned n, const double *x, double *grad, void *my_func_params) ;
static void kinematics_eq_constr(unsigned m, double *result, unsigned n,
		                  const double *x, double *grad, void *f_data);
static void joint_limits_ineq_constr(unsigned m, double *result,
		                      unsigned n, const double *x, double *grad, void *data);

static void calc_strike_poly_coeff(const double *q0, const double *q0dot, const double *x,
		                    double *a1, double *a2);
static void calc_return_poly_coeff(const double *q0, const double *q0dot,
		                    const double *x, const double time2return,
		                    double *a1, double *a2);
static void calc_strike_extrema_cand(const double *a1, const double *a2, const double T,
		                      const double *q0, const double *q0dot,
							  double *joint_max_cand, double *joint_min_cand);
static void calc_return_extrema_cand(const double *a1, const double *a2,
		                      const double *x, const double time2return,
							  double *joint_max_cand, double *joint_min_cand);
static void init_soln(const optim * params, double x[OPTIM_DIM])

static void first_order_hold(const cracket* racket, const double T, double racket_pos[NCART],
		               double racket_vel[NCART], double racket_n[NCART]);
#endif /* OPTIMPOLY_H_ */
