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
	double dt;
	int Nmax; // max column length
} racketdes;

typedef struct {
	double* q0;
	double* q0dot;
	double* qrest;
	double* lb;
	double* ub;
	double time2return;
} coptim; // optimization coparameters

typedef struct {
	double* qf;
	double* qfdot;
	double T;
	int running;
	int update;
} optim; // optimization variables

typedef struct {
	racketdes *racketdata;
	coptim * coparams;
	double **ballpred;
	double *qwait;
	double dt;
	int Nmax;
} lazy_data;

typedef struct {
	double * R_strike;
	double * R_return;
	double * R_wait;
	double * R_hit;
	double * R_land;
	double R_net;
} weights; // weights for optimization penalties

// interface for LAZY player
double nlopt_optim_lazy_run(double** ballpred,
		              coptim *coparams,
	                  racketdes *racketdata,
		              optim *params);

// interface for FIXED player
double nlopt_optim_fixed_run(coptim * coparams,
		                     racketdes * racketdata,
							 optim * params);

// interface for VHP player
double nlopt_vhp_run(coptim *coparams,
					 racketdes *racketdata,
					 optim *params);

#endif /* OPTIMPOLY_H_ */
