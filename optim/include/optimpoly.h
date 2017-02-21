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
} racket;

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
	int update;
} optim; // optimization variables

// interface
double nlopt_optim_poly_run(coptim * coparams, racket * racket, optim * params);

#endif /* OPTIMPOLY_H_ */
