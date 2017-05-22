/**
 * @file optim.h
 *
 * @brief Interface and declarations for the 3 optimization approaches.
 *
 * Exposes VHP, LP and FP optimizations to outside (e.g. Player class
 * can call them).
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
#include "constants.h"

// defines
const int OPTIM_DIM = 2*NDOF+1;
const int EQ_CONSTR_DIM = 3*NCART;
const int INEQ_CONSTR_DIM = 2*NDOF + 2*NDOF; // both strike and returning trajectories, min and max
const double MAX_VEL = 200;
const double MAX_ACC = 200;

/**
 * @brief Desired racket positions, vels and normals for dt*Nmax seconds.
 */
typedef struct {
	double** pos;
	double** vel;
	double** normal;
	double dt;
	int Nmax; // max column length
} racketdes;

/**
 * @brief Optimization coparameters passed to optimizer. (only LP changes them).
 */
typedef struct {
	double* q0;
	double* q0dot;
	double* qrest;
	double* lb;
	double* ub;
	double time2return;
	bool detach;
	bool moving;
	bool verbose;
} coptim; // optimization coparameters

/**
 * @brief Optimization parameters passed to optimizer.
 *
 * When running is TRUE, player does not update trajectories.
 * IF update is TRUE and running is FALSE, then player can update trajectories.
 */
typedef struct {
	double* qf;
	double* qfdot;
	double T;
	bool running;
	bool update;
} optim; // optimization variables

/**
 * @brief Lazy Player uses a different structure.
 */
typedef struct {
	racketdes *racketdata;
	coptim * coparams;
	double **ballpred;
	double dt;
	int Nmax;
} lazy_data;

/**
 * @brief Weights for optimization penalties used by LP.
 */
typedef struct {
	double * R_strike;
	double * R_hit;
	double * R_land;
	double R_net;
} weights; // weights for optimization penalties

class Optim {

private:
	virtual void predict() = 0;
	virtual void optim() = 0;

public:
	virtual void run() = 0;

};

class VHP : public Optim {

private:
	virtual void predict();
	virtual void optim();

public:
	virtual void run();
};

// interface for VHP player
double nlopt_vhp_run(coptim *coparams,
					 racketdes *racketdata,
					 optim *params);

// interface for LAZY player
double nlopt_optim_lazy_run(double** ballpred,
		              coptim *coparams,
	                  racketdes *racketdata,
		              optim *params);

// interface for FIXED player
double nlopt_optim_fixed_run(coptim * coparams,
		                     racketdes * racketdata,
							 optim * params);

#endif /* OPTIMPOLY_H_ */
