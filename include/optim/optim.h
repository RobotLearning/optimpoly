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
#include <thread>
#include <math.h>
#include <nlopt.h>
#include "string.h" //for bzero
#include "constants.h"

// defines
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

protected:
	static const int OPTIM_DIM = 2*NDOF;
	bool verbose = true;
	bool moving = false;
	bool update = false;
	bool running = false;
	bool detach = false;
	nlopt_opt opt;

	double qf[NDOF] = {0.0};
	double qfdot[NDOF] = {0.0};
	double T = 1.0;

	virtual void init_last_soln(double *x) const = 0;
	virtual void init_rest_soln(double *x) const = 0;
	virtual double test_soln(const double *x) const = 0;
	virtual void finalize_soln(const double *x, const double dt) = 0;
	virtual void optim();
public:
	racketdes *racket;
	double lb[NDOF];
	double ub[NDOF];
	double qrest[NDOF] = {0.0};
	double q0[NDOF] = {0.0};
	double q0dot[NDOF] = {0.0};
	double time2return = 1.0;
	virtual ~Optim() {};
	bool get_params(double qf_[NDOF], double qfdot_[NDOF], double T_);
	void update_init_state(double *j0, double *j0dot, double time_pred);
	void set_des_racket(racketdes *racket);
	void run();
};

class HittingPlane : public Optim {

protected:
	static const int OPTIM_DIM = 2*NDOF;
	//virtual bool predict(EKF & filter);
	virtual void init_last_soln(double x[]) const;
	virtual void init_rest_soln(double x[]) const;
	virtual double test_soln(const double x[]) const;
	virtual void finalize_soln(const double x[], const double time_elapsed);
public:
	double limit_avg[NDOF];
	HittingPlane(double qrest[NDOF], double lb[NDOF], double ub[NDOF]);
};

class FocusedOptim : public Optim {

protected:
	static const int OPTIM_DIM = 2*NDOF + 1;
	//virtual bool predict(EKF & filter);
	virtual void init_last_soln(double x[]) const;
	virtual void init_rest_soln(double x[]) const;
	virtual double test_soln(const double x[]) const;
	virtual void finalize_soln(const double x[], const double time_elapsed);
public:
	FocusedOptim() {}; // for lazy player
	FocusedOptim(double qrest[NDOF], double lb[NDOF], double ub[NDOF]);
};

class LazyOptim : public FocusedOptim {

protected:
	virtual double test_soln(const double x[]) const;
public:
	double **ballpred;
	weights *w;
	void set_ball_pred(double **ballpred);
	LazyOptim(double qrest[NDOF], double lb[NDOF], double ub[NDOF]);
};

// interface for LAZY player
double nlopt_optim_lazy_run(double** ballpred,
		              coptim *coparams,
	                  racketdes *racketdata,
		              optim *params);

#endif /* OPTIMPOLY_H_ */
