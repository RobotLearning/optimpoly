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
 * @brief Desired racket/ball positions, vels and racket normals for dt*Nmax seconds.
 */
typedef struct {
	double** pos;
	double** vel;
	double** normal;
	double dt;
	int Nmax; // max column length
} optim_des;

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
	optim_des *racketdata;
	coptim * coparams;
	double **ballpred;
	double dt;
	int Nmax;
} lazy_data;

/**
 * @brief Weights for optimization penalties used by LP.
 */
typedef struct {
	double R_strike[NDOF] = {0.0};
	double R_hit[NCART] = {0.0};
	double R_land[2] = {0.0};
	double R_net = 0.0;
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
	optim_des *param_des;
	double lb[NDOF];
	double ub[NDOF];
	double qrest[NDOF] = {0.0};
	double q0[NDOF] = {0.0};
	double q0dot[NDOF] = {0.0};
	double time2return = 1.0;
	virtual ~Optim() {};
	bool get_params(double qf_[NDOF], double qfdot_[NDOF], double T_);
	void update_init_state(double *j0, double *j0dot, double time_pred);
	void set_des_params(optim_des *params);
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
	weights w;
	LazyOptim(double qrest[NDOF], double lb[NDOF], double ub[NDOF]);
};

// functions that all players use
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

#endif /* OPTIMPOLY_H_ */
