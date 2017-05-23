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
#include <armadillo>

using namespace arma;

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

struct racket_des {
	double pos[NCART] = {0.0};
	double vel[NCART] = {0.0};
	double normal[NCART] = {0.0};
}; // racket pos, vel and normals

class Optim {

protected:
	bool verbose = false;
	bool moving = false;
	bool update = false;
	bool running = false;
	bool detach = false;
	nlopt_opt opt;
	double qrest[NDOF] = {0.0};
	double qf[NDOF] = {0.0};
	double qfdot[NDOF] = {0.0};
	double T = 1.0;
	virtual void init_last_soln() const = 0;
	virtual void init_rest_soln() const = 0;
	virtual double test_soln() const = 0;
	virtual void finalize_soln() = 0;
	virtual bool predict() = 0;
	virtual void optim() = 0;
	void calc_hit_traj() {
		if (predict()) {
			optim();
		}
	}
public:
	bool get_params();
	void run() {
		// run optimization in another thread
		std::thread t(calc_hit_traj);
		if (detach)
			t.detach();
		else
			t.join();
	}

};

class VHP : public Optim {

protected:
	virtual bool predict();
	virtual void optim();
	virtual void init_last_soln(double x[2*NDOF]) const;
	virtual void init_rest_soln(double x[2*NDOF]) const;
	virtual double test_soln(const double x[2*NDOF]) const;
	virtual void finalize_soln(const double x[2*NDOF], double time_elapsed);
private:
	double limit_avg[NDOF];
	racket_des racket;
	vec2 ball_land_des;
	double time_land_des;
	EKF & filter;
	game game_state;
public:
	VHP(double qrest_[NDOF], double lb[NDOF], double ub[NDOF]);
};

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
