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
const double MAX_VEL = 10;
const double MAX_ACC = 200;

using namespace arma;

/**
 * @brief Desired racket/ball positions, vels and racket normals for dt*Nmax seconds.
 *
 * This is the structure passed to Optimization class Optim
 * and its descendants.
 * For FP, we compute desired racket positions, velocities and normals
 * based on predicted ball path inside robot workspace.
 * For DP, we only use the ball predicted positions and velocities.
 */
struct optim_des {
	mat racket_pos = zeros<mat>(NCART,1); //!< racket desired pos
	mat racket_vel = zeros<mat>(NCART,1); //!< racket desired vel
	mat racket_normal = zeros<mat>(NCART,1); //!< racket desired normals
	mat ball_pos = zeros<mat>(NCART,1); //!< incoming ball predicted pos.
	mat ball_vel = zeros<mat>(NCART,1); //!< incoming ball predicted vels.
	double dt = DT; //!< time step between each prediction
	int Nmax = 1; //!< max number of time steps for prediction
};

/**
 * @brief 3rd order spline parameters returned from optimizer.
 *
 * If optimization is still running, player does not update trajectories.
 * IF optimization was successful (update is TRUE) and
 * it has terminated (running is FALSE), then player class can update trajectories.
 */
struct spline_params {
	double time2hit = 1.0; //!< free-final time (for hitting ball)
	mat a = zeros<mat>(NDOF,4); //!< strike poly params of 3rd order
	mat b = zeros<mat>(NDOF,4); //!< return poly params of 3rd order
};

/**
 * @brief Desired/actual joint positions, velocities, accelerations.
 *
 * Main Player function play() updates desired joint values
 * every 2ms for Barrett WAM.
 */
struct joint {
	vec7 q = zeros<vec>(NDOF); //!< desired joint pos
	vec7 qd = zeros<vec>(NDOF); //!< desired joint vels
	vec7 qdd = zeros<vec>(NDOF); //!< desired joint accs
};

/**
 * @brief Weights for optimization penalties used by DP only.
 *
 * These weights are NOT used in any interesting way so far.
 * It would be interesting to learn/update these weights based on
 * successful returns in real robot experiments.
 */
struct weights {
	double R_strike[NDOF] = {0.0}; //!< acceleration weights for running cost
	double R_hit = 0.0; //!< weight of dist from racket centre to hit location
	double R_land = 0.0; //!< weight of dist from landing pos to centre
	double R_net = 0.0; //!< weight of dist from net pos to a suitable pos above net
};

/**
 * @brief Data needed for resting state optimization
 *
 * Optimization tries to find a good resting posture close to the predicted ball locations ball_pred
 * and a given hitting state q_hit with low jacobian norm (Frobenius)
 */
struct rest_optim_data {
	mat ball_pred;
	vec7 q_hit;
};

/**
 * @brief Base class for all optimizers.
 *
 * Containts base class methods, members and virtual methods
 * to be implemented.
 */
class Optim {

protected:
	static const int OPTIM_DIM = 2*NDOF + 1; //!< dim. of optim problem
	bool lookup = false; //!< use lookup table methods to init. optim params.
	bool verbose = true; //!< verbose output printed
	bool moving = false; //!< robot is already moving so use last computed values to init.
	bool update = false; //!< optim finished and soln. seems valid
	bool running = false; //!< optim still RUNNING
	bool detach = false; //!< detach optim in another thread
	mat lookup_table; //!< lookup table used to init. optim values
	nlopt_opt opt; //!< optimizer from NLOPT library

	double qf[NDOF] = {0.0}; //!< saved joint positions after optim
	double qfdot[NDOF] = {0.0}; //!< saved joint velocities after optim
	double T = 1.0; //!< saved hitting time after optim terminates

	void init_lookup_soln(double *x);
	virtual void init_last_soln(double *x) const = 0;
	virtual void init_rest_soln(double *x) const = 0;
	virtual double test_soln(const double *x) const = 0;
	virtual void finalize_soln(const double *x, const double dt) = 0;
	void optim_rest_posture(vec7 & q_rest_des);
	void update_rest_state(const vec7 & q_rest_new);
	void optim();
public:
	optim_des *param_des; //!< Desired racket and/or ball predicted vals.
	double lb[2*NDOF+1]; //!< Joint lower limits, joint vel. lower limit and min. hitting time
	double ub[2*NDOF+1]; //!< Joint upper limits, joint vel. upper limit and max. hitting time
	double qrest[NDOF] = {0.0}; //!< FIXED Resting posture for optimizers to compute return traj.
	double q0[NDOF] = {0.0}; //!< Initial joint state needed to compute traj. acc.
	double q0dot[NDOF] = {0.0}; //!< Initial joint velocities needed to compute traj. acc.
	double time2return = 1.0; //!< FIXED Time to return
	virtual ~Optim() {};
	bool check_update();
	bool check_running();
	void set_moving(bool flag);
	void set_detach(bool flag);
	void set_verbose(bool flag);
	bool get_params(const joint & qact, spline_params & p);
	void run_qrest_optim(vec7 & q_rest_des);
	void update_init_state(const joint & qact);
	void set_des_params(optim_des *params);
	void run();
};

/**
 * @brief Optimizer for Virtual Hitting Plane player.
 *
 * Fixes hitting plane and ALSO fixes desired landing position and time.
 */
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
	void fix_hitting_time(double time_pred);
	HittingPlane(const vec7 & qrest_, double lb_[], double ub_[]);
};

/**
 * @brief Optimizer for Focused Player
 *
 * Focused Player fixes a desired landing position and time for the ball.
 */
class FocusedOptim : public Optim {

protected:
	static const int OPTIM_DIM = 2*NDOF + 1; //!< dim. of optim
	//virtual bool predict(EKF & filter);
	virtual void init_last_soln(double x[]) const; //!< see derived class
	virtual void init_rest_soln(double x[]) const; //!< see derived class
	virtual double test_soln(const double x[]) const; //!< see derived class
	virtual void finalize_soln(const double x[], const double time_elapsed); //!< see derived class
public:
	FocusedOptim() {}; // for lazy player
	FocusedOptim(const vec7 & qrest_, double lb_[], double ub_[]);
};

/**
 * @brief Optimizer for Defensive Player
 *
 * Defensive Player doesn't fix desired ball pos. or time!
 */
class DefensiveOptim : public FocusedOptim {

private:
	virtual double test_soln(const double x[]) const;
	void set_land_constr();
	void set_hit_constr();
	virtual void finalize_soln(const double x[], const double time_elapsed);
	std::vector<double> mult_vel = {0.9,0.8,0.83};
public:
	bool land; //!< compute strikes for landing if true or hitting only if false
	weights w; //!< weights of optimization
	double x_last[OPTIM_DIM] = {0.0}; //!< last iteration values
	double t_land = -1.0; //!< computed landing time
	double t_net = -1.0; //!< computed net passing time
	double x_land[NCART] = {0.0}; //!< computed ball landing pos.
	double x_net[NCART] = {0.0}; //!< computed ball net pass. pos.
	double dist_b2r_norm = 1.0; //!< normal dist. from ball to racket
	double dist_b2r_proj = 1.0; //!< dist. from ball to racket proj. to racket plane
	std::vector<double> penalty_loc = {0.0, 0.23, 0.0, -3.22}; //!< penalty locations for landing and net
	void set_velocity_multipliers(const std::vector<double> & mult);
	void set_weights(const std::vector<double> & weights);
	void set_penalty_loc(const std::vector<double> & penalty_loc_);
	double calc_punishment();
	void calc_times(const double x[]);
	void calc_hit_distance(const double bp[], const double rp[], const double n[]);
	DefensiveOptim(const vec7 & qrest_, double lb_[], double ub_[], bool land_ = true, bool lookup_ = false);
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
double calc_max_acc_violation(const double x[2*NDOF+1],
								const double q0[NDOF],
								const double q0dot[NDOF]);
#endif /* OPTIMPOLY_H_ */
