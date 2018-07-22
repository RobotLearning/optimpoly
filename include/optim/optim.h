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
#include "kalman.h" // for estimate_prior

using arma::mat;
using arma::vec;
using arma::zeros;
using arma::vec6;
using arma::vec7;
using const_tt::NCART;
using const_tt::NDOF;

// defines
const int EQ_CONSTR_DIM = 3*NCART;
const int INEQ_CONSTR_DIM = 2*NDOF + 2*NDOF; // both strike and returning trajectories, min and max
const double MAX_VEL = 10;
const double MAX_ACC = 200;

namespace optim {

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
	double dt = const_tt::DT; //!< time step between each prediction
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
	mat lookup_table = zeros<mat>(1,OPTIM_DIM+2*NCART); //!< lookup table used to init. optim values
	nlopt_opt opt = nullptr; //!< optimizer from NLOPT library

	double qf[NDOF] = {0.0}; //!< saved joint positions after optim
	double qfdot[NDOF] = {0.0}; //!< saved joint velocities after optim
	double T = 1.0; //!< saved hitting time after optim terminates

	/** @brief Initialize optimization parameters using a lookup table.
	 *
	 * Call this function AFTER setting desired BALL parameters.
	 * TODO: include k as a parameter
	 * @param x Array of robot parameters qf,qfdot,T to be updated
	 */
	void init_lookup_soln(double *x);

	virtual void init_last_soln(double *x) const = 0;
	virtual void init_rest_soln(double *x) const = 0;
	virtual double test_soln(const double *x) const = 0;
	virtual void finalize_soln(const double *x, const double dt) = 0;

	/**
	 * @brief Resting posture optimization that tries to find good resting joint positions
	 *
	 * The resting joint positions are optimized to be
	 * close to the predicted incoming ball trajectory and hitting joint states with low
	 * Frobenius-norm jacobian
	 */
	void optim_rest_posture(vec7 & q_rest_des);

	/**
	 * @brief Update the rest state from outside
	 *
	 * If there is an additional optimization somewhere else
	 * that optimizes for the resting posture, notify the optim classes
	 * @param q_rest_new
	 */
	void update_rest_state(const vec7 & q_rest_new);

	/** @brief NLOPT optimization happens here. */
	void optim();

public:
	optim_des *param_des = nullptr; //!< Desired racket and/or ball predicted vals.
	double lb[2*NDOF+1]; //!< Joint lower limits, joint vel. lower limit and min. hitting time
	double ub[2*NDOF+1]; //!< Joint upper limits, joint vel. upper limit and max. hitting time
	double qrest[NDOF] = {0.0}; //!< FIXED Resting posture for optimizers to compute return traj.
	double q0[NDOF] = {0.0}; //!< Initial joint state needed to compute traj. acc.
	double q0dot[NDOF] = {0.0}; //!< Initial joint velocities needed to compute traj. acc.
	double time2return = 1.0; //!< FIXED Time to return

	/** @brief Destroy nlopt structure */
	virtual ~Optim();

	/**
	 * @brief If the optimization was successful notify the Player class
	 *
	 * If the optimization was successful, update is turned ON and the table tennis
	 * player can launch/update the polynomial trajectories.
	 * @return update
	 */
	bool check_update();

	/**
	 * @brief Tells the player optimization thread is still BUSY.
	 *
	 * If the (detached) thread is still running then table tennis player does not
	 * update/launch new trajectories.
	 * @return running
	 */
	bool check_running();

	/**
	 * @brief If the robot starts moving the optimization is notified via this function
	 *
	 * If the robot is moving, this means last optimization was feasible, hence
	 * we can initialize the new optimization from the last optimized parameter values.
	 * @param flag_move
	 */
	void set_moving(bool flag);

	/**
	 * @brief Detach the optimization thread.
	 *
	 * If the optimization is performed on SL or REAL_ROBOT then the optimization
	 * thread should be detached.
	 * @param flag_detach
	 */
	void set_detach(bool flag);

	/** @brief Print verbose optimization output (detailed optimization results are printed) */
	void set_verbose(bool flag);

	/**
	 * @brief If optimization succeeded, update polynomial parameters p
	 *
	 * If the optimizers finished running and terminated successfully,
	 * then generates the striking and returning polynomial parameters
	 * from qf, qfdot and T, and given the actual joint states qact, qactdot
	 *
	 */
	bool get_params(const joint & qact, spline_params & p);

	/**
	 * @brief Run q_rest optimization to find a suitable resting posture
	 *
	 * Detaches the optimization if detach is set to true
	 */
	void run_qrest_optim(vec7 & q_rest_des);

	/**
	 * @brief Update the initial state of optimization to PLAYER's current joint states.
	 * @param qact Initial joint states acquired from sensors
	 */
	void update_init_state(const joint & qact);

	/** @brief Set the final time for the returning trajectory */
	void set_return_time(const double & time);

	/**
	 * @brief Set desired optimization parameters before running optim.
	 *
	 * @param params_ Desired optimization parameters are racket and/or ball values
	 * predicted or computed by player class.
	 */
	void set_des_params(optim_des *params);

	/**
	 * @brief Runs the optimization.
	 *
	 * Detaches the optimization if detach is set to TRUE. The
	 * optimization method is shared by all Optim class descendants (VHP,FP,DP).
	 */
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

	/**
	 * @brief Initialize solution using already computed soln. from prev. optim
	 *
	 * If robot is already moving, use the last optimized solution
	 * to initialize current optim
	 * @param x
	 */
	virtual void init_last_soln(double x[]) const;

	/**
	 * @brief Initialize solution using a rest posture
	 *
	 * Resting posture and zero-velocities and 0.5 hitting time
	 * is used to initialize optim params.
	 * @param x
	 */
	virtual void init_rest_soln(double x[]) const;

	/**
	 * @brief Check kinematics constraints
	 *
	 * Before updating trajectory, we need to make sure that
	 * the hitting constraints and joint limits are not violated.
	 * @param x
	 * @return Maximum violation
	 */
	virtual double test_soln(const double x[]) const;

	/**
	 * If predicted hitting time using the virtual hitting plane
	 * is more than 50 ms, then update!
	 * @param x
	 * @param time_elapsed
	 */
	virtual void finalize_soln(const double x[], const double time_elapsed);

public:
	double limit_avg[NDOF];

	/**
	 * @brief FIX desired hitting time for VHP trajectories!
	 * This is actually only used by Virtual Hitting Plane!
	 * @param time_pred If this is more than 50 ms, set desired
	 * hitting time at virtual hitting plane to this value.
	 */
	void fix_hitting_time(double time_pred);

	/**
	 * @brief Initialize VHP optimizer
	 * @param qrest_ FIXED rest posture
	 * @param lb_ Joint lower limits (pos,vel limits and time limits incl.)
	 * @param ub_ Joint upper limits (pos,vel limits and time limits incl.)
	 */
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

	/** @brief Initialize the optimization parameters the last optimized solution values.*/
	virtual void init_last_soln(double x[]) const;

	/**
	 * @brief Initialize the optim params to fixed resting posture.
	 *
	 * Initializes the optim params to qf fixed to q_rest, zero velocites,
	 * and 0.5 hitting time.
	 * @param x Optim params
	 */
	virtual void init_rest_soln(double x[]) const;

	/**
	 * @brief Test solution with hard kinematics constraints
	 *
	 * If constraints are violated then do not update/init. trajectories!
	 *
	 * @param x Optim params
	 * @return Maximum value of constraint violations.
	 */
	virtual double test_soln(const double x[]) const;

	/**
	 * @brief Finalize solution if more than 50 ms is available for hitting.
	 * @param x Optim params
	 * @param time_elapsed Time elapsed during optimization
	 */
	virtual void finalize_soln(const double x[], const double time_elapsed);

public:

	/** @brief Constructor useful for lazy player (subclass) */
	FocusedOptim() {};

	/**
	 * @brief Initialize the NLOPT optimization procedure here for FP
	 * @param qrest_ Fixed resting posture
	 * @param lb_ Fixed joint lower limits
	 * @param ub_ Fixed joint upper limits
	 */
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

	/**
	 * @brief Setting LANDING constraints, i.e., not just hitting the ball.
	 *
	 * Here we consider the full landing constraints where the ball has
	 * to pass over the net and land on the opponent's court.
	 * We can study the simpler optimization problem of hitting the ball
	 * by switching to pure hitting constraints.
	 */
	void set_land_constr();

	/**
	 * @brief Switch to only HITTING constraints, not landing.
	 *
	 * We can study the simpler optimization problem of hitting the ball
	 * by switching to pure hitting constraints.
	 */
	void set_hit_constr();

	/**
	 * @brief Finalize solution if more than 50 ms is available for hitting.
	 * @param x Optim params
	 * @param time_elapsed Time elapsed during optimization
	 */
	virtual void finalize_soln(const double x[],
	                            const double time_elapsed);
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

	/**
	 * @brief Ball velocity post-multiplier is set here
	 *
	 * This is useful for parameterizing Defensive Player for REAL ROBOT
	 * to compensate for any errors
	 * @param mult
	 */
	void set_velocity_multipliers(const std::vector<double> & mult);

	/**
	 * @brief Set weights for Defensive Player
	 * @param weights weights for hitting, net and landing penalties in cost function
	 */
	void set_weights(const std::vector<double> & weights);
	void set_penalty_loc(const std::vector<double> & penalty_loc_);

	/**
	 * @brief Calculate punishment for hitting, net and landing.
	 *
	 * Calculate the costs for hitting the ball, and then if land variable is set to true
	 * also calculate the costs for passing it over to a desired net location (x,z loc)
	 * and landing to a desired landing point (x,y loc.)
	 * @return hitting and landing penalties (if land is turned on combined, else only hitting)
	 */
	double calc_punishment();

	/**
	 * @brief Calculate the net hitting time and the landing time
	 * Assuming that there was an interaction (hitting) between racket and ball.
	 *
	 * Since we call this function multiple times within one optimization step
	 * we store the optimization variable x,
	 * landTime and netTime variables and return these when the optimization
	 * variable is the same (meaning optimization algorithm did not update it yet).
	 *
	 */
	void calc_times(const double x[]);


	/**
	 * @brief Calculate deviation and projection norms of ball to racket.
	 *
	 * Calculates normal and projected distance from ball to racket
	 * For satisfying hitting constraints.
	 *
	 */
	void calc_hit_distance(const double bp[],
	                        const double rp[],
	                        const double n[]);

	/**
	 * @brief Initialize Defensive Player
	 * @param qrest_ FIXED resting posture
	 * @param lb_ Lower joint pos,vel limits and min. hitting time
	 * @param ub_ Upper joint pos,vel limits and max. hitting time
	 * @param land_ Optimize for returning the ball (true) or only hitting (false)
	 * @param lookup_ Lookup optim params from table if true
	 */
	DefensiveOptim(const vec7 & qrest_,
	                double lb_[],
	                double ub_[],
	                bool land_ = true,
	                bool lookup_ = false);
};

/**
 * @brief Inequality constraint that makes sure we never exceed the
 * joint limits during the striking and returning motion
 *
 */
void joint_limits_ineq_constr(unsigned m,
                              double *result,
		                      unsigned n,
		                      const double *x,
		                      double *grad,
		                      void *data);
/**
 * @brief Calculate the polynomial coefficients from the optimized variables qf,qfdot,T
 * p(t) = a1*t^3 + a2*t^2 + a3*t + a4
 */
void calc_strike_poly_coeff(const double *q0,
                            const double *q0dot,
                            const double *x,
		                    double *a1,
		                    double *a2);

/**
 * @brief Calculate the returning polynomial coefficients from the optimized variables qf,qfdot
 * and time to return constant T
 * p(t) = a1*t^3 + a2*t^2 + a3*t + a4
 */
void calc_return_poly_coeff(const double *q0,
                            const double *q0dot,
		                    const double *x,
		                    const double time2return,
		                    double *a1,
		                    double *a2);

/**
 * @brief Calculate the extrema candidates for each joint (2*7 candidates in total)
 * For the striking polynomial
 * Clamp to [0,T]
 *
 */
void calc_strike_extrema_cand(const double *a1,
                              const double *a2,
                              const double T,
		                      const double *q0,
		                      const double *q0dot,
							  double *joint_max_cand,
							  double *joint_min_cand);

/**
 * @brief Calculate the extrema candidates for each joint (2*7 candidates in total)
 * For the return polynomial
 * Clamp to [0,TIME2RETURN]
 *
 */
void calc_return_extrema_cand(const double *a1,
                              const double *a2,
		                      const double *x,
		                      const double time2return,
							  double *joint_max_cand,
							  double *joint_min_cand);

/**
 * @brief Calculate the max acc throughout traj.
 *
 * Since the accelerations of third order polynomials
 * are linear functions of time we check the values
 * at start of traj, t = 0 and end of traj, t = T_hit
 * which are given by 6*a1*T + a2 and a2 respectively
 *
 */
double calc_max_acc_violation(const double x[2*NDOF+1],
						      const double q0[NDOF],
						      const double q0dot[NDOF]);

/**
 * @brief Set upper and lower bounds on the optimization.
 * First loads the joint limits and then puts some slack
 */
void set_bounds(double *lb, double *ub, double SLACK, double Tmax);

/**
 * @brief Compute desired racket pos,vel,normals and/or ball positions, vels.
 * Function that calculates a racket strategy : positions, velocities and racket normal
 * for each point on the predicted ball trajectory (ballMat)
 * to return it a desired point (ballLand) at a desired time (landTime)
 *
 */
optim_des calc_racket_strategy(const mat & balls_predicted,
                               const arma::vec2 & ball_land_des,
                               const double time_land_des,
                               optim_des & racket_params);

/**
 * @brief Compute desired racket pos,vel,normals and/or ball positions, vels. assuming spin model
 * Function that calculates a racket strategy : positions, velocities and racket normal
 * for each point on the predicted ball trajectory (ballMat)
 * to return it a desired point (ballLand) at a desired time (landTime)
 *
 * As opposed to calculating with spin-free models, this function
 * runs an optimization for each predicted ball to find desired outgoing ball velocities!
 *
 */
optim_des calc_spin_racket_strategy(const mat & balls_predicted,
                                    const double & topspin,
                                    const arma::vec3 & ball_land_des,
                                    const double time_land_des,
                                    optim_des & racket_params);

/**
 * @brief Estimates initial ball state + ball topspin
 *
 * @param observations Ball positions
 * @param times Ball time stamps for each position observation
 * @param verbose Verbose output for estimation if true
 * @param NLOPT_FINISHED If detached, the thread will communicate that it has finished
 * @param filter Filter state will be initialized after estimation
 */
void estimate_prior(const mat & observations,
                    const mat & times,
                    const int & verbose,
                    bool & init_ball,
                    player::EKF & filter);

/**
 * @brief Least squares to estimate prior given
 * matrix of observations Y of column length N, each row
 * corresponding to one
 *
 * If observations arrive as column vectors then we take
 * transpose of it.
 *
 */
void estimate_ball_linear(const mat & observations,
                          const vec & times,
                          const bool verbose,
                          vec6 & init_est);
}

#endif /* OPTIMPOLY_H_ */
