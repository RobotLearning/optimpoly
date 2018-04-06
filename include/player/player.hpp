/**
 * @file player.hpp
 *
 * @brief Main class for playing Table Tennis declared here.
 *
 *  Created on: Feb 9, 2017
 *      Author: okoc
 */

#ifndef PLAYER_HPP_
#define PLAYER_HPP_

#include "kalman.h"
#include "optim.h"

using namespace arma;

/**
 * @brief Optimizer for Player class.
 */
enum algo {
	VHP,  //!< VHP
	FOCUS,//!< FOCUSED PLAYER
	DP, //!< DEFENSIVE PLAYER
};

/**
 * Finite State machine for Table Tennis
 */
enum game { //trial state
	AWAITING,//!< AWAITING
	ILLEGAL, //!< ILLEGAL
	LEGAL,   //!< LEGAL
	HIT,     //!< HIT
};

/**
 * @brief Options passed to Player class (algorithm, saving, corrections, etc.).
 */
struct player_flags {
	bool detach = false; //!< detach optimizations in another thread
	bool check_bounce = false; //!< check bounce to detect legal incoming ball (turn off for real robot exp.)
	bool outlier_detection = false; //!< outlier detection for real robot
	bool mpc = false; //!< turn on/off corrections
	bool reset = true; //!< reinitializing player class
	bool save = false; //!< saving ball/robot data
	bool spin = false; //!< turn on and off spin-based prediction models
	bool optim_rest_posture = false; //!< turn on rest posture optimization
	algo alg = FOCUS; //!< algorithm for trajectory generation
	int verbosity = 0; //!< OFF, LOW, HIGH, ALL
	int freq_mpc = 1; //!< frequency of mpc updates if turned on
	int min_obs = 5; //!< number of observations to initialize filter
	double out_reject_mult = 2.0; //!< multiplier of variance for outlier detection
	double ball_land_des_offset[2] = {0.0}; //!< desired ball landing offsets (w.r.t center of opponent court)
	double time_land_des = 0.8; //!< desired ball landing time
	double optim_offset = 0.0; //!< offset after net for starting optim (if mpc is off)
	double time2return = 1.0; //!< time to return to starting posture after hit
	double var_noise = 0.001; //!< variance of noise process (R)
	double var_model = 0.001; //!< variance of process noise (Q)
	double t_reset_thresh = 0.3; //!< resetting Kalman filter after this many seconds pass without getting valid obs.
	double VHPY = -0.3; //!< location of hitting plane for VHP method
	std::vector<double> weights = {0.0, 0.0, 0.0}; //!< hit,net,land weights for DP (lazy player)
	std::vector<double> mult_vel = {0.9, 0.8, 0.83}; //!< vel. mult. for DP
	std::vector<double> penalty_loc = {0.0, 0.23, 0.0, -3.22}; //!< penalty locations for DP
};

/**
 *
 * @brief Table Tennis Player class for playing Table Tennis.
 *
 * The methods play() or cheat() must be called every DT milliseconds.
 */
class Player {

private:

	// data fields
	bool init_ball_state = false;
	EKF & filter; // filter for the ball estimation
	vec2 ball_land_des = zeros<vec>(2); // desired landing position
	vec7 q_rest_des; // desired resting joint state
	double t_obs = 0.0; // counting time stamps for resetting filter
	double t_poly = 0.0; // time passed on the hitting spline
	bool valid_obs = true; // ball observed is valid (new ball and not an outlier)
	int num_obs = 0; // number of observations received
	game game_state = AWAITING;
	player_flags pflags;
	optim_des pred_params;
	mat observations; // for initializing filter
	mat times; // for initializing filter
	spline_params poly;
	mat lookup_table;
	Optim *opt; // optimizer

	/**
	 * @brief Filter the blob information with a Kalman Filter.
	 * (Extended) KF is used both in simulation mode and for real robot.
	 *
	 * Checking for new ball that is at least 1 mm away from last observation
	 * Checking for also outliers.
	 * Resets if the ball suddenly appears on the opponent's court.
	 *
	 * Ball is valid if ball is a new ball and (in real robot mode)
	 * it is not an outlier!
	 *
	 * Note: we're assuming that time elasped dt = DT = 0.002 seconds every time!
	 *
	 */
	void estimate_ball_state(const vec3 & obs);

	/**
	 * @brief Run optimizer for FOCUSED PLAYER
	 *
	 * Calculate the optimization parameters using an NLOPT nonlinear optimization algorithm
	 * in another thread
	 *
	 * The optimized parameters are: qf, qf_dot, T
	 * assuming T_return and q0 are fixed
	 *
	 */
	void optim_fp_param(const joint & qact);

	/**
	 * @brief Run optimizer for DEFENSIVE PLAYER
	 *
	 * Calculate the optimization parameters using an NLOPT nonlinear optimization algorithm
	 * in another thread
	 *
	 * The optimized parameters are: qf, qf_dot, T
	 */
	void optim_dp_param(const joint & qact);

	/**
	 * @brief Calculate hitting parameters qf, qfdot
	 * on the Virtual Hitting Plane (VHP) by running Inverse Kinematics
	 *
	 * The inverse kinematics routine runs an optimization to minimize
	 * the distance to a rest posture
	 */
	void optim_vhp_param(const joint & qact);

	void calc_opt_params(const joint & qact);

	/**
	 * @brief Check MPC flag and update if possible
	 *
	 * IF MPC IS TURNED OFF
	 * if ball is incoming and robot is not moving consider optimization
	 *
	 * IF MPC IS TURNED ON
	 * then additionally consider (after running initial optimization)
	 * relaunching optimization if ball is valid (new ball and not an outlier)
	 * the frequency of updates is respected, and the ball has not passed the y-limit
	 *
	 */
	bool check_update(const joint & qact) const;

	/**
	 * @brief Unfold the next desired state of the 3rd order polynomials in joint space
	 * If movement finishes then the desired state velocities and accelerations are zeroed.
	 *
	 * Multithreading : if after initial lookup, the computation in the other thread terminates, then we
	 * synchonize the values between the threads and compute the next desired states given these new polynomial
	 * parameters (qf, qf_dot and T_hit)
	 *
	 */
	void calc_next_state(const joint & qact, joint & qdes);

	/** @brief Start moving pre-optim based on lookup if lookup flag is turned ON. */
	void lookup_soln(const vec6 & ball_state, const int k, const joint & qact);

public:

	/**
	 * @brief Initialize Table Tennis Player.
	 *
	 * Table Tennis Player class that can run 3 different trajectory
	 * generation algorithms.
	 * VHP and FP try to return the ball to the centre of the opponents court,
	 * LP tries to just return the ball to the opponents court.
	 *
	 * @param q0 Initial joint positions.
	 * @param filter_ Reference to an input filter (must be initialized in a separate line).
	 * @param flags Flags/options for player class, initialized with c++11 (see player.hpp)
	 *              or through player.cfg file (see sl_interface)
	 */
	Player(const vec7 & q0, EKF & filter, player_flags & flags);

	/**
	 * @brief Deconstructor for the Player class.
	 * Frees the pointer to Optimization classes.
	 */
	~Player();

	/**
	 * @brief Public interface for estimating ball state.
	 *
	 * This interface allows us to test/debug ball state estimation
	 * (which is private).
	 *
	 * @param obs Ball position observation as a 3-vector.
	 * @return Ball state as a 6-vector, if filter is not initialized,
	 * returns the observation as positions and zeroes as velocities.
	 */
	vec6 filt_ball_state(const vec3 & obs);

	/**
	 * @brief If filter is initialized returns true
	 */
	bool filter_is_initialized() const ;

	/**
	 * @brief Method useful for testing performance of different players.
	 *
	 * Using many balls in simulation requires fast resetting
	 * Setting a time threshold as a resetting condition won't work in this case.
	 *
	 */
	void reset_filter(double std_model, double std_noise);

	/** @brief Get players strategy (if exists) */
	void get_strategy(vec2 & ball_des, double & des_land_time);

	/**
	 * @brief Play Table Tennis.
	 *
	 * Main function for playing Table Tennis. Calls one of three different
	 * trajectory generation algorithms (depending on initialization) and
	 * updates the desired joint states when the optimization threads have finished.
	 *
	 * @param qact Actual joint positions, velocities, accelerations.
	 * @param ball_obs Ball observations (positions as 3-vector).
	 * @param qdes Desired joint positions, velocities, accelerations.
	 */
	void play(const joint & qact, const vec3 & ball_obs, joint & qdes);

	/**
	 * @brief Cheat Table Tennis by getting the exact ball state in simulation.
	 *
	 * Similar to play() method, this method receives from the simulator the
	 * exact ball states, which bypasses then the ball estimation method.
	 * Useful for debugging the internal filter used.
	 *
	 * @param qact Actual joint positions, velocities, accelerations.
	 * @param ballstate Ball state (positions AND velocities).
	 * @param qdes Desired joint positions, velocities, accelerations.
	 */
	void cheat(const joint & qact, const vec6 & ballstate, joint & qdes);

};

/**
 * @brief Initialize an EKF
 *
 * Called generally when ball state needs to be reset
 * Useful for passing to Player constructor.
 * @param var_model Process noise multiplier for identity matrix.
 * @param var_noise Observation noise mult. for identity matrix.
 * @param spin Use spin model if true
 * @param out_reject_mult Mult. for outlier rejection.
 * @param topspin Set topspin parameter (NOT state!) for kalman filter prediction
 * if spin is TRUE
 * @return EKF Extended Kalman Filter (state uninitialized!)
 */
EKF init_filter(const double var_model = 0.001,
                const double var_noise = 0.001,
		        const bool spin = false,
		        const double out_reject_mult = 2.0,
		        const double *topspin = nullptr);

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
					EKF & filter);

/**
 * @brief Checks to see if the observation is new (updated)
 *
 * The blobs need to be at least tol apart from each other in distance
 *
 */
bool check_new_obs(const vec3 & obs, double tol);

/**
 * @brief Check to see if we want to reset the filter.
 *
 * Basically if a new ball appears 300 ms later than the last new ball
 * we reset the filter.
 *
 */
bool check_reset_filter(const bool newball, const int verbose, const double threshold);

/**
 * @brief Generate BATCH 3rd order strike + return polynomials.
 *
 * Based on hitting and returning joint state parameters qf,qfdot
 * and hitting time T, calculates the relevant polynomial parameters
 * and generates BATCH polynomial values till time T.
 * @param qf Hitting joint pos.
 * @param qfdot Hitting joint vels.
 * @param T Hitting time
 * @param qact From actual joint state generate the joint des values
 * @param q_rest_des Desired resting posture
 * @param time2return Time to return after hit to rest posture
 * @param Q Generated joint pos.
 * @param Qd Generated joint vel.
 * @param Qdd Generated joint acc.
 */
void generate_strike(const vec7 & qf,
                     const vec7 & qfdot,
                     const double T,
                     const joint & qact,
		             const vec7 & q_rest_des,
		             const double time2return,
		             mat & Q,
		             mat & Qd,
		             mat & Qdd);

/**
 * @brief Generate strike and return traj. incrementally
 *
 * Given polynomial parameters saved in poly,
 * move on to the NEXT desired state only (joint pos,vel,acc).
 * @param poly Polynomial parameters updated in OPTIM classes
 * @param q_rest_des FIXED desired resting posture
 * @param time2return FIXED time to return to rest posture after hit
 * @param t The time passed already following trajectory
 * @param qdes Update pos,vel,acc values of this desired joints structure
 * @return
 */
bool update_next_state(const spline_params & poly,
		               const vec7 & q_rest_des,
		               const double time2return,
		               double & t_poly,
		               joint & qdes);

/**
 * @brief Generate matrix of joint angles, velocities and accelerations
 */
void gen_3rd_poly(const rowvec & times,
                  const vec7 & a3,
                  const vec7 & a2,
                  const vec7 & a1,
                  const vec7 & a0,
                  mat & Q,
                  mat & Qd,
                  mat & Qdd);


/**
 * @brief Set upper and lower bounds on the optimization.
 * First loads the joint limits and then puts some slack
 */
void set_bounds(double *lb, double *ub, double SLACK, double Tmax);

/**
 * @brief Predict ball with the models fed into the filter
 *
 * Number of prediction steps is given by Nmax in racket
 * parameters
 *
 */
void predict_ball(const double & time_pred, mat & balls_pred, EKF & filter);

/**
 * @brief Predict hitting point on the Virtual Hitting Plane (VHP)
 *
 * If the ball is legal (legal detected bounce or legal predicted bounce)
 * and there is enough time (50 ms threshold) predict loc. on VHP.
 *
 * The location of the VHP is given as an argument (vhp-y-location vhpy)
 *
 */
bool predict_hitting_point(const double & vhpy, const bool & check_b,
		                   vec6 & ball_pred, double & time_pred,
		                   EKF & filter, game & game_state);

/**
 * @brief Compute desired racket pos,vel,normals and/or ball positions, vels.
 * Function that calculates a racket strategy : positions, velocities and racket normal
 * for each point on the predicted ball trajectory (ballMat)
 * to return it a desired point (ballLand) at a desired time (landTime)
 *
 */
optim_des calc_racket_strategy(const mat & balls_predicted,
		                       const vec2 & ball_land_des,
		                       const double time_land_des,
							   optim_des & racket_params);

/**
 * @brief Check if the table tennis trial is LEGAL (hence motion planning can be started).
 *
 * If it is exactly one predicted bounce when in awaiting mode
 * (before an actual legal bounce was detected) trial will be legal!
 *
 * If ball has bounced legally bounce, then there should be no more bounces.
 *
 * TODO: no need to check after ball passes table
 * We can turn this check off in the configuration file
 *
 */
bool check_legal_ball(const vec6 & ball_est,
                        const mat & balls_predicted,
                        game & game_state);

/**
 * @brief Checks for legal ball bounce
 * If an incoming ball has bounced before
 * it is declared ILLEGAL (legal_bounce as DATA MEMBER of player class)
 *
 * Bounce variable is static variable of estimate_ball_state method of player class
 * which is reset each time an incoming ball from ball gun is detected.
 *
 * TODO: also consider detecting HIT by robot racket
 * We can turn this off in configuration file
 *
 */
void check_legal_bounce(const vec6 & ball_est, game & game_state);

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
								    const vec3 & ball_land_des,
								    const double time_land_des,
								    optim_des & racket_params);
#endif /* PLAYER_HPP_ */
