/*
 * player.cpp
 *
 * Table Tennis Player Class
 *
 *  Created on: Feb 8, 2017
 *      Author: okoc
 */

#include <armadillo>
#include <thread>
#include <string>
#include "stdlib.h"
#include "player.hpp"
#include "constants.h"
#include "kalman.h"
#include "tabletennis.h"
#include "kinematics.h"
#include "utils.h"
#include "optim.h"

using namespace arma;

/*
 * Initialize Table Tennis Player
 *
 * Incorporates 3 different algorithms
 * VHP and FP try to return the ball to the centre of the opponents court
 * LP tries to just return the ball to the opponents court
 *
 * TODO: remove callocs and replace with new!
 *
 */
Player::Player(const vec7 & q0, EKF & filter_, algo alg_, bool mpc_, int verbose_)
               : filter(filter_), alg(alg_), mpc(mpc_), verbose(verbose_) {

	time_land_des = 0.8;
	time2return = 1.0;
	num_obs = 0;
	validball = true;

	ball_land_des(X) = 0.0;
	ball_land_des(Y) = dist_to_table - 3*table_length/4;
	q_rest_des = q0;

	int N = 500;
	double** pos = my_matrix(0,NCART,0,N);
	double** vel = my_matrix(0,NCART,0,N);
	double** normal = my_matrix(0,NCART,0,N);
	racket_params = {pos, vel, normal, 0.002, N};

	double* qzerodot = (double*)calloc(NDOF,sizeof(double));
	double* qzerodot2 = (double*)calloc(NDOF,sizeof(double));
	double* qzero = (double*)calloc(NDOF, sizeof(double));
	double* qinit = (double*)calloc(NDOF, sizeof(double));
	double* qrest = (double*)calloc(NDOF, sizeof(double));
	double *lb = (double*)calloc(OPTIM_DIM,sizeof(double));
	double *ub = (double*)calloc(OPTIM_DIM,sizeof(double));
	double SLACK = 0.02;
	double Tmax = 1.0;
	set_bounds(lb,ub,SLACK,Tmax);

	for (int i = 0; i < NDOF; i++) {
		qinit[i] = qrest[i] = qzero[i] = q0(i);
	}

	optim_params = {qzero, qzerodot, 0.5, false, false, verbose > 1};
	coparams = {qinit, qzerodot2, qrest, lb, ub, time2return};
	moving = false;

}

/*
 *
 * Deconstructor for the Player class.
 * Frees the memory using free() as in C-style since calloc() was called.
 *
 */
Player::~Player() {

	free(optim_params.qf);
	free(optim_params.qfdot);
	free(coparams.lb);
	free(coparams.ub);
	free(coparams.q0dot);
	free(coparams.qrest);
	free(coparams.q0);
	my_free_matrix(racket_params.normal,0,NCART,0,racket_params.Nmax);
	my_free_matrix(racket_params.pos,0,NCART,0,racket_params.Nmax);
	my_free_matrix(racket_params.vel,0,NCART,0,racket_params.Nmax);
}

/*
 * Filter the blob information with a Kalman Filter.
 * (Extended) KF is used both in simulation mode and for real robot.
 *
 * Checking for new ball that is at least 1 mm away from last observation
 * Checking for also outliers.
 * Resets if the ball suddenly appears on the opponent's court.
 *
 * Returns the valid ball flag.
 *
 * TODO: we're assuming that time elasped dt = DT = 0.002 seconds every time!
 *
 */
void Player::estimate_ball_state(const vec3 & obs) {

	// observation matrix
	static const int min_obs = 5;
	static mat OBS = zeros<mat>(3,min_obs);
	static vec TIMES = zeros<vec>(min_obs);
	static double t_cum;
	bool newball = check_new_obs(obs,1e-3);
	validball = false;

	if (num_obs == 0) { // firsttime
		t_cum = 0.0; // t_cumulative
	}

	if (check_reset_filter(newball,verbose,filter)) {
		num_obs = 0;
		t_cum = 0.0; // t_cumulative
	}

	if (num_obs < min_obs) {
		if (newball) {
			TIMES(num_obs) = t_cum;
			OBS.col(num_obs) = obs;
			num_obs++;
			if (num_obs == min_obs) {
				if (verbose >= 1)
					cout << "Estimating initial ball state\n";
				estimate_prior(OBS,TIMES,filter);
				//cout << OBS << TIMES << filter.get_mean() << endl;
			}
		}
	}
	else { // comes here if there are enough balls to start filter
		filter.predict(DT);
		if (newball) { // && !filter.check_outlier(obs,verbose > 0)) {
			validball = true;
			filter.update(obs);
			//cout << "Updating...\n" << "OBS\n" << obs << "FILT\n" << filter.get_mean() << endl;
		}
	}
	t_cum += DT;
}

/*
 * Public interface, suitable for testing
 *
 */
vec6 Player::filt_ball_state(const vec3 & obs) {

	estimate_ball_state(obs);
	try {
		return filter.get_mean();
	}
	catch (const char * exception) {
		return join_vert(obs,zeros<vec>(3));
	}
}

/*
 * Play Table Tennis
 */
void Player::play(const joint & qact,
		           const vec3 & ball_obs,
				   joint & qdes) {

	estimate_ball_state(ball_obs);

	// initialize optimization and get the hitting parameters
	switch (alg) {
		case FIXED:
			optim_fixedp_param(qact);
			break;
		case VHP:
			optim_vhp_param(qact);
			break;
		case LAZY:
			optim_lazy_param(qact);
			break;
		default:
			throw ("Algorithm is not recognized!\n");
	}

	// generate movement or calculate next desired step
	calc_next_state(qact, qdes);

}

/*
 * Cheat Table Tennis by getting the exact ball state in simulation
 * Useful for debugging filter
 *
 */
void Player::cheat(const joint & qact, const vec6 & ballstate, joint & qdes) {


	filter.set_prior(ballstate,0.01*eye<mat>(6,6));

	switch (alg) {
		case FIXED:
			optim_fixedp_param(qact);
			break;
		case VHP:
			optim_vhp_param(qact);
			break;
		case LAZY:
			optim_lazy_param(qact);
			break;
		default:
			throw ("Algorithm is not recognized!\n");
	}

	// generate movement or calculate next desired step
	calc_next_state(qact, qdes);
}

/*
 * Calculate hitting parameters qf, qfdot
 * on the Virtual Hitting Plane (VHP) by running Inverse Kinematics
 *
 * The inverse kinematics routine runs an optimization to minimize
 * the distance to a rest posture
 *
 *
 */
void Player::optim_vhp_param(const joint & qact) {

	double time_pred;
	vec6 balls_pred;

	// if ball is fast enough and robot is not moving consider optimization
	if (check_update()) {
		if (predict_hitting_point(balls_pred,time_pred)) { // ball is legal and reaches VHP
			calc_racket_strategy(balls_pred,ball_land_des,
					time_land_des,racket_params);
			for (int i = 0; i < NDOF; i++) {
				coparams.q0[i] = qact.q(i);
				coparams.q0dot[i] = qact.qd(i);
			}
			optim_params.T = time_pred;
			// run optimization in another thread
			std::thread t(&nlopt_vhp_run,
					&coparams,&racket_params,&optim_params);
			t.detach();
		}
	}

}

/*
 * FIXED PLAYER
 *
 * Calculate the optimization parameters using an NLOPT nonlinear optimization algorithm
 * in another thread
 *
 * The optimized parameters are: qf, qf_dot, T
 * assuming T_return and q0 are fixed
 *
 */
void Player::optim_fixedp_param(const joint & qact) {

	mat balls_pred;

	// if ball is fast enough and robot is not moving consider optimization
	if (check_update()) {
		predict_ball(balls_pred);
		if (check_legal_ball(balls_pred)) { // ball is legal
			calc_racket_strategy(balls_pred,ball_land_des,time_land_des,racket_params);
			for (int i = 0; i < NDOF; i++) {
				coparams.q0[i] = qact.q(i);
				coparams.q0dot[i] = qact.qd(i);
			}
			//cout << state_est << endl;
			// run optimization in another thread
			std::thread t(&nlopt_optim_fixed_run,
					&coparams,&racket_params,&optim_params);
			t.detach();
		}
		else {
			if (verbose > 1)
				cout << "ball is not legal!" << endl;
		}
	}
}

/*
 * LAZY PLAYER
 *
 * Calculate the optimization parameters using an NLOPT nonlinear optimization algorithm
 * in another thread
 *
 * The optimized parameters are: qf, qf_dot, T
 *
 * TODO: delete ballpred in destructor?
 *
 */
void Player::optim_lazy_param(const joint & qact) {

	static double** ballpred = my_matrix(0,2*NCART,0,racket_params.Nmax);
	mat balls_pred;

	// if ball is fast enough and robot is not moving consider optimization
	if (check_update()) {
		predict_ball(balls_pred);
		if (check_legal_ball(balls_pred)) { // ball is legal
			calc_racket_strategy(balls_pred,ball_land_des,time_land_des,racket_params);
			for (int i = 0; i < NDOF; i++) {
				coparams.q0[i] = qact.q(i);
				coparams.q0dot[i] = qact.qd(i);
			}
			for (int i = 0; i < racket_params.Nmax; i++) {
				for (int j = 0; j < 2*NCART; j++) {
					ballpred[j][i] = balls_pred(j,i);
				}
			}
			// run optimization in another thread
			std::thread t(&nlopt_optim_lazy_run,
					ballpred,&coparams,&racket_params,&optim_params);
			t.detach();
		}
	}


}

/*
 * Check MPC flag and update if possible
 *
 * IF MPC IS TURNED OFF
 * if ball is fast enough and robot is not moving consider optimization
 *
 * IF MPC IS TURNED ON
 * then additionally consider (after running initially optimization)
 * relaunching optimization if ball is valid (new ball and not an outlier)
 * the frequency of updates is respected, and the ball has not passed the y-limit
 *
 */
bool Player::check_update() const {

	vec6 state_est;
	bool update;
	//static int num_updates;
	static const double FREQ_MPC = 20.0;
	static wall_clock timer;
	bool activate, passed_lim = false;

	try {
		state_est = filter.get_mean();
		update = !optim_params.update && !optim_params.running
				&& state_est(Y) > (dist_to_table - table_length/2) ;
		if (mpc && moving) {
			activate = (timer.toc() > (1.0/FREQ_MPC));
			//passed_lim = state_est(Y) > -0.2; //cart_state(Y);
			update = update && validball && activate && !passed_lim;
		}
		else {
			update = update && !moving;
		}
	}
	catch (const char * not_init_error) {
		update = false;
	}
	if (update) {
		//cout << num_updates++ << endl;
		timer.tic();
	}
	return update;
}

/*
 * Predict ball with the models fed into the filter
 *
 * Number of prediction steps is given by Nmax in racket
 * parameters
 *
 */
void Player::predict_ball(mat & balls_pred) const {

	//static wall_clock timer;
	//timer.tic();
	int N = racket_params.Nmax;
	balls_pred = filter.predict_path(racket_params.dt,N);
	//cout << timer.toc() << endl;
}

/*
 * Unfold the next desired state of the 3rd order polynomials in joint space
 * If movement finishes then the desired state velocities and accelerations are zeroed.
 *
 * Multithreading : if after initial lookup, the computation in the other thread terminates, then we
 * synchonize the values between the threads and compute the next desired states given these new polynomial
 * parameters (qf, qf_dot and T_hit)
 *
 */
void Player::calc_next_state(const joint & qact, joint & qdes) {

	static unsigned idx = 0;
	static mat Q_des, Qd_des, Qdd_des;

	// this should be only for MPC?
	if (optim_params.update) {
		if (verbose > 0) {
			std::cout << "Launching/updating strike" << std::endl;
			//for (int i = 0; i < 7; i++)
			//	std::cout << optim_params.qf[i] << "\t" << optim_params.qfdot[i] << "\t";
			//std::cout << optim_params.T << std::endl;
		}
		moving = true;
		idx = 0;
		optim_params.update = false;
		if (alg == LAZY) {
			for (int i = 0; i < NDOF; i++)
				q_rest_des(i) = coparams.qrest[i];
		}
		generate_strike(optim_params,qact,q_rest_des,time2return,Q_des,Qd_des,Qdd_des);
		// call polynomial generation
	}

	// make sure we update after optim finished
	if (moving) {
		qdes.q = Q_des.col(idx);
		qdes.qd = Qd_des.col(idx);
		qdes.qdd = Qdd_des.col(idx);
		idx++;
		if (idx == Q_des.n_cols) {
			// hitting process will finish
			moving = false;
			qdes.q = q_rest_des;
			qdes.qd = zeros<vec>(7);
			qdes.qdd = zeros<vec>(7);
			idx = 0;
		}
	}

}

/*
 * Predict hitting point on the Virtual Hitting Plane (VHP)
 *
 * The location of the VHP is defined as a constant (constants.h)
 *
 */
bool Player::predict_hitting_point(vec6 & ball_pred, double & time_pred) const {

	bool valid_hp = false;
	mat balls_path;
	predict_ball(balls_path);
	uvec vhp_index;
	unsigned idx;

	if (check_legal_ball(balls_path)) {
		vhp_index = find(balls_path.row(Y) >= VHPY, 1);
		if (vhp_index.n_elem == 1) {
			idx = as_scalar(vhp_index);
			ball_pred = balls_path.col(idx);
			time_pred = racket_params.dt * (idx + 1);
			valid_hp = true;
		}
	}

	return valid_hp;
}

/*
 * Function that calculates a racket strategy : positions, velocities and racket normal
 * for each point on the predicted ball trajectory (ballMat)
 * to return it a desired point (ballLand) at a desired time (landTime)
 *
 */
racketdes calc_racket_strategy(const mat & balls_predicted,
		                       const vec2 & ball_land_des, const double time_land_des,
							   racketdes & racket_params) {

	//static wall_clock timer;
	//timer.tic();

	int N = balls_predicted.n_cols;
	mat balls_out_vel = zeros<mat>(3,N);
	mat racket_des_normal = zeros<mat>(3,N);
	mat racket_des_vel = zeros<mat>(3,N);
	calc_des_ball_out_vel(ball_land_des,time_land_des,balls_predicted,balls_out_vel);
	calc_des_racket_normal(balls_predicted.rows(DX,DZ),balls_out_vel,racket_des_normal);
	calc_des_racket_vel(balls_predicted.rows(DX,DZ),balls_out_vel,racket_des_normal,racket_des_vel);

	// place racket centre on the predicted ball

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < NCART; j++) {
			racket_params.pos[j][i] = balls_predicted(j,i);
			racket_params.vel[j][i] = racket_des_vel(j,i);
			racket_params.normal[j][i] = racket_des_normal(j,i);
		}
	}
	//cout << timer.toc() << endl;
	return racket_params;
}

/*
 * Create batch hitting and returning joint state 3rd degree polynomials
 *
 */
void generate_strike(const optim & params, const joint & qact,
		             const vec7 & q_rest_des, const double time2return,
		            mat & Q, mat & Qd, mat & Qdd) {

	// first create hitting polynomials
	vec7 a2, a3;
	vec7 b2, b3; // for returning
	double T = params.T;
	vec7 qf(params.qf);
	vec7 qfdot(params.qfdot);
	vec7 qnow = qact.q;
	vec7 qdnow = qact.qd;
	a3 = 2.0 * (qnow - qf) / pow(T,3) + (qfdot + qdnow) / pow(T,2);
	a2 = 3.0 * (qf - qnow) / pow(T,2) - (qfdot + 2.0*qdnow) / T;
	b3 = 2.0 * (qf - q_rest_des) / pow(time2return,3) + (qfdot) / pow(time2return,2);
	b2 = 3.0 * (q_rest_des - qf) / pow(time2return,2) - (2.0*qfdot) / time2return;

	int N_hit = T/DT;
	rowvec times_hit = linspace<rowvec>(DT,T,N_hit);
	int N_return = time2return/DT;
	rowvec times_ret = linspace<rowvec>(DT,time2return,N_return);

	mat Q_hit, Qd_hit, Qdd_hit, Q_ret, Qd_ret, Qdd_ret;
	Q_hit = Qd_hit = Qdd_hit = zeros<mat>(NDOF,N_hit);
	Q_ret = Qd_ret = Qdd_ret = zeros<mat>(NDOF,N_return);

	gen_3rd_poly(times_hit,a3,a2,qdnow,qnow,Q_hit,Qd_hit,Qdd_hit);
	gen_3rd_poly(times_ret,b3,b2,qfdot,qf,Q_ret,Qd_ret,Qdd_ret);
	Q = join_horiz(Q_hit,Q_ret);
	Qd = join_horiz(Qd_hit,Qd_ret);
	Qdd = join_horiz(Qdd_hit,Qdd_ret);
}

/*
 * Initialize an Extended Kalman Filter
 * useful for passing to Player constructor
 *
 */
EKF init_filter() {

	double std_model = 0.3;
	double std_noise = 0.1;
	mat C = eye<mat>(3,6);
	mat66 Q = std_model * eye<mat>(6,6);
	mat33 R = std_noise * eye<mat>(3,3);
	return EKF(calc_next_ball,C,Q,R);
}

/*
 * Generate matrix of joint angles, velocities and accelerations
 */
void gen_3rd_poly(const rowvec & times, const vec7 & a3, const vec7 & a2, const vec7 & a1, const vec7 & a0,
		     mat & Q, mat & Qd, mat & Qdd) {

	// IN MATLAB:
	//	qStrike(m,:) = a(1)*t.^3 + a(2)*t.^2 + a(3)*t + a(4);
	//	qdStrike(m,:) = 3*a(1)*t.^2 + 2*a(2)*t + a(3);
	//	qddStrike(m,:) = 6*a(1)*t + 2*a(2);

	for(int i = 0; i < NDOF; i++) {
		Q.row(i) = a3(i) * pow(times,3) + a2(i) * pow(times,2) + a1(i) * times + a0(i);
		Qd.row(i) = 3*a3(i) * pow(times,2) + 2*a2(i) * times + a1(i);
		Qdd.row(i) = 6*a3(i) * times + 2*a2(i);
	}
}

/*
 * Calculate desired racket normal assuming mirror law
 */
void calc_des_racket_normal(const mat & v_in, const mat & v_out, mat & normal) {

	normal = v_out - v_in;
	// normalize
	normal = normalise(normal);
}

/*
 *
 *  Computes the desired outgoing velocity of the ball after contact
 *	to hit the goal on a desired landing position on the
 *	opponents court
 *
 *
 */
void calc_des_ball_out_vel(const vec2 & ball_land_des,
						   const double time_land_des,
						   const mat & balls_predicted, mat & balls_out_vel) {

	static double z_table = floor_level - table_height + ball_radius;

	// elementwise division
	balls_out_vel.row(X) = (ball_land_des(X) - balls_predicted.row(X)) / time_land_des;
	balls_out_vel.row(Y) = (ball_land_des(Y) - balls_predicted.row(Y)) / time_land_des;
	balls_out_vel.row(Z) = (z_table - balls_predicted.row(Z) -
			                0.5 * gravity * pow(time_land_des,2)) / time_land_des;

	//TODO: consider air drag, hack for now
	balls_out_vel.row(X) *= 1.1;
	balls_out_vel.row(Y) *= 1.1;
	balls_out_vel.row(Z) *= 1.2;
}

/*
 * Calculate desired racket velocity given ball incoming and
 * outgoing velocities
 * Assuming a mirror law
 * Assumes no desired spin, i.e. racket velocity along the racket will be set to zero
 *
 * Output is the last parameter: racketVel
 *
 */
void calc_des_racket_vel(const mat & vel_ball_in, const mat & vel_ball_out,
		                 const mat & racket_normal, mat & racket_vel) {

	int N = vel_ball_in.n_cols;
	for (int i = 0; i < N; i++) {
		racket_vel.col(i) = dot(((vel_ball_out.col(i) + CRR * vel_ball_in.col(i)) / (1 + CRR)),
								racket_normal.col(i)) * racket_normal.col(i);
	}
}

/*
 * Check if there is a bounce among the predicted ball states
 * If it is exactly one ball will be legal!
 *
 * TODO: no need to check after ball passes table
 *
 */
bool check_legal_ball(const mat & balls_predicted) {

	int num_bounces = 0;
	int N = balls_predicted.n_cols;

	// if sign of z-velocity changes then the ball bounces
	for (int i = 0; i < N-1; i++) {
		if (balls_predicted(DZ,i) < 0 && balls_predicted(DZ,i+1) > 0) {
			num_bounces++;
		}
	}

	// multiple bounces are predicted
	if (num_bounces > 1) {
		//std::cout << "Multiple bounces predicted. Not moving" << endl;
		return false;
	}
	// no bounce is predicted
	if (num_bounces == 0) {
		//std::cout << "No bounce predicted. Not moving\n" << endl;
		return false;
	}

	return true;
}

/*
 * Checks to see if the observation is new (updated)
 *
 * The blobs need to be at least tol = 1e-3 apart from each other in distance
 *
 */
bool check_new_obs(const vec3 & obs, double tol) {

	static vec3 last_obs = zeros<vec>(3);

	if (norm(obs - last_obs) > tol) {
		last_obs = obs;
		return true;
	}
	return false;
}

/*
 * Empirical Bayes procedure to estimate prior given
 * matrix of observations Y of column length N, each row
 * corresponding to one
 *
 * If observations arrive as column vectors then we take
 * transpose of it.
 *
 * Velocity estimation is biased, we multiply velocities by 1.1
 * since they often underestimate actual velocities.
 *
 *
 */
void estimate_prior(const mat & observations,
		            const vec & times,
					EKF & filter) {

	vec6 x; mat66 P;
	int num_samples = times.n_elem;
	mat M = zeros<mat>(num_samples,3);
	vec3 vel_multiplier = {1.1, 1.1, 1.1};

	// and create the data matrix
	for (int i = 0; i < num_samples; i++) {
		M(i,0) = 1.0;
		M(i,1) = times(i);
		M(i,2) = times(i) * times(i);
	}
	// solving for the parameters
	//cout << "Data matrix:" << endl << M << endl;
	mat Beta = solve(M,observations.t());
	//cout << "Parameters:" << endl << Beta << endl;
	x = join_horiz(Beta.row(0),Beta.row(1)).t(); //vectorise(Beta.rows(0,1));
	P.eye(6,6);
	//P *= 0.1;
	filter.set_prior(x,P);
	filter.update(observations.col(0));

	double dt;
	for (unsigned i = 1; i < times.n_elem; i++) {
		dt = times(i) - times(i-1);
		filter.predict(dt);
		filter.update(observations.col(i));
	}
	//x = filter.get_mean();
	//x(span(DX,DZ)) = x(span(DX,DZ)) % vel_multiplier;
	//filter.set_prior(x,P);
}

/*
 * Check to see if we want to reset the filter.
 *
 * Basically if a new ball appears 300 ms later than the last new ball
 * we reset the filter.
 *
 * TODO: check if it works on real ball data!
 *
 */
bool check_reset_filter(const bool newball, const int verbose, EKF & filter) {

	bool reset = false;
	static int reset_cnt = 0;
	static double threshold = 0.3; // 300 miliseconds
	static bool firsttime = true;
	static wall_clock timer;

	if (firsttime) {
		firsttime = false;
		timer.tic();
	}

	if (newball) {
		if (timer.toc() > threshold) {
			reset = true;
			if (verbose > 0) {
				std::cout << "Resetting filter! Count: " << ++reset_cnt << std::endl;
			}
			filter = init_filter();
		}
		timer.tic();
	}

	return reset;
}

/*
 * Set upper and lower bounds on the optimization.
 * First loads the joint limits and then
 */
void set_bounds(double *lb, double *ub, double SLACK, double Tmax) {

	read_joint_limits(lb,ub);
	// lower bounds and upper bounds for qf are the joint limits
	for (int i = 0; i < NDOF; i++) {
		ub[i] -= SLACK;
		lb[i] += SLACK;
		ub[i+NDOF] = MAX_VEL;
		lb[i+NDOF] = -MAX_VEL;
	}
	// constraints on final time
	ub[2*NDOF] = Tmax;
	lb[2*NDOF] = 0.0;
}
