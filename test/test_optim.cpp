/*
 * test_optim.cpp
 *
 * Unit Tests for polynomial optimization
 *
 *  Created on: Feb 17, 2017
 *      Author: okoc
 */

#ifndef BOOST_TEST_MODULE
#define BOOST_TEST_MODULE test_table_tennis
#endif

#include <boost/test/unit_test.hpp>
#include <iostream>
#include <armadillo>
#include <thread>
#include "kinematics.h"
#include "utils.h"
#include "optim.h"
#include "lookup.h"
#include "player.hpp"
#include "tabletennis.h"
#include "kinematics.hpp"

using namespace arma;
using namespace optim;
using namespace player;
int randval;

static double cost_fnc(unsigned n, const double *x,
		                   double *grad, void *data);
static void interp_ball(const mat & ballpred,
		                const double T, vec3 & ballpos);
/**
 * @brief Initial ball data time stamps and position observations
 */
struct des_ball_data {
	vec3 ball_incoming;
	vec3 ball_land_des;
	double time_land_des;
	double topspin;
};

static double calc_landing_res(unsigned n, const double *x, double *grad, void *data);
static void optim_spin_outgoing_ball_vel(const des_ball_data & data, const bool verbose, vec3 & est); // spin based optimization


/*
 * Initialize robot posture on the right size of the robot
 */
inline void init_right_posture(vec7 & q0) {

	q0(0) = 1.0;
	q0(1) = -0.2;
	q0(2) = -0.1;
	q0(3) = 1.8;
	q0(4) = -1.57;
	q0(5) = 0.1;
	q0(6) = 0.3;
}

/*
 * Initialize robot posture
 */
void init_posture(vec7 & q0, int posture, bool verbose) {

	rowvec qinit;
	switch (posture) {
	case 2: // right
		if (verbose)
			cout << "Initializing robot on the right side.\n";
		qinit << 1.0 << -0.2 << -0.1 << 1.8 << -1.57 << 0.1 << 0.3 << endr;
		break;
	case 1: // center
		if (verbose)
			cout << "Initializing robot on the center.\n";
		qinit << 0.0 << 0.0 << 0.0 << 1.5 << -1.75 << 0.0 << 0.0 << endr;
		break;
	case 0: // left
		if (verbose)
			cout << "Initializing robot on the left side\n";
		qinit << -1.0 << 0.0 << 0.0 << 1.5 << -1.57 << 0.1 << 0.3 << endr;
		break;
	default: // default is the right side
		qinit << 1.0 << -0.2 << -0.1 << 1.8 << -1.57 << 0.1 << 0.3 << endr;
		break;
	}
	q0 = qinit.t();
}

/*
 *
 * Here testing NLOPT optimization for VHP player
 *
 */
BOOST_AUTO_TEST_CASE(test_vhp_optim) {

	cout << "Testing VHP Trajectory Optimizer...\n";
	const double VHPY = -0.3;
	double lb[2*NDOF+1], ub[2*NDOF+1];
	double SLACK = 0.01;
	double Tmax = 1.0;
	joint qact;
	spline_params poly;

	// update initial parameters from lookup table
	std::cout << "Looking up a random ball entry..." << std::endl;
	arma_rng::set_seed_random();
	randval = (randi(1).at(0));
	arma_rng::set_seed(randval);
	vec::fixed<15> strike_params;
	vec6 ball_state;
	lookup_random_entry(ball_state,strike_params);
	//std::cout << ball_state << std::endl;
	init_right_posture(qact.q);
	set_bounds(lb,ub,SLACK,Tmax);

	EKF filter = init_filter();
	mat66 P; P.eye();
	filter.set_prior(ball_state,P);

	double time_pred;
	vec6 ball_pred;
	game game_state = AWAITING;
	vec2 ball_land_des = {0.0, dist_to_table - 3*table_length/4};
	double time_land_des = 0.8;
	BOOST_TEST(predict_hitting_point(VHPY,true,ball_pred,time_pred,filter,game_state));
	//cout << ball_pred << endl;
	optim_des racket_params;
	calc_racket_strategy(ball_pred,ball_land_des,time_land_des,racket_params);

	vec3 normal_example = racket_params.racket_normal(span(X,Z),0);
	BOOST_TEST(arma::norm(normal_example) == 1.0, boost::test_tools::tolerance(0.01));

	HittingPlane opt = HittingPlane(qact.q.memptr(),lb,ub);
	opt.set_des_params(&racket_params);
	opt.fix_hitting_time(time_pred);
	opt.update_init_state(qact);
	opt.run();
	bool update = opt.get_params(qact,poly);

	BOOST_TEST(update);
}

/*
 * Testing Fixed Player (or Focused Player)
 */
BOOST_AUTO_TEST_CASE(test_fp_optim) {

	cout << "Testing FP Trajectory Optimizer...\n";
	double lb[2*NDOF+1], ub[2*NDOF+1];
	double SLACK = 0.01;
	double Tmax = 1.0;
	joint qact;
	spline_params poly;

	// update initial parameters from lookup table
	std::cout << "Looking up a random ball entry..." << std::endl;
	arma_rng::set_seed(randval);
	//arma_rng::set_seed_random();
	vec::fixed<15> strike_params;
	vec6 ball_state;
	lookup_random_entry(ball_state,strike_params);
	init_right_posture(qact.q);
	set_bounds(lb,ub,SLACK,Tmax);
	optim_des racket_params;
	int N = 1000;
	racket_params.Nmax = 1000;

	EKF filter = init_filter();
	mat66 P; P.eye();
	filter.set_prior(ball_state,P);

	vec6 ball_pred;
	double time_land_des = 0.8;
	mat balls_pred = filter.predict_path(DT,N);
	vec2 ball_land_des = {0.0, dist_to_table - 3*table_length/4};
	racket_params = calc_racket_strategy(balls_pred,ball_land_des,time_land_des,racket_params);

	FocusedOptim opt = FocusedOptim(qact.q.memptr(),lb,ub);
	opt.set_des_params(&racket_params);
	opt.update_init_state(qact);
	opt.run();
	bool update = opt.get_params(qact,poly);

	BOOST_TEST(update);
}

/*
 * Testing Lazy Player (or Defensive Player)
 */
BOOST_AUTO_TEST_CASE(test_dp_optim) {

	cout << "Testing LAZY Trajectory Optimizer...\n";
	double lb[2*NDOF+1], ub[2*NDOF+1];
	double SLACK = 0.01;
	double Tmax = 1.0;
	joint qact;
	spline_params poly;

	// update initial parameters from lookup table
	std::cout << "Looking up a random ball entry..." << std::endl;
	arma_rng::set_seed(randval);
	//arma_rng::set_seed_random();
	vec::fixed<15> strike_params;
	vec6 ball_state;
	lookup_random_entry(ball_state,strike_params);
	init_right_posture(qact.q);
	set_bounds(lb,ub,SLACK,Tmax);

	int N = 1000;
	EKF filter = init_filter();
	mat66 P; P.eye();
	filter.set_prior(ball_state,P);
	mat balls_pred = filter.predict_path(DT,N);
	optim_des ball_params;
	ball_params.ball_pos = balls_pred.rows(X,Z);
	ball_params.ball_vel = balls_pred.rows(DX,DZ);
	ball_params.Nmax = N;
	bool land = true;
	bool lookup = true;
	DefensiveOptim opt = DefensiveOptim(qact.q.memptr(),lb,ub,land,lookup); //only touch the ball if false!
	opt.set_des_params(&ball_params);
	opt.update_init_state(qact);
	opt.set_verbose(true);
	opt.run();
	bool update = opt.get_params(qact,poly);

	BOOST_TEST(update);
}

/**
 * Evaluate time efficiency of FP and DP optimizations over hundreds of tests
 */
/*
BOOST_AUTO_TEST_CASE( test_time_efficiency ) {

	// For FP or DP
	// initialize ball randomly each time
	// measure time elapsed for each and record in a histogram
	int num_trials = 500;
	BOOST_TEST_MESSAGE("Testing time efficiency of FP & DP on " << num_trials << " instances.");
	double Tmax = 1.0, lb[2*NDOF+1], ub[2*NDOF+1];
	set_bounds(lb,ub,0.01,Tmax);
	joint qact;
	spline_params poly;
	TableTennis tt = TableTennis(true,false);
	arma_rng::set_seed_random();
	double std_noise = 0.0001;
	double std_model = 0.3;
	int ball_launch_side;
	int joint_init_pose;
	optim_des racket_params;
	int N = 1000;
	racket_params.Nmax = 1000;
	double time_land_des = 0.8;
	vec2 ball_land_des = {0.0, dist_to_table - 3*table_length/4};
	mat balls_pred = zeros<mat>(6,N);
	vec time_elapsed_init = zeros<vec>(num_trials);
	vec time_elapsed_run_fp = zeros<vec>(num_trials);
	vec time_elapsed_dp_lookup = zeros<vec>(num_trials);
	vec time_elapsed_dp_no_lookup = zeros<vec>(num_trials);
    wall_clock timer;

	for (int n = 0; n < num_trials; n++) { // for each trial
		std::cout << "Trial: " << n+1 << std::endl;
		ball_launch_side = (randi(1,distr_param(0,2)).at(0));
		joint_init_pose = (randi(1,distr_param(0,2)).at(0));
		init_posture(qact.q,joint_init_pose,false);
		tt.reset_stats();
		tt.set_ball_gun(0.05,ball_launch_side);
		for (int i = 0; i < N; i++) {
			tt.integrate_ball_state(DT);
			balls_pred.col(i) = tt.get_ball_state();
		}
		timer.tic();
		racket_params = calc_racket_strategy(balls_pred,ball_land_des,time_land_des,racket_params);
		FocusedOptim fp = FocusedOptim(qact.q.memptr(),lb,ub);
		fp.set_verbose(false);
		fp.set_des_params(&racket_params);
		fp.update_init_state(qact);
		time_elapsed_init(n) = timer.toc() * 1e3;
		timer.tic();
		fp.run();
		time_elapsed_run_fp(n) = timer.toc() * 1e3;
		DefensiveOptim dp = DefensiveOptim(qact.q.memptr(),lb,ub,true,false);
		dp.set_verbose(false);
		dp.set_des_params(&racket_params);
		dp.update_init_state(qact);
		timer.tic();
		dp.run();
		time_elapsed_dp_no_lookup(n) = timer.toc() * 1e3;
		DefensiveOptim dp2 = DefensiveOptim(qact.q.memptr(),lb,ub,true,true);
		dp2.set_verbose(false);
		dp2.set_des_params(&racket_params);
		dp2.update_init_state(qact);
		timer.tic();
		dp2.run();
		time_elapsed_dp_lookup(n) = timer.toc() * 1e3;
	}
	umat hists = zeros<umat>(15,4);
	hists.col(0) = hist(time_elapsed_init, linspace<vec>(0,9,15));
	hists.col(1) = hist(time_elapsed_run_fp, linspace<vec>(5,75,15));
	hists.col(2) = hist(time_elapsed_dp_no_lookup, linspace<vec>(5,75,15));
	hists.col(3) = hist(time_elapsed_dp_lookup, linspace<vec>(5,75,15));
	hists.save("histograms.txt", raw_ascii);
}
*/

/*
 * Find a qf and t such that J(qf) has minimal Frobenius norm
 * while being close to predicted ball state b(t) and being close to q_hit
 */
BOOST_AUTO_TEST_CASE(find_rest_posture) {

	cout << "\nOptimizing resting posture based on jacobian...\n";
	double lb[2*NDOF+1], ub[2*NDOF+1];
	double SLACK = 0.01;
	double Tmax = 1.0;
	vec7 q_rest_des = zeros<vec>(7);
	joint qact;
	// update initial parameters from lookup table
	std::cout << "Looking up a random ball entry..." << std::endl;
	arma_rng::set_seed_random();
	vec::fixed<15> strike_params;
	vec6 ball_state;
	lookup_random_entry(ball_state,strike_params);
	vec7 q_hit = strike_params.head(NDOF);

	set_bounds(lb,ub,SLACK,Tmax);
	int N = 1000;
	EKF filter = init_filter();
	mat66 P; P.eye();
	filter.set_prior(ball_state,P);
	mat balls_pred = filter.predict_path(DT,N);
	optim_des ball_params;
	ball_params.ball_pos = balls_pred.rows(X,Z);
	ball_params.ball_vel = balls_pred.rows(DX,DZ);
	ball_params.Nmax = N;
	DefensiveOptim opt = DefensiveOptim(q_rest_des,lb,ub,true,true);
	opt.set_verbose(false);
	opt.set_des_params(&ball_params);
	opt.update_init_state(qact);
	opt.run();
	opt.run_qrest_optim(q_rest_des);
}

/**
 * For the Focused Player with a spin model and spin estimation:
 *
 * Check accuracy of optimization based racket calculation
 * solving a BVP with spin model
 *
 * comparing with the ballistic model
 *
 */
BOOST_AUTO_TEST_CASE( check_accuracy_spin_based_racket_calc ) {

	// first test spin-less case
	//static wall_clock timer;
	//timer.tic();
	cout << "\nSolving BVP for one particular ball's desired outgoing velocity using a spin model\n";
	cout << "Optimization should produce accurate inversion...\n";
	arma_rng::set_seed_random();
	double topspin = -50.0;
	TableTennis tt = TableTennis(true,false);
	tt.set_topspin(topspin);
	int ball_launch_side = (randi(1,distr_param(0,2)).at(0));
	tt.set_ball_gun(0.05,ball_launch_side);
	int N = 50;
	double dt = 0.02;
	//mat balls_pred = zeros<mat>(6,N);
	for (int i = 0; i < N; i++) {
		tt.integrate_ball_state(dt);
		//balls_pred.col(i) = tt.get_ball_state();
	}
	vec3 ball_land_des = {0.0, dist_to_table - 3*table_length/4.0, floor_level - table_height + ball_radius};
	double time_land_des = 0.6;
	vec3 ball_vel_out;
	vec6 ballin = tt.get_ball_state();
	tt.calc_des_ball_out_vel(ball_land_des.head(2),time_land_des,true,ballin,ball_vel_out);
	cout << "Incoming ball: " << ballin.t();
	cout << "Estimated outgoing ball velocity with ballistic model: " << ball_vel_out.t();

	int N_horizon = time_land_des/dt;
	tt.set_ball_state(join_vert(ballin.head(3),ball_vel_out));
	for (int i = 0; i < N_horizon; i++) {
		tt.integrate_ball_state(dt);
	}
	cout << "Resulting landing position: " << tt.get_ball_position().t();
	double res1 = norm(tt.get_ball_position() - ball_land_des,2);
	cout << "Norm of error :" << res1 << endl;
	// initialize optimization with adjusted values (e.g. to reduce optim time)
	ball_vel_out(X) *= 1.1;
	ball_vel_out(X) *= 1.2;
	ball_vel_out(X) *= 0.9;

	des_ball_data data;
	data.ball_land_des = ball_land_des;
	data.time_land_des = time_land_des;
	data.topspin = topspin;
	data.ball_incoming = ballin.head(3);

	optim_spin_outgoing_ball_vel(data,true,ball_vel_out);
	cout << "Solution to BVP: " << ball_vel_out.t() << endl;

	tt.set_ball_state(join_vert(ballin.head(3),ball_vel_out));
	for (int i = 0; i < N_horizon; i++) {
		tt.integrate_ball_state(dt);
	}
	cout << "Resulting landing position: " << tt.get_ball_position().t();
	double res2 = norm(tt.get_ball_position() - ball_land_des,2);
	cout << "Norm of error :" << res2 << endl;

	BOOST_TEST(res2 < res1);

}


/**
 * @brief Solve BVP for a particular predicted spinning ball's outgoing desired ball velocity
 *
 * BVP is solved using optimization
 */
void optim_spin_outgoing_ball_vel(const des_ball_data & data, const bool verbose, vec3 & est) {

	static double x[3];  /* some initial guess */
	static double minf; /* the minimum objective value, upon return */
	static double init_time;
	static int res; // error code
	static nlopt_opt opt;

	opt = nlopt_create(NLOPT_LD_MMA, 3);
	nlopt_set_min_objective(opt, calc_landing_res, (void*)&data);
	nlopt_set_xtol_rel(opt, 1e-2);

	for(int i = 0; i < 3; i++) {
		x[i] = est(i);
	}

	init_time = get_time();
	if ((res = nlopt_optimize(opt, x, &minf)) < 0) {
		if (verbose)
			printf("NLOPT failed!\n");
	}
	else {
		if (verbose) {
			printf("NLOPT took %f ms\n", (get_time() - init_time)/1e3);
			printf("Found minimum at f = %0.10g\n", minf);
		}
		for(int i = 0; i < 3; i++) {
			est(i) = x[i];
		}
	}
	//nlopt_destroy(opt);
}

/**
 * Cost function for computing the residual (norm squared)
 * of the outgoing ball landing error
 * Calculates also the gradient if grad is TRUE
 *
 */
double calc_landing_res(unsigned n, const double *x, double *grad, void *data) {

	static double dt = 0.02;
    static TableTennis tt = TableTennis(true,false,false); // no contact checking
    static vec3 vel_out;
    static vec6 init_state;
    static vec3 out_pos;

    des_ball_data *mydata = (des_ball_data*) data;
    tt.set_topspin(mydata->topspin);
    vel_out(X) = x[0];
    vel_out(Y) = x[1];
    vel_out(Z) = x[2];
    init_state = join_vert(mydata->ball_incoming,vel_out);
    tt.set_ball_state(init_state);
	for (int i = 0; i < mydata->time_land_des/dt; i++)
    	tt.integrate_ball_state(dt);
	//for (int i = 0; i < 5; i++)
	//	tt.symplectic_int_fourth(mydata->time_land_des/5.0);

	out_pos = tt.get_ball_position();

    if (grad) {

    	grad[0] = mydata->time_land_des * (out_pos(X) - mydata->ball_land_des(X));
    	grad[1] = mydata->time_land_des * (out_pos(Y) - mydata->ball_land_des(Y));
    	grad[2] = mydata->time_land_des * (out_pos(Z) - mydata->ball_land_des(Z));
    	// Finite difference
		/*static double h = 1e-6;
		static double val_plus, val_minus;
		static double xx[3];
		for (unsigned i = 0; i < n; i++)
			xx[i] = x[i];
		for (unsigned i = 0; i < n; i++) {
			xx[i] += h;
			val_plus = calc_landing_res(n, xx, NULL, data);
			xx[i] -= 2*h;
			val_minus = calc_landing_res(n, xx, NULL, data);
			grad[i] = (val_plus - val_minus) / (2*h);
			xx[i] += h;
		}*/
    }

	return pow(norm(out_pos - mydata->ball_land_des),2);
}


/**
 * Testing speed of racket calculations with spin-based BVP solver (optimization)
 */
BOOST_AUTO_TEST_CASE( check_speed_spin_based_racket_calc ) {

	BOOST_TEST_MESSAGE("Checking speed of spin based racket calculation on predicted ball traj.");

	arma_rng::set_seed_random();
	double topspin = -50.0;
	TableTennis tt = TableTennis(true,false);
	tt.set_topspin(topspin);
	int ball_launch_side = (randi(1,distr_param(0,2)).at(0));
	tt.set_ball_gun(0.05,ball_launch_side);
	int N = 50;
	double dt = 0.02;
	mat balls_pred = zeros<mat>(6,N);
	for (int i = 0; i < N; i++) {
		tt.integrate_ball_state(dt);
		balls_pred.col(i) = tt.get_ball_state();
	}
	vec3 ball_land_des = {0.0, dist_to_table - 3*table_length/4.0, floor_level - table_height + ball_radius};
	double time_land_des = 0.6;

	optim_des racket_des;
	racket_des.Nmax = N;
	static wall_clock timer;
	timer.tic();
	calc_racket_strategy(balls_pred,ball_land_des.head(2),time_land_des,racket_des);
	double time1 = timer.toc() * 1e3;
	cout << "Elapsed time in ms: " << time1 << endl;
	timer.tic();
	calc_spin_racket_strategy(balls_pred,topspin,ball_land_des,time_land_des,racket_des);
	double time2 = timer.toc() * 1e3;
	cout << "Elapsed time in ms: " << time2 << endl;
	BOOST_TEST(time2 < time1 * 100);
}

/*
 * Testing 4th order Symplectic Integration
 * It should be more accurate with fewer cycles
 *
 */
BOOST_AUTO_TEST_CASE( test_symplectic_int_4th) {

	BOOST_TEST_MESSAGE("Checking 4th order Symplectic Integration");
	arma_rng::set_seed_random();
	double topspin = -50.0;
	TableTennis tt = TableTennis(true,false);
	tt.set_topspin(topspin);
	int ball_launch_side = (randi(1,distr_param(0,2)).at(0));

	tt.set_ball_gun(0.05,ball_launch_side);
	vec6 init_state = tt.get_ball_state();
	int N = 25; // it should not bounce
	double dt = 0.02;
	for (int i = 0; i < N; i++) {
		tt.integrate_ball_state(dt);
	}
	vec6 ball_pred1 = tt.get_ball_state();

	tt.set_ball_state(init_state);
	for (int i = 0; i < 5; i++)
		tt.symplectic_int_fourth(0.1);
	//tt.symplectic_int_fourth(0.5);
	vec6 ball_pred2 = tt.get_ball_state();

	cout << ball_pred1.t();
	cout << ball_pred2.t();
	BOOST_TEST_MESSAGE("Nothing enforced here!");
}
