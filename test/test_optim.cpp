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
int randval;

static double cost_fnc(unsigned n, const double *x,
		                   double *grad, void *data);
static void interp_ball(const mat & ballpred,
		                const double T, vec3 & ballpos);

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
 * Measure elapsed timee for spin based racket calculation
 * for a full predicted trajectory. It should be not much more than
 * a spin-less racket strategy computation
 *
 * Compute spin-based desired racket pos,vel,normals and/or ball positions, vels.
 *
 * Outgoing ball velocity is computed using a boundary value problem
 * for the spin model, which is initialized using ballistic model.
 *
 */
BOOST_AUTO_TEST_CASE( test_spin_based_racket_calc ) {

	// first test spin-less case
	//static wall_clock timer;
	//timer.tic();
	arma_rng::set_seed_random();
	double topspin = -50.0;
	TableTennis tt = TableTennis(true,false);
	tt.set_topspin(topspin);
	int ball_launch_side = (randi(1,distr_param(0,2)).at(0));
	tt.set_ball_gun(0.05,ball_launch_side);
	int N = 50;
	mat balls_pred = zeros<mat>(6,N);
	for (int i = 0; i < N; i++) {
		tt.integrate_ball_state(0.02);
		balls_pred.col(i) = tt.get_ball_state();
	}
	vec2 ball_land_des = {0.0, dist_to_table - 3*table_length/4.0};
	double time_land_des = 0.6;
	vec3 ball_vel_out;
	vec6 ballin = tt.get_ball_state();
	tt.calc_des_ball_out_vel(ball_land_des,time_land_des,ballin,ball_vel_out);
	cout << "Incoming ball: " << ballin.t();
	cout << "Estimated outgoing ball velocity: " << ball_vel_out.t();
	optim_outgoing_ball_vel(ball_vel_out);
	cout << "Solution to BVP: " << ball_vel_out.t();

	BOOST_TEST(false);

}
