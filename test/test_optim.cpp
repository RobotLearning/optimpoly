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

using namespace arma;

/*
 * Initialize robot posture on the right size of the robot
 */
inline void init_right_posture(double* q0) {

	q0[0] = 1.0;
	q0[1] = -0.2;
	q0[2] = -0.1;
	q0[3] = 1.8;
	q0[4] = -1.57;
	q0[5] = 0.1;
	q0[6] = 0.3;
}

/*
 *
 * Here testing NLOPT optimization for VHP player
 *
 */
BOOST_AUTO_TEST_CASE(test_vhp_optim) {

	cout << "Testing VHP Trajectory Optimizer...\n";
	double qf[NDOF], qfdot[NDOF], q0[NDOF], q0dot[NDOF], T;
	double lb[2*NDOF+1], ub[2*NDOF+1];
	double SLACK = 0.01;
	double Tmax = 1.0;

	// update initial parameters from lookup table
	std::cout << "Looking up a random ball entry..." << std::endl;
	arma_rng::set_seed(2);
	//arma_rng::set_seed_random();
	vec::fixed<15> strike_params;
	vec6 ball_state;
	lookup_random_entry(ball_state,strike_params);
	//std::cout << ball_state << std::endl;
	init_right_posture(q0);

	// initialize ball and racket //
	// predict for T_pred seconds
	set_bounds(lb,ub,SLACK,Tmax);

	int N = 1; //1000;
	double** pos = my_matrix(0,NCART,0,N);
	double** vel = my_matrix(0,NCART,0,N);
	double** normal = my_matrix(0,NCART,0,N);
	racketdes racket_params = {pos, vel, normal, DT, N};

	EKF filter = init_filter();
	mat66 P; P.eye();
	filter.set_prior(ball_state,P);

	double time_pred;
	vec6 ball_pred;
	game game_state = AWAITING;
	vec2 ball_land_des = {0.0, dist_to_table - 3*table_length/4};
	double time_land_des = 0.8;
	BOOST_TEST(predict_hitting_point(ball_pred,time_pred,filter,game_state));
	//cout << ball_pred << endl;
	calc_racket_strategy(ball_pred,ball_land_des,time_land_des,racket_params);

	vec3 normal_example;
	for (int i = 0; i < NCART; i++) { // check for normal orthonormality
		normal_example(i) = racket_params.normal[i][0];
	}
	BOOST_TEST(arma::norm(normal_example) == 1.0, boost::test_tools::tolerance(0.01));

	Optim *opt = new HittingPlane(q0,lb,ub);
	opt->fill(&racket_params,q0,q0dot,time_pred);
	opt->run();
	bool update = opt->get_params(qf,qfdot,T);

	BOOST_TEST(update);
	delete opt;
	my_free_matrix(normal,0,NCART,0,racket_params.Nmax);
	my_free_matrix(pos,0,NCART,0,racket_params.Nmax);
	my_free_matrix(vel,0,NCART,0,racket_params.Nmax);
}

/*
 * Testing Fixed Player (or Focused Player)
 */
BOOST_AUTO_TEST_CASE(test_fp_optim) {

	cout << "Testing FP Trajectory Optimizer...\n";
	double qf[NDOF], qfdot[NDOF], q0[NDOF], q0dot[NDOF], T;
	double lb[2*NDOF+1], ub[2*NDOF+1];
	double SLACK = 0.01;
	double Tmax = 1.0;

	// update initial parameters from lookup table
	std::cout << "Looking up a random ball entry..." << std::endl;
	arma_rng::set_seed(2);
	//arma_rng::set_seed_random();
	vec::fixed<15> strike_params;
	vec6 ball_state;
	lookup_random_entry(ball_state,strike_params);
	//std::cout << ball_state << std::endl;
	init_right_posture(q0);

	// initialize ball and racket //
	// predict for T_pred seconds
	set_bounds(lb,ub,SLACK,Tmax);

	int N = 1000;
	double** pos = my_matrix(0,NCART,0,N);
	double** vel = my_matrix(0,NCART,0,N);
	double** normal = my_matrix(0,NCART,0,N);
	racketdes racket_params = {pos, vel, normal, DT, N};

	EKF filter = init_filter();
	mat66 P; P.eye();
	filter.set_prior(ball_state,P);

	vec6 ball_pred;
	double time_land_des = 0.8;
	mat balls_pred = filter.predict_path(DT,N);
	vec2 ball_land_des = {0.0, dist_to_table - 3*table_length/4};
	racket_params = calc_racket_strategy(balls_pred,ball_land_des,time_land_des,racket_params);

	Optim *opt = new FocusedOptim(q0,lb,ub);
	opt->fill(&racket_params,q0,q0dot,0.5);
	opt->run();
	bool update = opt->get_params(qf,qfdot,T);

	BOOST_TEST(update);
	delete opt;
	my_free_matrix(normal,0,NCART,0,racket_params.Nmax);
	my_free_matrix(pos,0,NCART,0,racket_params.Nmax);
	my_free_matrix(vel,0,NCART,0,racket_params.Nmax);
}
