/*
 * optim.cpp
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
#include "utils.h"
#include "kinematics.h"
#include "optim.h"
#include "lookup.h"
#include "player.hpp"
#include "tabletennis.h"

using namespace arma;

/*
 * Sending to C optimization routine the right data
 */
inline void init_coptim_params(const vec7 & qinit, double *q0) {

	for (int i = 0; i < NDOF; i++) {
		q0[i] = qinit(i);
	}
}

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
 * Here testing NLOPT optimization for FIXED player
 *
 */
BOOST_AUTO_TEST_CASE(test_nlopt_optim) {

	std::cout << "**************Testing NLOPT Optimization***********" << std::endl;
	double *q0dot = (double*)calloc(NDOF,sizeof(double));
	double *q0 = (double*)calloc(NDOF,sizeof(double));
	// initial guess for optim //
	double *lb = (double*)calloc(OPTIM_DIM,sizeof(double));
	double *ub = (double*)calloc(OPTIM_DIM,sizeof(double));
	double SLACK = 0.01;
	double Tmax = 1.0;

	// update initial parameters from lookup table
	std::cout << "Looking up a random entry..." << std::endl;
	arma_rng::set_seed(2);
	//arma_rng::set_seed_random();
	vec::fixed<15> strike_params;
	vec6 ball_state;
	lookup_random_entry(ball_state,strike_params);
	//std::cout << ball_state << std::endl;

	vec7 init_joint_state;
	init_right_posture(init_joint_state);

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
	mat balls_pred = filter.predict_path(DT,N);
	vec2 ball_land_des = {0.0, dist_to_table - 3*table_length/4};
	racket_params = calc_racket_strategy(balls_pred,ball_land_des,0.8,racket_params);
	vec3 normal_example;
	for (int i = 0; i < NCART; i++) { // check for normal orthonormality
		normal_example(i) = racket_params.normal[i][5];
	}
	BOOST_TEST(arma::norm(normal_example) == 1.0, boost::test_tools::tolerance(0.01));

	init_coptim_params(init_joint_state, q0);
	coptim coparams = {q0, q0dot, q0, lb, ub, 1.0};
	optim opt_params = {q0, q0dot, 0.5, false, false};

	// run NLOPT opt algorithm here //
	/*std::thread t(&nlopt_optim_fixed_run,
			&coparams,&racket_params,&opt_params);
	t.join();
	std::cout << "***************************************" << std::endl;*/

	/*auto future = std::async(nlopt_optim_fixed_run,
			&coparams,&racket_params,&opt_params);
    double max_violation = future.get();*/
	double max_violation = nlopt_optim_fixed_run
			(&coparams,&racket_params,&opt_params);

	// test to see if kinematics constraints are violated
	BOOST_TEST(max_violation < 0.01);
}
