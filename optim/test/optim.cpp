/*
 * optim.cpp
 *
 * Unit Tests for polynomial optimization
 *
 *  Created on: Feb 17, 2017
 *      Author: okoc
 */

#include <boost/test/unit_test.hpp>
#include <iostream>
#include <armadillo>
#include "constants.h"
#include "kinematics.h"
#include "optimpoly.h"
#include "lookup.h"
#include "player.hpp"

using namespace std;
using namespace arma;

/*
 * Set upper and lower bounds on the optimization.
 * First loads the joint limits and then
 */
inline void set_bounds(double *lb, double *ub, double SLACK, double Tmax) {

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

//BOOST_AUTO_TEST_CASE(test_predict_path) {
//
//	cout << "Testing Robot racket calculations..." << endl;
//	static double pos[NCART] = {1.0, -2.0, -0.5};
//	static double vel[NCART] = {3.0, 5.0, 4.0};
//	EKF filter = init_filter();
//	vec3 ballpos(pos);
//	vec3 ballvel(vel);
//	mat66 P; P.eye();
//	filter.set_prior(join_vert(ballpos,ballvel),P);
//	mat balls_pred = filter.predict_path(dt,10);
//	//cout << "Balls predicted:" << endl << balls_pred << endl;
//
//}

BOOST_AUTO_TEST_CASE(test_nlopt_optim) {

	cout << "Testing NLOPT Optimization" << endl;
	double *q0dot = (double*)calloc(NDOF,sizeof(double));
	double *q0 = (double*)calloc(NDOF,sizeof(double));
	// initial guess for optim //
	double *lb = (double*)calloc(OPTIM_DIM,sizeof(double));
	double *ub = (double*)calloc(OPTIM_DIM,sizeof(double));
	double SLACK = 0.01;
	double Tmax = 1.0;

	// update initial parameters from lookup table
	cout << "Looking up a random entry..." << endl;
	arma_rng::set_seed_random();
	vec::fixed<15> strike_params;
	vec6 ball_state;
	lookup_random_entry(ball_state,strike_params);

	vec7 init_joint_state;
	init_right_posture(init_joint_state);

	// initialize ball and racket //
	// predict for T_pred seconds
	set_bounds(lb,ub,SLACK,Tmax);

	double time2return = 1.0;
	EKF filter = init_filter();
	mat66 P; P.eye();
	filter.set_prior(ball_state,P);
	Player robot = Player(init_joint_state,filter);
	racket racket_params = send_racket_strategy(robot);
	vec3 normal_example;
	int example = 5;
	for (int i = 0; i < NCART; i++) {
		normal_example(i) = racket_params.normal[i][example];
	}
	BOOST_TEST(arma::norm(normal_example) == 1.0, boost::test_tools::tolerance(0.01));

	init_coptim_params(init_joint_state, q0);
	coptim coparams = {q0, q0dot, q0, lb, ub, time2return};
	optim opt_params = {q0, q0dot, 0.5, false};

	// run NLOPT opt algorithm here //
	double max_violation = nlopt_optim_poly_run(&coparams,&racket_params,&opt_params);

	// test to see if kinematics constraints are violated
	BOOST_TEST(max_violation < 0.01);
}
