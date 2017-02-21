/*
 * optim.cpp
 *
 * Unit Tests for polynomial optimization
 *
 *  Created on: Feb 17, 2017
 *      Author: okoc
 */

#define BOOST_TEST_MODULE test_optim
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <armadillo>
#include "constants.h"
#include "kinematics.h"
#include "optimpoly.h"
#include "utils.h"
#include "stdlib.h"
#include "player.hpp"

using namespace std;
using namespace arma;

#define LOOKUP_TABLE_SIZE 4002 //769
#define LOOKUP_COLUMN_SIZE 2*NDOF + 1 + 2*NCART // ball state and optimization parameters (6 + 15)
#define LOOKUP_TABLE_NAME "LookupTable-16-May-2016" //"LookupTable-March-2016" //"LookupTable-April-2016"

/*
 * Load the lookup table of coparameters and main parameters of interest
 */
inline void load_lookup_table(mat & lookup) {

	lookup = zeros<mat>(LOOKUP_TABLE_SIZE * LOOKUP_COLUMN_SIZE, 1);
	string env = getenv("HOME");
	string filename = env + "/robolab/barrett/saveData/" +
			          LOOKUP_TABLE_NAME + ".txt";

	lookup.load(filename);
	lookup.reshape(LOOKUP_TABLE_SIZE,LOOKUP_COLUMN_SIZE);
	//cout << lookup(span(0,5),span(0,5)) << endl;
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
 *
 * K-nearest-neighbours method for looking up optimization values
 *
 * Find the closest ball states in the lookup table
 * and lookup the corresponding qf,qfdot,T values
 * and average them
 *
 * INPUT:
 *
 * 	ball_state:     ball position and velocity estimates
 * 	params    :     optimization parameters to be loaded
 *
 *
 */
inline void knn(const mat & lookupt, const vec6 & ballstate, vec::fixed<15> & params, int k) {

	// find the closest entry
	static bool firsttime = true;
	static vec dots;
	static mat A;

	if (firsttime) {
		firsttime = false;
		A = lookupt.cols(span(X,DZ));
		for (unsigned i = 0; i < A.n_rows; i++) {
			dots(i) = dot(A.row(i), A.row(i));
		}
	}

	uvec idx = sort_index(dots - 2*A.t()*ballstate, "descend");

	vec fullvec = zeros<vec>(OPTIM_DIM + 2*NCART);
	for (int i = 0; i < k; i++) {
		fullvec += lookupt.row(idx(i));
	}
	params = fullvec(span(DZ+1,DZ+OPTIM_DIM)) / k;
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
 * Lookup a random entry with optimization coparameters (b0,v0) and optimization
 * main parameters (qf,qfdot,T)
 */
inline void lookup_random_entry(vec & coparams, vec & params) {

	mat lookup;
	load_lookup_table(lookup);
	int entry = as_scalar(randi<vec>(1,distr_param(0,LOOKUP_TABLE_SIZE-1)));
	vec lookup_state = lookup.row(entry).t();
	coparams = lookup_state(span(X,DZ));
	params = lookup_state(span(6,LOOKUP_COLUMN_SIZE-1));

}

BOOST_AUTO_TEST_CASE(test_nlopt_optim) {

	cout << "Testing NLOPT Optimization" << endl;
	double *q0dot = (double*)calloc(NDOF,sizeof(double));
	double *q0 = (double*)calloc(NDOF,sizeof(double));
	// initial guess for optim //
	double *lb = (double*)calloc(OPTIM_DIM,sizeof(double));
	double *ub = (double*)calloc(OPTIM_DIM,sizeof(double));
	double SLACK = 0.01;
	double Tmax = 1.0;

	cout << "Looking up a random entry..." << endl;
	vec::fixed<15> strike_params;
	vec6 ball_state;
	lookup_random_entry(ball_state,strike_params);

	vec7 init_joint_state;
	init_right_posture(init_joint_state);

	// initialize ball and racket //
	// predict for T_pred seconds
	set_bounds(lb,ub,SLACK,Tmax);

	double time2return = 1.0;
	racket* racket_params = send_racket_strategy(init_joint_state,ball_state,Tmax);
	vec3 normal_example;
	int example = 5;
	for (int i = 0; i < NCART; i++) {
		normal_example(i) = racket_params->normal[i][example];
	}
	BOOST_TEST(arma::norm(normal_example) == 1);

	init_coptim_params(init_joint_state, q0);
	coptim coparams = {q0, q0dot, q0, lb, ub, time2return};
	optim opt_params = {q0, q0dot, 0.5};

	// run NLOPT opt algorithm here //
	double max_violation = nlopt_optim_poly_run(&coparams,racket_params,&opt_params);

	// test to see if kinematics constraints are violated
	//double max_violation = test_optim(x,FALSE);
	BOOST_TEST(max_violation < 0.01);
}
