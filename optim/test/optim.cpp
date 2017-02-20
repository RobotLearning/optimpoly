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
#include "SL.h"
#include "constants.h"
#include "kinematics.h"
#include "optimpoly.h"
#include "utils.h"
#include "stdlib.h"
#include "player.hpp"

using namespace std;
using namespace arma;

/*
 * Set upper and lower bounds on the optimization.
 * First loads the joint limits and then
 */
inline void set_bounds(double *lb,
		               double *ub, double SLACK, double Tmax) {

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
 * Set the ball values to a reasonable value
 * SO FAR setting it to the first lookup table entry from May 2016
 *
 */
inline void init_ball_state(double *b0, double *v0) {

	// initialize the ball
	b0[0] = 0.1972;
	b0[1] = -2.4895;
	b0[2] = -0.5040;
	v0[0] = -1.7689;
	v0[1] = 4.7246;
	v0[2] = -1.0867;
}

/*
 * Set the initial posture of the robot
 */
inline void init_joint_state(double *q0) {

	// initialize the variables
	//q0 = [1.0; -0.2; -0.1; 1.8; -1.57; 0.1; 0.3];
	q0[0] = 1.0;
	q0[1] = -0.2;
	q0[2] = -0.1;
	q0[3] = 1.8;
	q0[4] = -1.57;
	q0[5] = 0.1;
	q0[6] = 0.3;
}

/*
 * Create the desired racket structure for the optimization
 * to be performed in C
 */
inline cracket make_c_strategy(const racket & strategy) {

	int N = strategy.pos.n_cols;
	Matrix pos = my_matrix(0,NCART,0,N);
	Matrix vel = my_matrix(0,NCART,0,N);
	Matrix normal = my_matrix(0,NCART,0,N);

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < NCART; j++) {
			pos[j][i] = strategy.pos(j,i);
			vel[j][i] = strategy.vel(j,i);
			normal[j][i] = strategy.normal(j,i);
		}
	}
	cracket params = {pos, vel, normal, N};
	return params;
}

/*
 * Find the closest ball state in the lookup table
 * and lookup the corresponding qf,qfdot,T values
 *
 * INPUT:
 * 	b0,v0: ball position and velocity estimates
 * 	x:     optimization parameters to be loaded
 *
 * TODO: expand with more sophisticated interpolation methods (kNN with k = 1 at the moment)
 */
int lookup(const Matrix lookupTable, const double* b0, const double* v0, double *x) {

	int i,j;
	double minVal = 0.0;
	double bestVal = 1e6;
	int bestIdx = 1;
	Vector ballVec = my_vector(1,2*NCART);

	for (i = 1; i <= NCART; i++) {
		ballVec[i] = b0[i-1];
		ballVec[i+NCART] = v0[i-1];
	}

	for (i = 1; i <= LOOKUP_TABLE_SIZE; i++) {
		minVal = 0.0;
		for (j = 1; j <= 2*NCART; j++) {
			minVal += pow(ballVec[j] - lookupTable[i][j],2);
		}
		if (minVal < bestVal) {
			bestIdx = i;
			bestVal = minVal;
		}
	}

	for (i = 1; i <= NDOF; i++) {
		x[i-1] = lookupTable[bestIdx][2*NCART+i];
		x[i-1+NDOF] = lookupTable[bestIdx][2*NCART+NDOF+i];
	}
	x[2*NDOF] = lookupTable[bestIdx][2*NCART+2*NDOF+1];

	return TRUE;
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
	double *b0 = (double*)calloc(NCART,sizeof(double));
	double *v0 = (double*)calloc(NCART,sizeof(double));
	double *x = (double*)calloc(OPTIM_DIM,sizeof(double));
	// initial guess for optim //
	double *lb = (double*)calloc(OPTIM_DIM,sizeof(double));
	double *ub = (double*)calloc(OPTIM_DIM,sizeof(double));
	double SLACK = 0.01;
	double Tmax = 1.0;

	Matrix lookup_table = my_matrix(1, LOOKUP_TABLE_SIZE, 1, LOOKUP_COLUMN_SIZE);
	load_lookup_table(lookup_table);
	//set_des_land_param(ballLand,&landTime);

	// initialize ball and racket //
	// predict for T_pred seconds
	set_bounds(lb,ub,SLACK,Tmax);
	init_joint_state(q0);
	init_ball_state(b0,v0);
	lookup(lookup_table,b0,v0,x); // initialize solution
	double time2return = 1.0;
	racket strategy = send_racket_strategy(q0,b0,v0,Tmax);
	cracket racket = make_c_strategy(strategy);

	/*BOOST_TEST(racket.pos[Z][5] != strategy.pos(Z,5));
	BOOST_TEST(racket.normal[X][5] != strategy.normal(X,5));
	BOOST_TEST(racket.vel[Y][5] != strategy.vel(Y,5));*/

	/*cout << strategy.pos << endl
		 << strategy.vel << endl
		 << strategy.normal	<< endl;*/


	BOOST_TEST(arma::norm(strategy.normal.col(5)) == 1);

	coptim params = {q0, q0dot, q0, lb, ub, time2return};

	// run NLOPT opt algorithm here //
	double max_violation = nlopt_optim_poly_run(&params,&racket);

	// test to see if kinematics constraints are violated
	//double max_violation = test_optim(x,FALSE);
	BOOST_TEST(max_violation < 0.01);
}
