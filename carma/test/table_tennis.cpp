/*
 * Unit tests for table tennis related functionalities
 * residing in src/table_tennis.cpp
 *
 * The aims are
 *
 * 1. Check if ball functions in Table Tennis class are working well:
 *    - check if ball touches ground, hits table, etc.
 * 2. Check if the Table Tennis Player works as it should:
 *    - ball estimation
 *    - ball prediction
 *    - trajectory generation
 *
 * table_tennis.cpp
 *
 *  Created on: Feb 2, 2017
 *      Author: okoc
 */

#define BOOST_TEST_MODULE test_table_tennis
#include <boost/test/unit_test.hpp>
#include <armadillo>
#include "player.hpp"
#include "tabletennis.h"
#include "kalman.h"

using namespace std;
using namespace arma;


/*
 * TODO:
 */
BOOST_AUTO_TEST_CASE(test_player) {

	cout << "Testing Robot racket calculations..." << endl;
	static double pos[NCART] = {1.0, -2.0, -0.5};
	static double vel[NCART] = {3.0, 5.0, 4.0};
	EKF filter = init_filter();
	vec3 ballpos(pos);
	vec3 ballvel(vel);
	mat66 P; P.eye();
	filter.set_prior(join_vert(ballpos,ballvel),P);
	mat balls_pred = filter.predict_path(dt,10);
	//cout << "Balls predicted:" << endl << balls_pred << endl;

}


/*
 * Testing whether table tennis ball bounces on table and touches the ground
 */
BOOST_AUTO_TEST_CASE( test_touch_ground ) {

	TableTennis tt = TableTennis();

	int N = 200;
	double dt = 0.01;
	tt.set_ball_state(0.2);
	for (int i = 0; i < N; i++) {
		tt.integrate_ball_state(dt);
	}

	vec3 ball_pos = tt.get_ball_position();
	BOOST_TEST(ball_pos(Z) == floor_level, boost::test_tools::tolerance(0.01));
}

/*
 * Testing whether the errors in the KF filter estimate
 * are shrinking
 *
 */
/*BOOST_AUTO_TEST_CASE( test_ball_kf ) {

	cout << endl << "Running KF table tennis estimator" << endl;
	// initialize TableTennis and Filter classes
	TableTennis tt = TableTennis();
	double std = 0.001;
	mat C = eye<mat>(3,7);
	mat77 Q = zeros<mat>(7,7);
	mat33 R = std * eye<mat>(3,3);
	KF filter = KF(C,Q,R);
	mat77 Ac = zeros<mat>(7,7);  // continuous
	mat77 B = zeros<mat>(7,7); // not used
	// fill the continuous matrix
	Ac(span(X,Z),span(DX,DZ)) = eye<mat>(3,3);
	Ac(DZ,DZ+1) = gravity;

	// set table tennis ball and filter
	tt.set_ball_state(0.2);
	vec3 init_pos = tt.get_ball_position() + 0.0 * ones<vec>(3);
	vec3 init_vel = tt.get_ball_velocity();
	mat77 P0;
	P0.eye(7,7);
	vec7 x0 = join_vert(join_vert(init_pos,init_vel),vec(1,fill::ones));
	filter.set_prior(x0,P0);

	int N = 20;
	double dt = 0.001;
	vec3 ball_pos;
	vec3 ball_vel;
	vec err = zeros<vec>(N);

	for (int i = 0; i < N; i++) {
		tt.integrate_ball_state(dt);
		ball_pos = tt.get_ball_position();
		ball_vel = tt.get_ball_velocity();
		filter.discretize(Ac,B,dt);
		filter.predict();
		filter.update(ball_pos);
		err(i) = norm(filter.get_mean()(span(X,DZ)) - join_vert(ball_pos,ball_vel),2);
	}

	cout << "Error of state estimate" << endl << err.t() << endl;
	BOOST_TEST(err(N-1) < err(0));
}*/

/*
 * Testing whether the errors in the EKF filter estimate
 * are shrinking
 *
 */
BOOST_AUTO_TEST_CASE( test_ball_ekf ) {

	cout << endl << "Running EKF table tennis estimator..." << endl;
	// initialize TableTennis and Filter classes
	TableTennis tt = TableTennis();
	double std = 0.001;
	mat C = eye<mat>(3,6);
	mat66 Q = zeros<mat>(6,6);
	mat33 R = std * eye<mat>(3,3);
	EKF filter = EKF(calc_next_ball,C,Q,R);

	// set table tennis ball and filter
	tt.set_ball_state(0.2);
	vec3 init_pos = tt.get_ball_position() + 0.5 * randu<vec>(3);
	vec3 init_vel = tt.get_ball_velocity() + 0.2 * randu<vec>(3);
	mat66 P0;
	P0.eye(6,6);
	vec6 x0 = join_vert(init_pos,init_vel);
	filter.set_prior(x0,P0);

	int N = 20;
	double dt = 0.01;
	vec3 ball_pos;
	vec3 ball_vel;
	vec err = zeros<vec>(N);

	for (int i = 0; i < N; i++) {
		tt.integrate_ball_state(dt);
		ball_pos = tt.get_ball_position();
		ball_vel = tt.get_ball_velocity();
		filter.predict(dt);
		filter.update(ball_pos);
		err(i) = norm(filter.get_mean()(span(X,DZ)) - join_vert(ball_pos,ball_vel),2);
	}

	//cout << "Error of state estimate" << endl << err << endl;
	BOOST_TEST(err(N-1) <= err(0), boost::test_tools::tolerance(0.0001));
}


/*
 * Testing the EKF filter of Player class
 *
 * We expect the error to be decreasing
 * at some point
 *
 */
BOOST_AUTO_TEST_CASE( test_player_ekf_filter ) {

	cout << endl << "Testing Player class's Filtering performance" << endl;
	TableTennis tt = TableTennis();
	EKF filter = init_filter();
	vec3 init_pos = tt.get_ball_position() + 0.5 * randu<vec>(3);
	vec3 init_vel = tt.get_ball_velocity() + 0.2 * randu<vec>(3);
	mat66 P0;
	P0.eye(6,6);
	vec6 x0 = join_vert(init_pos,init_vel);
	filter.set_prior(x0,P0);
	vec7 q0 = zeros<vec>(NDOF);
	Player cp = Player(q0,filter);

	int N = 20;
	double dt = 0.01;
	tt.set_ball_state(0.2);
	vec6 ball_state;
	double racket_pos[3] = {0.0};
	double racket_orient[4] = {0.0};
	double filter_est[6] = {0.0};
	vec3 obs;
	vec err = zeros<vec>(N);
	int reset = 0;

	for (int i = 0; i < N; i++) {
		tt.integrate_ball_state(dt);
		ball_state = join_vert(tt.get_ball_position(),tt.get_ball_velocity());
		obs = ball_state(span(X,Z)); // TODO: add noise!
		err(i) = norm(ball_state - cp.filt_ball_state(obs),2);
		usleep(10e3);

	}
	cout << "Error of state estimate"; //<< endl << err.t() << endl;
	BOOST_TEST(err(N-1) < err(0), boost::test_tools::tolerance(0.1));
	cout << " decreases." << endl;
	//BOOST_TEST(filter_est[Z] == floor_level, boost::test_tools::tolerance(0.1));
}
