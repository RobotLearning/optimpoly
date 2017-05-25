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

#ifndef BOOST_TEST_MODULE
#define BOOST_TEST_MODULE test_table_tennis
#endif

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/data/monomorphic.hpp>
#include <armadillo>
#include "player.hpp"
#include "constants.h"
#include "tabletennis.h"
#include "kinematics.hpp"
#include "kalman.h"

using namespace arma;
namespace data = boost::unit_test::data;

/*
 * Initialize robot posture
 */
inline void init_posture(vec7 & q0, int posture) {

	rowvec qinit;
	switch (posture) {
	case 0: // right
		qinit << 1.0 << -0.2 << -0.1 << 1.8 << -1.57 << 0.1 << 0.3 << endr;
		break;
	case 1: // center
		qinit << 0.0 << 0.0 << 0.0 << 1.5 << -1.75 << 0.0 << 0.0 << endr;
		break;
	case 2: // left
		qinit << -1.0 << 0.0 << 0.0 << 1.5 << -1.57 << 0.1 << 0.3 << endr;
		break;
	default: // default is the right side
		qinit << 1.0 << -0.2 << -0.1 << 1.8 << -1.57 << 0.1 << 0.3 << endr;
		break;
	}
	q0 = qinit.t();
}

algo algs[] = {FIXED, VHP, LAZY};

/*
 * Testing whether the ball can be returned to the opponents court
 *
 * Comment from Raffi:
 *

    so if you know what those numbers should be for each algorithm

    you can add their expected values in the algs as parameters of a std::tuple
    (or a combination of boost.test datasets with the join operator
    http://www.boost.org/doc/libs/1_64_0/libs/test/doc/html/boost_test/tests_organization/
    test_cases/test_case_generation/operations.html
    #boost_test.tests_organization.test_cases.test_case_generation.operations.joins).
    and add a check on the excepted value inside the unit test

    It would mean that the seed is fixed.
 *
 *
 */
/*
BOOST_DATA_TEST_CASE(test_land_mpc, data::make(algs), alg) {

	std::cout << "Running MPC Test..." << std::endl;
	double Tmax = 1.0, lb[2*NDOF+1], ub[2*NDOF+1];
	set_bounds(lb,ub,0.01,Tmax);
	vec7 lbvec(lb); vec7 ubvec(ub);
	TableTennis tt;
	int num_trials = 1;
	int num_lands = 0;
	int num_misses = 0;
	int num_not_valid = 0;
	//arma_rng::set_seed_random();
	arma_rng::set_seed(0);
	vec7 q0;
	double std_noise = 0.001;
	double std_model = 0.3;
	init_posture(q0,2);
	joint qact = {q0, zeros<vec>(7), zeros<vec>(7)};
	vec3 obs;
	EKF filter = init_filter(std_model,std_noise);
	Player robot = Player(q0,filter,alg,true,0);
	int N = 2000;
	joint qdes = qact;
	racket robot_racket;

	for (int n = 0; n < num_trials; n++) { // for each trial
		tt = TableTennis(true,true);
		std::cout << "Trial: " << n+1 << std::endl;
		tt.set_ball_gun(0.05,2);
		robot.reset_filter(std_model,std_noise);
		for (int i = 0; i < N; i++) { // one trial
			obs = tt.get_ball_position() + std_noise * randn<vec>(3);
			robot.play(qact, obs, qdes);
			calc_racket_state(qdes,robot_racket);
			tt.integrate_ball_state(robot_racket,DT);
			qact.q = qdes.q;
			qact.qd = qdes.qd;
		}
		if (tt.has_landed()) {
			num_lands++;
		}
		else if (!tt.is_legal_ball())
			num_not_valid++;
		else
			num_misses++;
	}
	std::cout << "======================================================" << endl;
	std::cout << "Out of " << num_trials << " trials, "
			<< num_lands << " lands, " << num_not_valid <<
			" not valid balls, " << num_misses << " misses!" <<std::endl;
	std::cout << "======================================================" << endl;
}*/

/*
 * Testing whether the ball can be returned to the opponents court
 */
BOOST_DATA_TEST_CASE(test_land, data::make(algs), alg) {

	std::cout << "Testing Robot Ball Landing" << std::endl;

	double Tmax = 1.0, lb[2*NDOF+1], ub[2*NDOF+1];
	set_bounds(lb,ub,0.01,Tmax);
	vec7 lbvec(lb); vec7 ubvec(ub);
	TableTennis tt = TableTennis(false,true);
	//arma_rng::set_seed_random();
	arma_rng::set_seed(5);
	tt.set_ball_gun(0.05,1); // init ball on the centre

	vec7 q0;
	double std_obs = 0.000; // std of the noisy observations
	init_posture(q0,0);
	joint qact = {q0, zeros<vec>(7), zeros<vec>(7)};
	vec3 obs;
	EKF filter = init_filter(0.03,std_obs);
	Player robot = Player(q0,filter,alg,false,2);

	int N = 2000;
	joint qdes = {q0, zeros<vec>(NDOF), zeros<vec>(NDOF)};
	racket robot_racket;
	mat Qdes = zeros<mat>(NDOF,N);
	for (int i = 0; i < N; i++) {
		obs = tt.get_ball_position() + std_obs * randn<vec>(3);
		robot.play(qact, obs, qdes);
		//robot.cheat(qact, tt.get_ball_state(), qdes);
		Qdes.col(i) = qdes.q;
		calc_racket_state(qdes,robot_racket);
		//cout << "robot ball dist\t" << norm(robot_racket.pos - tt.get_ball_position()) << endl;
		tt.integrate_ball_state(robot_racket,DT);
		//usleep(DT*1e6);
		qact.q = qdes.q;
		qact.qd = qdes.qd;
	}
	//cout << max(Qdes,1) << endl;
	std::cout << "Testing joint limits as well...\n";
	BOOST_TEST(all(max(Qdes,1) < ubvec));
	BOOST_TEST(all(min(Qdes,1) > lbvec));
	BOOST_TEST(tt.has_landed());
	std::cout << "******************************************************" << std::endl;
}

/*
 * Testing whether table tennis ball bounces on table and touches the ground
 */
BOOST_AUTO_TEST_CASE( test_touch_ground ) {

	std::cout << std::endl << "Testing table tennis ball hitting ground..." << std::endl;
	TableTennis tt = TableTennis();

	int N = 200;
	double dt = 0.01;
	tt.set_ball_gun(0.2);
	for (int i = 0; i < N; i++) {
		tt.integrate_ball_state(dt);
	}

	vec3 ball_pos = tt.get_ball_position();
	BOOST_TEST(ball_pos(Z) == floor_level, boost::test_tools::tolerance(0.01));
}

/*
 * Testing whether the errors in the EKF filter estimate
 * are shrinking
 *
 */
BOOST_AUTO_TEST_CASE( test_ball_ekf ) {

	std::cout << std::endl <<
			"Testing EKF table tennis estimator with initial estimate error..." << std::endl;
	// initialize TableTennis and Filter classes
	TableTennis tt = TableTennis(false,true);
	const double std_noise = 0.001;
	const double std_model = 0.03;
	mat C = eye<mat>(3,6);
	mat66 Q = std_model * eye<mat>(6,6);
	mat33 R = std_noise * eye<mat>(3,3);
	EKF filter = EKF(calc_next_ball,C,Q,R);

	// set table tennis ball and filter
	tt.set_ball_gun(0.2);
	vec3 init_pos = tt.get_ball_position() + 0.5 * randu<vec>(3);
	vec3 init_vel = tt.get_ball_velocity() + 0.2 * randu<vec>(3);
	mat66 P0;
	P0.eye(6,6);
	P0 *= 1e6;
	vec6 x0 = join_vert(init_pos,init_vel);
	filter.set_prior(x0,P0);
	const int N = 100;
	const double dt = 0.01;
	vec3 obs;
	vec err = zeros<vec>(N);

	for (int i = 0; i < N; i++) {
		tt.integrate_ball_state(dt);
		obs = tt.get_ball_position() + std_noise * randn<vec>(3);
		filter.predict(dt);
		filter.update(obs);
		err(i) = norm(filter.get_mean() - tt.get_ball_state(),2);
	}
	//cout << err << endl;
	cout << "Error of state estimate start: " << err(0) << " end: " << err(N-1) << endl;
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

	std::cout << std::endl << "Testing Player class's Filtering performance..." << std::endl;

	const int N = 50;
	vec3 obs;
	vec err = zeros<vec>(N);
	const double std_noise = 0.001;
	const double std_model = 0.001;
	TableTennis tt = TableTennis(false,true);
	EKF filter = init_filter(std_model,std_noise);
	Player cp = Player(zeros<vec>(NDOF),filter);
	tt.set_ball_gun(0.2);

	for (int i = 0; i < N; i++) {
		tt.integrate_ball_state(DT);
		obs = tt.get_ball_position() + std_noise * randn<vec>(3);
		err(i) = norm(tt.get_ball_state() - cp.filt_ball_state(obs),2);
		//usleep(10e3);
	}
	//cout << err << endl;
	cout << "Error of state estimate start: " << err(0) << " end: " << err(N-1) << endl;
	BOOST_TEST(err(N-1) < err(0), boost::test_tools::tolerance(0.01));
	//cout << err;
	//BOOST_TEST(filter_est[Z] == floor_level, boost::test_tools::tolerance(0.1));
}
