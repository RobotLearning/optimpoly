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

algo algs[] = {FIXED};
bool mpcs[] = {true,false};

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
 * Testing whether the ball can be returned to the opponents court
 */
BOOST_DATA_TEST_CASE(test_land, data::make(mpcs), mpc) {

	std::cout << "Running MPC test..." << std::endl;
	double Tmax = 1.0, lb[OPTIM_DIM], ub[OPTIM_DIM];
	set_bounds(lb,ub,0.01,Tmax);
	vec7 lbvec(lb); vec7 ubvec(ub);
	TableTennis tt;
	int num_trials = 5;
	int num_lands = 0;
	//arma_rng::set_seed_random();
	arma_rng::set_seed(2);
	vec7 q0;
	double std_obs = 1e-4; // std of the noisy observations
	init_right_posture(q0);
	joint qact = {q0, zeros<vec>(7), zeros<vec>(7)};
	vec3 obs;
	EKF filter = init_filter(0.0,std_obs*std_obs);
	Player* robot = new Player(q0,filter,FIXED,mpc,1);
	int N = 2000;
	joint qdes = {q0, zeros<vec>(NDOF), zeros<vec>(NDOF)};
	racket robot_racket;

	for (int n = 0; n < num_trials; n++) { // for each trial
		tt = TableTennis(false,true);
		std::cout << "New ball coming!" << std::endl;
		tt.set_ball_state(0.2);
		tt.load_params("test/ball_params_mismatch");
		robot->reset_filter(0.0,std_obs*std_obs);
		for (int i = 0; i < N; i++) { // one trial
			obs = tt.get_ball_position() + std_obs * randn<vec>(3);
			robot->play(qact, obs, qdes);
			calc_racket_state(qdes,robot_racket);
			tt.integrate_ball_state(robot_racket,DT);
			qact.q = qdes.q;
			qact.qd = qdes.qd;
		}
		if (tt.has_landed()) {
			num_lands++;
		}
	}
	std::cout << "======================================================" << endl;
	std::cout << "Out of " << num_trials << " attemps, "
			<< num_lands << " lands!" << std::endl;
	std::cout << "======================================================" << endl;
	delete(robot);
}

///*
// * Testing whether the ball can be returned to the opponents court
// */
//BOOST_DATA_TEST_CASE(test_land, data::make(algs), alg) {
//
//	std::cout << "*************** Testing Robot Ball Landing ************" << std::endl;
//
//	double Tmax = 1.0, lb[OPTIM_DIM], ub[OPTIM_DIM];
//	set_bounds(lb,ub,0.01,Tmax);
//	vec7 lbvec(lb); vec7 ubvec(ub);
//	TableTennis tt = TableTennis(false,true);
//	//arma_rng::set_seed_random();
//	arma_rng::set_seed(5);
//	tt.set_ball_state(0.2);
//
//	vec7 q0;
//	double std_obs = 0.0001; // std of the noisy observations
//	init_right_posture(q0);
//	joint qact = {q0, zeros<vec>(7), zeros<vec>(7)};
//	vec3 obs;
//	EKF filter = init_filter(0.0,std_obs*std_obs);
//	Player *robot = new Player(q0,filter,alg,true,2);
//
//	int N = 2000;
//	joint qdes = {q0, zeros<vec>(NDOF), zeros<vec>(NDOF)};
//	racket robot_racket;
//	mat Qdes = zeros<mat>(NDOF,N);
//	for (int i = 0; i < N; i++) {
//		//if (i % 20 == 0)
//		obs = tt.get_ball_position() + std_obs * randn<vec>(3);
//		robot->play(qact, obs, qdes);
//		//robot->cheat(qact, tt.get_ball_state(), qdes);
//		Qdes.col(i) = qdes.q;
//		calc_racket_state(qdes,robot_racket);
//		//cout << "robot ball dist\t" << norm(robot_racket.pos - tt.get_ball_position()) << endl;
//		tt.integrate_ball_state(robot_racket,DT);
//		usleep(DT*1e6);
//	}
//	std::cout << "Testing joint limits as well...\n";
//	BOOST_TEST(all(max(Qdes,1) < ubvec));
//	BOOST_TEST(all(min(Qdes,1) > lbvec));
//	BOOST_TEST(tt.has_landed());
//	delete(robot);
//	std::cout << "******************************************************" << std::endl;
//}

///*
// * Testing whether the ball can be returned to the opponents court
// */
//BOOST_DATA_TEST_CASE(test_mismatch, data::make(algs), alg) {
//
//	std::cout << "*************** Testing Prediction Mismatch *****************" << std::endl;
//
//	double Tmax = 1.0, lb[OPTIM_DIM], ub[OPTIM_DIM];
//	set_bounds(lb,ub,0.01,Tmax);
//	vec7 lbvec(lb); vec7 ubvec(ub);
//	TableTennis tt = TableTennis(false,true);
//	tt.load_params("test/ball_params_mismatch");
//	//arma_rng::set_seed_random();
//	arma_rng::set_seed(2);
//	tt.set_ball_state(0.2);
//
//	vec7 q0;
//	init_right_posture(q0);
//	joint qact = {q0, zeros<vec>(7), zeros<vec>(7)};
//	EKF filter = init_filter();
//	Player *robot = new Player(q0,filter,alg,true,2);
//
//	int N = 2000;
//	joint qdes = {q0, zeros<vec>(NDOF), zeros<vec>(NDOF)};
//	racket robot_racket;
//	for (int i = 0; i < N; i++) {
//		//robot->play(qact, tt.get_ball_position(), qdes);
//		robot->cheat(qact, tt.get_ball_state(), qdes);
//		calc_racket_state(qdes,robot_racket);
//		//cout << "robot ball dist\t" << norm(robot_racket.pos - tt.get_ball_position()) << endl;
//		//tt.integrate_ball_state(dt);
//		tt.integrate_ball_state(robot_racket,DT);
//		usleep(DT*1e6);
//	}
//	BOOST_TEST(tt.has_landed());
//	delete(robot);
//	std::cout << "******************************************************" << std::endl;
//}

/*
 * Testing whether table tennis ball bounces on table and touches the ground
 */
BOOST_AUTO_TEST_CASE( test_touch_ground ) {

	std::cout << std::endl << "Testing table tennis ball hitting ground..." << std::endl;
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
 * Testing whether the errors in the EKF filter estimate
 * are shrinking
 *
 */
BOOST_AUTO_TEST_CASE( test_ball_ekf ) {

	std::cout << std::endl <<
			"Testing EKF table tennis estimator with initial estimate error..." << std::endl;
	// initialize TableTennis and Filter classes
	TableTennis tt = TableTennis();
	double std = 0.000000001;
	mat C = eye<mat>(3,6);
	mat66 Q = zeros<mat>(6,6);
	mat33 R = std * eye<mat>(3,3);
	EKF filter = EKF(calc_next_ball,C,Q,R);

	// set table tennis ball and filter
	tt.set_ball_state(0.2);
	vec3 init_pos = tt.get_ball_position() + 0.5 * randu<vec>(3);
	vec3 init_vel = tt.get_ball_velocity() + 0.2 * randu<vec>(3);
	mat66 P0;
	P0.eye(6,6); P0 *= 1e6;
	vec6 x0 = join_vert(init_pos,init_vel);
	filter.set_prior(x0,P0);
	int N = 20;
	double dt = 0.01;
	vec6 ball_state;
	vec3 obs;
	vec err = zeros<vec>(N);

	for (int i = 0; i < N; i++) {
		tt.integrate_ball_state(dt);
		ball_state = tt.get_ball_state();
		obs = tt.get_ball_position() + 0.0 * randn<vec>(3);
		filter.predict(dt);
		filter.update(obs);
		err(i) = norm(filter.get_mean() - ball_state,2);
	}
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

	int N = 50;
	vec6 ball_state;
	vec3 obs;
	vec err = zeros<vec>(N);
	double std = 1e-8;
	TableTennis tt = TableTennis();
	EKF filter = init_filter(0.0,std);
	Player *cp = new Player(zeros<vec>(NDOF),filter);
	double std_obs = 0.001; // std of the noisy observations
	tt.set_ball_state(0.2);

	for (int i = 0; i < N; i++) {
		tt.integrate_ball_state(DT);
		ball_state = join_vert(tt.get_ball_position(),tt.get_ball_velocity());
		obs = tt.get_ball_position() + 0.001 * randn<vec>(3); // TODO: add noise!
		err(i) = norm(ball_state - cp->filt_ball_state(obs),2);
		//usleep(10e3);
	}
	cout << "Error of state estimate start: " << err(0) << " end: " << err(N-1) << endl;
	BOOST_TEST(err(N-1) < err(0), boost::test_tools::tolerance(0.01));
	delete(cp);
	//std::cout << " decreases." << std::endl;
	//BOOST_TEST(filter_est[Z] == floor_level, boost::test_tools::tolerance(0.1));
}
