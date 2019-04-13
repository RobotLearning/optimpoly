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

#include "gtest/gtest.h"
#include <armadillo>
#include "player.hpp"
#include "constants.h"
#include "tabletennis.h"
#include "kinematics.hpp"
#include "kalman.h"

using namespace arma;
using namespace const_tt;
using namespace player;
using namespace optim;

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

TEST(TableTennisTests, TestBallTouchesGround) {
	TableTennis tt = TableTennis();
	int N = 200;
	double dt = 0.01;
	tt.set_ball_gun(0.2);
	for (int i = 0; i < N; i++) {
		tt.integrate_ball_state(dt);
	}
	vec3 ball_pos = tt.get_ball_position();
	EXPECT_EQ(ball_pos(Z), floor_level); // tol = 0.01
}

TEST(TableTennisTests, TestEKFBallEstimateImprovesOverTime) {

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
    std::cout << "Error of state estimate start: " << err(0) << " end: " << err(N-1);
    EXPECT_LE(err(N-1),err(0)); // tol = 0.0001
}


TEST(TableTennisTests, TestPlayersBallEstimateImprovesOverTime) {

	const int N = 50;
	vec3 obs;
	vec err = zeros<vec>(N);
	const double std_noise = 0.001;
	TableTennis tt = TableTennis(false,true);
	player_flags flags;
	Player cp = Player(zeros<vec>(NDOF),flags);
	tt.set_ball_gun(0.2);

	for (int i = 0; i < N; i++) {
		tt.integrate_ball_state(DT);
		obs = tt.get_ball_position() + std_noise * randn<vec>(3);
		err(i) = norm(tt.get_ball_state() - cp.filt_ball_state(obs),2);
		//usleep(10e3);
	}
	//cout << err << endl;
	cout << "Error of state estimate start: " << err(0) << " end: " << err(N-1) << endl;
	EXPECT_LT(err(N-1), err(0));
}
