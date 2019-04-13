/*
 * Test cases for the Kalman Filter C++ implementation.

 * Kalman Filtering unit test suite
 * for testing various use cases.
 *
 * We can compare with MATLAB classes (on another repo)
 * by using ARMA generated random matrices with same seeds on both sides!
 *
 *  Created on: Jan 25, 2017
 *      Author: okoc
 */

#include "gtest/gtest.h"
#include <armadillo>
#include "kalman.h"
#include "player.hpp"
#include "tabletennis.h"

using namespace std;
using namespace arma;
using namespace const_tt;
using namespace player;

static bool test_uninit_exception(const KF & filter);

TEST(KalmanFilterTests, TestKalmanFilterInitialization) {

	// init Kalman Filter
	int dimx = 2;
	int dimu = 1;
	int dimy = 1;
	// model matrices
	mat A(dimx,dimx,fill::randu);
	mat B(dimx,dimu,fill::randu);
	mat C(dimy,dimx,fill::zeros);
	C(0,0) = 1.0;
	// diagonal covariances
	mat Q(dimx,dimx,fill::eye);
	mat R(dimy,dimy,fill::eye);
	double s2x = 0.01;
	double s2y = 0.10;
	Q *= sqrt(s2x);
	R *= sqrt(s2y);
	KF filter = KF(A,B,C,Q,R);
	EXPECT_TRUE(test_uninit_exception(filter));
}

/*
 * Testing the discretize funtion that
 * discretizes two continous model matrices A and B
 *
 */
TEST(KalmanFilterTests, CompareContinuousModelDiscretizationWithMATLAB) {

	int dimx = 2;
	int dimy = 1;
	int dimu = 1;
	int seed = 1;
	arma_rng::set_seed(seed);
	mat C(dimy,dimx,fill::zeros);
	C(0,0) = 1.0;
	// diagonal covariances
	mat Q(dimx,dimx,fill::eye);
	mat R(dimy,dimy,fill::eye);
	KF filter = KF(C, Q, R);

	mat A_matlab;
	A_matlab << 1.0138 << 0.0455 << endr
	         << 0.0137 << 1.0024 << endr;
	mat B_matlab;
	B_matlab << 0.0374 << endr
			 << 0.0915 << endr;

	mat Ac(dimx,dimx,fill::randu);
	mat Bc(dimx,dimu,fill::randu);
	double dt = 0.1;
	filter.discretize(Ac,Bc,dt);

	EXPECT_TRUE(approx_equal(filter.get_model(1),A_matlab,"absdiff",1e-3));
	EXPECT_TRUE(approx_equal(filter.get_model(2),B_matlab,"absdiff",1e-3));

}

TEST(RandomizationTests, CompareRandomMatrixWithMATLABValues) {

	int dimx = 2;
	int seed = 5;
	arma_rng::set_seed(seed);
	mat A(dimx,dimx,fill::randu);

	mat A_matlab;
	A_matlab << 0.6731 << 0.2253 << endr
			 << 0.0385 << 0.6759 << endr;

	EXPECT_TRUE(approx_equal(A,A_matlab,"absdiff",1e-3));
}

/*
 * Testing for a simple
 * predict/update performance
 */
TEST(KalmanFilterTests, CompareStateAfterPredictUpdateToMATLABValues) {

	// init Kalman Filter
	int dimx = 2;
	int dimu = 1;
	int dimy = 1;
	// variance
	double eps = 0.1;
	int N = 10;
	arma_rng::set_seed(1);

	// model matrices
	mat A(dimx,dimx,fill::randu);
	mat B(dimx,dimu,fill::zeros);
	mat C(dimy,dimx,fill::zeros);
	C(0,1) = 1;
	// diagonal covariances
	mat Q(dimx,dimx,fill::eye);
	mat R(dimy,dimy,fill::eye);
	R *= eps;
	Q *= 1e-10;

	vec x0;
	x0 << 10 << 1;
	mat P0(dimx,dimx,fill::eye);
	P0 *= eps;
	KF filter = KF(x0,P0,A,B,C,Q,R);

	mat Ys = filter.sample_observations(N);

	for (int i = 0; i < N-1; i++) {
		filter.update(Ys.col(i));
		filter.predict();
	}
	filter.update(Ys.col(N-1));

	vec answer;
	answer << 1e-3 * 0.3441 << 1e-3 * 0.1516;

	EXPECT_TRUE(approx_equal(answer,filter.get_mean(),"absdiff",1e-3));
}

TEST(KalmanFilterTests, RunExtendedKalmanFilterInitializationAndPrediction) {

	int dimx = 6;
	int dimy = 3;
	double eps = 0.01;

	mat C(dimy,dimx,fill::zeros);
	C.cols(1,dimy) = eye<mat>(dimy,dimy);
	// diagonal covariances
	mat Q(dimx,dimx,fill::zeros);
	mat R(dimy,dimy,fill::eye);
	R *= eps;
	EKF filter = EKF(calc_next_ball,C,Q,R);

	vec x0 = zeros<vec>(dimx);
	mat P0(dimx,dimx,fill::eye);
	P0 *= eps;
	filter.set_prior(x0,P0);
	filter.predict(0.01);
}

/*
 * Test predict path function of EKF with table tennis
 *
 */
TEST(KalmanFilterTests, ComparePathPredictionToActualValuesForTableTennisBall) {

	TableTennis tt = TableTennis();
	arma_rng::set_seed(5);

	// initialize a filter to predict
	EKF filter = init_ball_filter();

	int N = 1000;
	double dt = 0.002;
	tt.set_ball_gun(0.2);

	vec6 ball_state = tt.get_ball_state();
	filter.set_prior(ball_state,0.01 * eye<mat>(6,6));

	mat M = zeros<mat>(6,N);
	for (int i = 0; i < N; i++) {
		tt.integrate_ball_state(dt);
		M.col(i) = join_vert(tt.get_ball_position(),tt.get_ball_velocity());
	}
	mat Mfilt = filter.predict_path(dt,N);

	EXPECT_TRUE(approx_equal(M, Mfilt, "absdiff", 0.002));
}


static bool test_uninit_exception(const KF & filter) {

    bool flag = false;
    try {
        vec mean = filter.get_mean();
    }
    catch (const runtime_error & exception) {
        cout << "Caught the exception!" << endl;
        flag = true;
    }

    return flag;
}
