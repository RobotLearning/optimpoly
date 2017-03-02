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

#include <boost/test/unit_test.hpp>
#include <armadillo>
#include "kalman.h"
#include "tabletennis.h"

using namespace std;
using namespace arma;

/*
 * Test for successful initialization
 * Check if uninitialized filter sends an exception
 */
inline bool test_uninit_exception(KF filter) {

	bool flag = false;
	try {
		vec mean = filter.get_mean();
	}
	catch (const char * exception) {
		cout << "Caught the exception!" << endl;
		flag = true;
	}

	return flag;
}

/*
 * Testing whether we can initialize the KF class
 */
BOOST_AUTO_TEST_CASE( test_init ) {

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

	BOOST_TEST(test_uninit_exception(filter));
}

/*
 * Testing the discretize funtion that
 * discretizes two continous model matrices A and B
 *
 */
BOOST_AUTO_TEST_CASE( test_discretize ) {

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

	BOOST_TEST(approx_equal(filter.get_model(1),A_matlab,"absdiff",1e-3));
	BOOST_TEST(approx_equal(filter.get_model(2),B_matlab,"absdiff",1e-3));

}

/*
 * Testing whether the random matrices are the same as in MATLAB
 */
BOOST_AUTO_TEST_CASE( test_random ) {

	int dimx = 2;
	int seed = 5;
	arma_rng::set_seed(seed);
	mat A(dimx,dimx,fill::randu);

	mat A_matlab;
	A_matlab << 0.6731 << 0.2253 << endr
			 << 0.0385 << 0.6759 << endr;

	//cout << A_matlab << endl;
	//cout << A << endl;

	BOOST_TEST(approx_equal(A,A_matlab,"absdiff",1e-3));
}

/*
 * Testing for a simple
 * predict/update performance
 */
BOOST_AUTO_TEST_CASE( test_update_predict ) {

	// init Kalman Filter
	int dimx = 2;
	int dimu = 1;
	int dimy = 1;
	// variance
	double eps = 0.1;
	int N = 10;

	int seed = 1;
	arma_rng::set_seed(seed);

	// model matrices
	mat A(dimx,dimx,fill::randu);
	mat B(dimx,dimu,fill::zeros);
	mat C(dimy,dimx,fill::zeros);
	C(0,1) = 1;
	// diagonal covariances
	mat Q(dimx,dimx,fill::zeros);
	mat R(dimy,dimy,fill::eye);
	R *= eps;

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

	//cout << answer << endl;
	//cout << filter.get_mean() << endl;

	BOOST_TEST(approx_equal(answer,filter.get_mean(),"absdiff",1e-3));
}

/*
 * Test the Extended Kalman Filter
 *
 * Using the table tennis prediction function for EKF
 */
BOOST_AUTO_TEST_CASE( test_ekf ) {

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
	//cout << "x0" << endl << filter.get_mean() << endl;
	//cout << "P0" << endl << filter.get_covar() << endl;
	filter.predict(0.01);
	//cout << "x1" << endl << filter.get_mean() << endl;
	//cout << "P1" << endl << filter.get_covar() << endl;

	//BOOST_TEST(true);
}

/*
 * Test predict path function of EKF with table tennis
 */
BOOST_AUTO_TEST_CASE( test_predict_path ) {

	TableTennis tt = TableTennis();
	arma_rng::set_seed(5);

	// initialize a filter to predict
	mat C(3,6,fill::zeros);
	C.cols(1,3) = eye<mat>(3,3);
	// diagonal covariances
	mat Q(6,6,fill::zeros);
	mat R(3,3,fill::eye);
	R *= 0.001;
	EKF filter = EKF(calc_next_ball,C,Q,R);

	int N = 50;
	double dt = 0.01;
	tt.set_ball_state(0.2);

	vec6 ball_state = join_vert(tt.get_ball_position(),tt.get_ball_velocity());
	filter.set_prior(ball_state,0.01 * eye<mat>(6,6));

	mat M = zeros<mat>(6,N);
	for (int i = 0; i < N; i++) {
		tt.integrate_ball_state(dt);
		M.col(i) = join_vert(tt.get_ball_position(),tt.get_ball_velocity());
	}
	mat Mfilt = filter.predict_path(dt,N);

	BOOST_TEST(approx_equal(M, Mfilt, "absdiff", 0.002));

}
