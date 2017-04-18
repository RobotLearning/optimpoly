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

#ifndef BOOST_TEST_MODULE
#define BOOST_TEST_MODULE test_table_tennis
#endif

#include <boost/test/unit_test.hpp>
#include <armadillo>
#include "kalman.h"
#include "player.hpp"
#include "tabletennis.h"

using namespace std;
using namespace arma;

static bool fuse_blobs(const vec3 & blob1, const vec3 & blob3,
		               const bool & status1, const bool & status3, vec3 & obs);
static bool check_blob_validity(const vec3 & blob, const bool & status);

/*
 * Test for successful initialization
 * Check if uninitialized filter sends an exception.
 */
inline bool test_uninit_exception(const KF & filter) {

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
 *
 */
BOOST_AUTO_TEST_CASE( test_predict_path ) {

	TableTennis tt = TableTennis();
	arma_rng::set_seed(5);

	// initialize a filter to predict
	EKF filter = init_filter();

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

	BOOST_TEST(approx_equal(M, Mfilt, "absdiff", 0.002));

}

/*
 * Test mismatch in ball prediction
 * by using a spin model and a normal drag model for prediction
 *
 */
BOOST_AUTO_TEST_CASE( test_mismatch_pred ) {

	std::cout << "Testing ball prediction accuracy in mismatch case (spin)" << std::endl;
	TableTennis tt = TableTennis(true, false); //init with spin
	arma_rng::set_seed(5);
	EKF filter = init_filter(0.3,0.001);
	int N = 100;
	tt.set_ball_gun(0.2);
	vec6 ball_init = tt.get_ball_state();
	filter.set_prior(ball_init,eye<mat>(6,6));
	vec pos_errs, vel_errs;
	pos_errs = vel_errs = zeros<vec>(N);
	for (int i = 0; i < N; i++) {
		tt.integrate_ball_state(DT);
		filter.predict(DT,true);
		filter.update(tt.get_ball_position());
		pos_errs(i) = norm(filter.get_mean()(span(X,Z)) - tt.get_ball_position());
		vel_errs(i) = norm(filter.get_mean()(span(DX,DZ)) - tt.get_ball_velocity());
	}
	std::cout << "Pos error max: " << max(pos_errs) << std::endl;
	std::cout << "Vel error max: " << max(vel_errs) << std::endl;

}

/*
 * Test outlier prediction with Extended Kalman Filter
 *
 * Here we want to test/debug the Kalman Filter
 * that is to be used realtime with our vision system
 * for ball tracking
 *
 * Things to watch out for:
 * 1. It should reset correctly everytime ball is launched from ball gun
 * 2. Outliers should be cleaned properly
 * 3. New balls should be updated
 *
 */
BOOST_AUTO_TEST_CASE( test_outlier_detection ) {

	std::cout << "**************************************\n";
	std::cout << "Testing filtering on REAL BALL DATA!\n";
	bool status1, status3;
	//double time_data = 0.0;
	int head = 0; //3300;
	static vec3 blob1, blob3, obs;
	mat real_ball_data;
	std::string home = std::getenv("HOME");
	try {
		real_ball_data.load(home + "/Dropbox/data/realBallData_030516.txt");
	}
	catch (const char * exception) {
		std::cout << "Problem accessing/finding real ball data on Dropbox!" << std::endl;
	}
	int N = real_ball_data.n_rows; //4000;
	mat ball_states = zeros<mat>(N-head,6);
	EKF filter = init_filter(0.3,0.0001);
	Player cp = Player(zeros<vec>(7),filter,FIXED,false,1);
	for (int i = head; i < N; i++) {
		status1 = real_ball_data(i,1);
		blob1 = real_ball_data(i,span(2,4)).t();
		status3 = real_ball_data(i,6);
		blob3 = real_ball_data(i,span(7,9)).t();
		fuse_blobs(blob1,blob3,status1,status3,obs);
		//time_data = real_ball_data(i,10);
		ball_states.row(i-head) = cp.filt_ball_state(obs).t();
		usleep(2000);
	}
	ball_states.save(home + "/Dropbox/data/realBallData_filtered.txt",raw_ascii);
}

/*
 *
 * Fusing the blobs
 * If both blobs are valid blob3 is preferred
 * Only updates if the blobs are valid, i.e. not obvious outliers
 *
 */
static bool fuse_blobs(const vec3 & blob1, const vec3 & blob3,
		               const bool & status1, const bool & status3, vec3 & obs) {

	static bool status = false;

	// if ball is detected reliably
	// Here we hope to avoid outliers and prefer the blob3 over blob1
	if (check_blob_validity(blob3,status3) ||
			check_blob_validity(blob1,status1)) {
		status = true;
		obs = (status3) ? blob3 : blob1;
	}
	return status;
}

/*
 *
 * Checks for the validity of blob ball data using obvious table tennis checks.
 * Returns TRUE if valid.
 *
 * Does not use uncertainty estimates to assess validity
 * so do not rely on it as the sole source of outlier detection!
 *
 */
static bool check_blob_validity(const vec3 & blob, const bool & status) {

	bool valid;
	static double zMax = 0.5;
	static double zMin = floor_level - table_height;
	static double xMax = table_width/2.0;
	static double yMax = 0.5;
	static double yMin = dist_to_table - table_length - 0.5;
	static double yCenter = dist_to_table - table_length/2.0;

	if (!status) {
		valid = false;
	}
	else if (blob(Z) > zMax) {
		valid = false;
	}
	else if (blob(Y) > yMax || blob(Y) < yMin) {
		valid = false;
	}
	// between the table if ball appears outside the x limits
	else if (fabs(blob(Y) - yCenter) < table_length/2.0 &&
			fabs(blob(X)) > xMax) {
		valid = false;
	}
	// on the table blob should not appear under the table
	else if (fabs(blob(X)) < xMax && fabs(blob(Y) - yCenter) < table_length/2.0
			&& blob(Z) < zMin) {
		valid = false;
	}
	else {
		valid = true;
	}
	return valid;
}
