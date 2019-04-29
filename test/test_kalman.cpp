/*
 * Test cases for the Kalman Filter C++ implementation.

 * Kalman Filtering unit test suite
 * for testing various use cases.
 *
 *  Created on: Jan 25, 2017
 *      Author: okoc
 */

#include "kalman.h"
#include "player.hpp"
#include "tabletennis.h"
#include "gtest/gtest.h"
#include <armadillo>

using namespace std;
using namespace arma;
using namespace const_tt;
using namespace player;

static bool test_uninit_exception(const KF &filter);

TEST(KalmanFilterTests, TestKalmanFilterInitialization) {

  // init Kalman Filter
  int dimx = 2;
  int dimu = 1;
  int dimy = 1;
  // model matrices
  mat A(dimx, dimx, fill::randu);
  mat B(dimx, dimu, fill::randu);
  mat C(dimy, dimx, fill::zeros);
  C(0, 0) = 1.0;
  // diagonal covariances
  mat Q(dimx, dimx, fill::eye);
  mat R(dimy, dimy, fill::eye);
  double s2x = 0.01;
  double s2y = 0.10;
  Q *= sqrt(s2x);
  R *= sqrt(s2y);
  KF filter = KF(A, B, C, Q, R);
  EXPECT_TRUE(test_uninit_exception(filter));
}

TEST(KalmanFilterTests, RunExtendedKalmanFilterInitializationAndPrediction) {

  int dimx = 6;
  int dimy = 3;
  double eps = 0.01;

  mat C(dimy, dimx, fill::zeros);
  C.cols(1, dimy) = eye<mat>(dimy, dimy);
  // diagonal covariances
  mat Q(dimx, dimx, fill::zeros);
  mat R(dimy, dimy, fill::eye);
  R *= eps;
  EKF filter = EKF(calc_next_ball, C, Q, R);

  vec x0 = zeros<vec>(dimx);
  mat P0(dimx, dimx, fill::eye);
  P0 *= eps;
  filter.set_prior(x0, P0);
  filter.predict(0.01);
}

/*
 * Test predict path function of EKF with table tennis
 *
 */
TEST(KalmanFilterTests, ComparePathPredictionToActualValuesForTableTennisBall) {

  TableTennis tt = TableTennis();

  // initialize a filter to predict
  EKF filter = init_ball_filter();

  int N = 1000;
  double dt = 0.002;
  tt.set_ball_gun(0.2);

  vec6 ball_state = tt.get_ball_state();
  filter.set_prior(ball_state, 0.01 * eye<mat>(6, 6));

  mat M = zeros<mat>(6, N);
  for (int i = 0; i < N; i++) {
    tt.integrate_ball_state(dt);
    M.col(i) = join_vert(tt.get_ball_position(), tt.get_ball_velocity());
  }
  mat Mfilt = filter.predict_path(dt, N);

  EXPECT_TRUE(approx_equal(M, Mfilt, "absdiff", 0.002));
}

static bool test_uninit_exception(const KF &filter) {

  bool flag = false;
  try {
    vec mean = filter.get_mean();
  } catch (const runtime_error &exception) {
    cout << "Caught the exception!" << endl;
    flag = true;
  }

  return flag;
}
