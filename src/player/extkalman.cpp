/**
 * @file extkalman.cpp
 *
 * @brief (Discrete) Extended Kalman filtering class EKF in C++
 *
 * Kalman Filter class for basic filtering in C++ using ARMADILLO linear
 * algebra library.
 * Compatible with continuous models by using discretize() method.
 *
 *  Created on: Feb 2, 2017
 *      Author: okoc
 */
#include "kalman.h"
#include "tabletennis.h"
#include <armadillo>

using namespace arma;

namespace player {

EKF::EKF(model fp, mat &Cin, mat &Qin, mat &Rin, double rej_mult)
    : KF(Cin, Qin, Rin), fparams_(nullptr), outlier_reject_mult_(rej_mult),
      f_(fp) {}

EKF::EKF(const EKF &ekf_)
    : KF(ekf_), fparams_(ekf_.fparams_),
      outlier_reject_mult_(ekf_.outlier_reject_mult_), f_(ekf_.f_) {}

EKF &EKF::operator=(const EKF &ekf_) {
  if (this != &ekf_) {
    fparams_ = ekf_.fparams_;
    f_ = ekf_.f_;
    outlier_reject_mult_ = ekf_.outlier_reject_mult_;
  }
  return *this;
}

mat EKF::linearize(const double dt, const double h) const {

  int dimx = x_.n_elem;
  static mat delta = h * eye<mat>(dimx, dimx);
  static mat dfdx = zeros<mat>(dimx, dimx);
  static mat fPlus = zeros<mat>(dimx, dimx);
  static mat fMinus = zeros<mat>(dimx, dimx);
  for (int i = 0; i < dimx; i++) {
    fPlus.col(i) = this->f_(x_ + delta.col(i), dt, fparams_);
    fMinus.col(i) = this->f_(x_ - delta.col(i), dt, fparams_);
  }
  return (fPlus - fMinus) / (2 * h);
}

void EKF::predict(const double dt, const bool lin_flag) {

  x_ = this->f_(x_, dt, fparams_);
  if (lin_flag) {
    // cout << "A = \n" << linearize(dt,0.01);
    A_ = linearize(dt, 0.0001);
    P_ = A_ * P_ * A_.t() + Q_;
    // cout << "P = \n" << P << "A = \n" << A;
  }
}

mat EKF::predict_path(const double dt, const int N) const {

  int dimx = x_.n_elem;
  mat XX(dimx, N);
  EKF filter_pred = *this;

  for (int i = 0; i < N; i++) {
    filter_pred.predict(dt, false);
    XX.col(i) = filter_pred.get_mean();
  }

  return XX;
}

bool EKF::check_outlier(const vec &y, const bool verbose) const {

  static int dim_y = y.n_elem;
  bool outlier = true;
  vec threshold = outlier_reject_mult_ * arma::sqrt(P_.diag());
  vec inno = clamp(y - C_ * x_, -10, 10);
  // std::cout << "INNO:" << inno.t() << "THRESH:" << threshold.t();

  if (all(abs(inno) < threshold.head(dim_y))) {
    outlier = false;
  } else {
    if (verbose) {
      std::cout << "Outlier:"
                << y.t()
                //<< "State: \t" << x_.t()
                << "Inno:  " << inno.t()
                << "Thresh: " << threshold.head(dim_y).t();
    }
  }
  return outlier;
}

bool check_new_obs(const vec3 &obs, double tol) {

  static vec3 last_obs = zeros<vec>(3);

  if (norm(obs - last_obs) > tol) {
    last_obs = obs;
    return true;
  }
  return false;
}

bool check_reset_filter(const bool newball, const int verbose,
                        const double threshold) {

  bool reset = false;
  static int reset_cnt = 0;
  static bool firsttime = true;
  static wall_clock timer;

  if (firsttime) {
    firsttime = false;
    timer.tic();
  }

  if (newball) {
    if (timer.toc() > threshold) {
      reset = true;
      if (verbose > 0) {
        std::cout << "Resetting filter! Count: " << ++reset_cnt << std::endl;
      }
    }
    timer.tic();
  }
  return reset;
}

EKF init_ball_filter(const double var_model, const double var_noise,
                     const bool spin, const double out_reject_mult,
                     const double *topspin) {

  mat C = eye<mat>(3, 6);
  mat66 Q = var_model * eye<mat>(6, 6);
  mat33 R = var_noise * eye<mat>(3, 3);
  if (spin) {
    EKF filter = EKF(calc_spin_ball, C, Q, R, out_reject_mult);
    filter.set_fun_params((void *)topspin);
    return filter;
  } else {
    // cout << "Init spin-free model\n";
    EKF filter = EKF(calc_next_ball, C, Q, R, out_reject_mult);
    return filter;
  }
}

} // namespace player
