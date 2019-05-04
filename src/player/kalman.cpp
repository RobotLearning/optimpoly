/**
 * @file kalman.cpp
 *
 * @brief (Discrete) Kalman filtering class KF in C_++
 *
 * Kalman Filter class for basic filtering in C_++ using ARMADILLO linear
 * algebra library.
 * Compatible with continuous models by using discretize() method.
 * Can be extended easily (e.g. see EKF).
 *
 * TODO: add square-root form for stability!
 *
 *  Created on: Jan 25, 2017
 *      Author: okoc
 */

#include "kalman.h"
#include <armadillo>
#include <iostream>
#include <string>

using namespace arma;

namespace player {

KF::KF(mat &Cin, mat &Qin, mat &Rin)
    : A_(datum::inf * ones<mat>(Cin.n_cols, Cin.n_cols)),
      B_(datum::inf * ones<mat>(Cin.n_cols, Cin.n_cols)), C_(Cin), Q_(Qin),
      R_(Rin),                                             // init model
      x_(datum::inf * ones<vec>(Cin.n_cols)),              // init state
      P_(datum::inf * ones<mat>(Cin.n_cols, Cin.n_cols)) { // init covar

  check_models(A_, B_, C_);
  // checking for correct noise covariances
  check_spd(Qin);
  check_spd(Rin);
}

KF::KF(mat &Ain, mat &Bin, mat &Cin, mat &Qin, mat &Rin)
    : A_(Ain), B_(Bin), C_(Cin), Q_(Qin), R_(Rin),         // init model
      x_(datum::inf * ones<vec>(Ain.n_rows)),              // init state
      P_(datum::inf * ones<mat>(Ain.n_rows, Ain.n_rows)) { // init covar
  check_models(Ain, Bin, Cin);
  // checking for correct noise covariances
  check_spd(Qin);
  check_spd(Rin);
}

KF::KF(vec &x0, mat &P0, mat &Ain, mat &Bin, mat &Cin, mat &Qin, mat &Rin)
    : A_(Ain), B_(Bin), C_(Cin), Q_(Qin), R_(Rin), // init model
      x_(x0), P_(P0) {                             // init state

  check_models(Ain, Bin, Cin);
  // checking for correct noise covariances
  check_spd(Qin);
  check_spd(Rin);
}

void KF::set_prior(const vec &x0, const mat &P0) {

  x_ = x0;
  P_ = P0;
}

void KF::check_models(const mat &Ain, const mat &Bin, const mat &Cin) const {

  // checking for correct model matrix sizes
  if (Ain.n_rows != Ain.n_cols)
    cerr << "A_ must be square!" << endl;
  if (Bin.n_rows != Ain.n_rows)
    cerr << "A_ and B_ must have same row size!" << endl;
  if (Cin.n_cols != Ain.n_cols)
    cerr << "C_ and A_ must have same column size!" << endl;
  // if (Cin.n_rows != Din.n_rows)
  //	cerr << "C_ and D must have same row size!" << endl;
}

void KF::check_spd(const mat &M) const {

  if (M.n_cols != M.n_rows)
    cerr << "Covariance matrix must be square!" << endl;
  if (any(vectorise(M - M.t())))
    cerr << "Covariance matrix must be symmetric!" << endl;
  const std::string form = "sa";
  vec eigvals = eig_sym(M);
  if (eigvals(0) < 0.0) {
    throw "Covariance matrix must be positive semidefinite!";
  }
}

mat KF::chol_semi(const mat &M) const {

  int size;
  mat Y;
  try {
    size = M.n_cols;
    check_spd(M);
    Y = chol(M);
  } catch (char *err1) {
    std::cout << "Covariance is reset to zero!" << std::endl;
    Y = zeros<mat>(size, size);
  } catch (std::runtime_error &err2) {
    vec eigval;
    mat eigvec;
    eig_sym(eigval, eigvec, M);
    Y = eigvec * diagmat(sqrt(eigval));
  }
  return Y;
}

void KF::discretize(const mat &Ac, const mat &Bc, double dt) {

  /*! MATLAB code
      trick to get discrete time versions
      Mat = [obj.A_, obj.B_; zeros(dimu, dimx + dimu)];
      MD = expm(h * Mat);
      Ad = MD(1:dimx,1:dimx);
      Bd = MD(1:dimx,dimx+1:end); */

  int dimu = Bc.n_cols;
  int dimx = x_.n_elem;
  mat M = join_vert(join_horiz(Ac, Bc), zeros<mat>(dimu, dimx + dimu));
  mat MD = expmat(dt * M);

  A_ = MD(span(0, dimx - 1), span(0, dimx - 1));
  B_ = MD(span(0, dimx - 1), span(dimx, dimx + dimu - 1));
}

vec KF::get_mean() const {

  // quick and dirty check for initialization
  if (!x_.is_finite()) {
    throw std::runtime_error("KF not initialized! Please set prior!");
  }

  return x_;
}

mat KF::get_covar() const {

  // quick and dirty check for initialization
  if (P_(0, 0) == datum::inf) {
    throw std::runtime_error("KF not initialized! Please set prior!");
  }

  return P_;
}

mat KF::get_model(int idx) const {

  if (idx > 4 || idx < 1) {
    throw std::runtime_error("Index must be 1 to 4 only!");
  }
  mat out;
  switch (idx) {
  case 1:
    out = A_;
    break;
  case 2:
    out = B_;
    break;
  case 3:
    out = C_;
    break;
  case 4:
    throw std::runtime_error("D matrix unimplemented!");
    break;
  default:
    throw std::runtime_error("This shouldn't happen!");
  }
  if (out(0, 0) == datum::inf) {
    throw std::runtime_error("Matrix is not initialized!");
  }
  return out;
}

void KF::predict() {

  x_ = A_ * x_;
  P_ = A_ * P_ * A_ + Q_;
}

void KF::predict(const vec &u) {

  if (u.n_elem != B_.n_cols)
    cerr << "Control input u should have the right size!" << endl;

  x_ = A_ * x_ + B_ * u;
  P_ = A_ * P_ * A_ + Q_;
}

void KF::update(const vec &y) {

  if (y.n_elem != C_.n_rows)
    cerr << "Observation y should have the right size!" << endl;

  // innovation sequence
  vec z = y - C_ * x_;
  // innovation covariance
  mat Theta = C_ * P_ * C_.t() + R_;
  // optimal Kalman gain
  mat K = P_ * C_.t() * inv(Theta);

  // update state mean and covariance
  P_ = P_ - K * C_ * P_;
  // cout << "obs:" << "\t" << y.t();
  // cout << "x_pre:" << "\t" << x_.t();
  x_ = x_ + K * z;
  // cout << "x_post:" << "\t" << x_.t();
}

mat KF::sample_observations(int N) const {

  int dimy = C_.n_rows;
  int dimx = x_.n_elem;
  vec x_next = x_;
  mat Y(dimy, N, fill::zeros);
  vec w(dimy, fill::randn);

  Y.col(0) = C_ * x_ + chol(R_) * w;
  for (int i = 1; i < N; i++) {
    w = randn<vec>(dimy);
    Y.col(i) = C_ * x_next + chol_semi(R_) * w;
    w = randn<vec>(dimx);
    x_next = A_ * x_next + chol_semi(Q_) * w;
  }

  return Y;
}

} // namespace player
