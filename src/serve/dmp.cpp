/*
 * Dynamic Movement Primitives
 *
 * @author Okan Koc
 * @update Fri, Jun 8, 2018
 */

#include "dmp.h"
#include "constants.h"
#include "json.hpp"
#include <armadillo>

using namespace arma;
using json = nlohmann::json;

namespace serve {

void Canonical::reset() { phase = 1.0; }

void Canonical::step(const double &dt) { phase -= dt * a * tau * phase; }

DMP::DMP(std::vector<double> weights, std::vector<double> centers,
         std::vector<double> heights, double alpha, double beta, double goal,
         double init_pos)
    : x0_(zeros<vec>(3)), x_(zeros<vec>(3)), w_(zeros<vec>(10)),
      h_(zeros<vec>(10)), c_(zeros<vec>(10)), alpha_(alpha), beta_(beta),
      goal_(1.0), SAFE_ACC_(false) {

  set_weights(weights, centers, heights);
  set_goal_state(goal);
  set_init_state(init_pos);
  reset();
}

void DMP::reset() { x_ = x0_; }

void DMP::set_weights(const std::vector<double> &w_stdvec,
                      const std::vector<double> &c_stdvec,
                      const std::vector<double> &h_stdvec) {
  w_ = conv_to<vec>::from(w_stdvec);
  c_ = conv_to<vec>::from(c_stdvec);
  h_ = conv_to<vec>::from(h_stdvec);
}

void DMP::set_goal_state(const double &goal) { goal_ = goal; }

void DMP::set_init_state(const double &init_pos) {
  x0_(0) = init_pos;
  x0_(1) = 0.0;
  x0_(2) = 0.0;
}

void DMP::get_goal_state(double &goal) const { goal = goal_; }

void DMP::get_init_state(double &init_pos) const { init_pos = x0_(0); }

vec3 DMP::step(const Canonical &can_, const double &dt) {

  const double amp = 1.0;
  double f = forcing(can_.phase);

  // OLD VERSION written as a multivar. system
  /*mat22 A = {{0,1}, {-alpha_*beta_*can_.tau*can_.tau, -alpha_*can_.tau}};
  vec2 b = {0,can_.tau*can_.tau*(alpha_*beta_*goal_ + amp*f)};
  vec2 xdot = A * x_.head(2) + b;
  x_(0) += dt * xdot(0);
  x_(1) += dt * xdot(1);
  x_(2) = xdot(1);*/

  // compute accelerations
  x_(2) = can_.tau * can_.tau * alpha_ * beta_ * (goal_ - x_(0)) -
          can_.tau * alpha_ * x_(1);
  if (SAFE_ACC_)
    x_(2) *= (1 - can_.phase);
  x_(2) += can_.tau * can_.tau * amp * f;
  x_(0) += dt * x_(1);
  x_(1) += dt * x_(2);

  return x_;
}

vec DMP::basis(const double &x_) const {
  return exp(-h_ % ((x_ - c_) % (x_ - c_)));
}

double DMP::forcing(const double &phase) const {

  // double f = 0.0;
  vec psi = basis(phase);
  /*for (unsigned int i = 0; i < w_.n_elem; i++) {
      f += psi(i) * w_(i) * phase;
  }*/
  double f = phase * sum(psi % w_);
  double scale = sum(psi);
  f /= scale + w_.n_elem * 1e-10;
  return f;
}

void Joint_DMPs::reset() {
  can_.reset();
  for (unsigned int i = 0; i < dmps_.size(); i++) {
    dmps_[i].reset();
  }
}

Joint_DMPs::Joint_DMPs(const std::string &param_file) {

  // load from json file
  // cout << param_file << endl;
  std::ifstream stream(param_file);
  json jobs;
  stream >> jobs;
  can_.tau = jobs.at("tau");
  for (auto elem : jobs.at("joints")) {
    DMP dmp = DMP(elem.at("weights"), elem.at("centers"), elem.at("heights"),
                  jobs.at("alpha"), jobs.at("beta"), elem.at("goal"),
                  elem.at("init_pos"));
    dmps_.push_back(dmp);
  }
}

void Joint_DMPs::step(const double &dt, joint &Q) {

  for (unsigned int i = 0; i < dmps_.size(); i++) {
    vec x_ = dmps_[i].step(can_, dt);
    Q.q(i) = x_(0);
    Q.qd(i) = x_(1);
    Q.qdd(i) = x_(2);
  }
  can_.step(dt);
}

mat Joint_DMPs::evolve(const double &T) {

  using namespace const_tt;
  unsigned int N = T / DT;
  reset();
  mat q_evolve = zeros<mat>(NDOF, N);
  joint Q;
  for (unsigned int i = 0; i < N; i++) {
    step(DT, Q);
    q_evolve.col(i) = Q.q;
  }
  reset();
  return q_evolve;
}

void Joint_DMPs::evolve(const double &T, mat &Q, mat &Qd, mat &Qdd) {
  using namespace const_tt;
  unsigned N = T / DT;
  reset();
  Q = zeros<mat>(NDOF, N);
  Qd = zeros<mat>(NDOF, N);
  Qdd = zeros<mat>(NDOF, N);
  joint J;
  for (unsigned i = 0; i < N; i++) {
    step(DT, J);
    Q.col(i) = J.q;
    Qd.col(i) = J.qd;
    Qdd.col(i) = J.qdd;
  }
  reset();
}

void Joint_DMPs::get_init_pos(vec &pos) const {

  using namespace std;
  try {
    for (unsigned int i = 0; i < dmps_.size(); i++) {
      dmps_[i].get_init_state(pos(i));
    }
  } catch (std::exception &ex) {
    cerr << "Array length incorrect: " << ex.what() << endl;
  }
}

void Joint_DMPs::set_init_pos(const vec &pos) {

  using namespace std;
  try {
    for (unsigned int i = 0; i < dmps_.size(); i++) {
      dmps_[i].set_init_state(pos(i));
    }
  } catch (std::exception &ex) {
    cerr << "Array length incorrect: " << ex.what() << endl;
  }
}

void Joint_DMPs::set_goal_pos(const vec &pos) {

  using namespace std;
  try {
    for (unsigned int i = 0; i < dmps_.size(); i++) {
      dmps_[i].set_goal_state(pos(i));
    }
  } catch (std::exception &ex) {
    cerr << "Array length incorrect: " << ex.what() << endl;
  }
}

void Joint_DMPs::get_goal_pos(vec &pos) const {

  using namespace std;
  try {
    for (unsigned int i = 0; i < dmps_.size(); i++) {
      dmps_[i].get_goal_state(pos(i));
    }
  } catch (std::exception &ex) {
    cerr << "Array length incorrect: " << ex.what() << endl;
  }
}

double Joint_DMPs::get_time_constant() const { return can_.tau; }

void Joint_DMPs::set_time_constant(const double &tau) { can_.tau = tau; }

void Joint_DMPs::turn_on_safe_acc() {
  for (unsigned int i = 0; i < dmps_.size(); i++) {
    dmps_[i].SAFE_ACC_ = true;
  }
}

} // namespace serve
