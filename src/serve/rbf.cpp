/*
 * rbf.cpp
 *
 *  Created on: Sep 1, 2018
 *      Author: robolab
 */

#include "rbf.h"
#include "json.hpp"
#include "optim.h"
#include "utils.h"
#include <armadillo>

using json = nlohmann::json;
using namespace optim;
using namespace arma;

namespace serve {

RBF::RBF()
    : t_(0.0), Tmax_(1.0), intercepts_(zeros<vec>(7)), widths_(zeros<vec>(10)),
      centers_(zeros<vec>(10)), params_(zeros<mat>(10, 7)) {}

RBF::RBF(const std::string &param_file)
    : t_(0.0), Tmax_(1.0), intercepts_(zeros<vec>(7)), widths_(zeros<vec>(10)),
      centers_(zeros<vec>(10)), params_(zeros<mat>(10, 7)) {

  std::ifstream stream(param_file);
  json jobs;
  stream >> jobs;
  centers_ = json2vec(jobs.at("centers"));
  widths_ = json2vec(jobs.at("widths"));
  intercepts_ = json2vec(jobs.at("intercepts"));
  int idx = 0;
  params_ = zeros<mat>(centers_.n_elem, 7);
  for (auto elem : jobs.at("joints")) {
    idx = elem.at("ID");
    params_.col(idx) = json2vec(elem.at("params"));
  }
}

void RBF::basis_fnc(const double T, joint &Q) const {
  vec xt = exp(-pow(T - centers_, 2) / widths_);
  Q.q = params_.t() * xt;
  Q.qd = params_.t() * ((-2 * (T - centers_) / widths_) % xt);
  Q.qdd =
      params_.t() * ((-2 / widths_ + 4 * pow((T - centers_) / widths_, 2)) % xt);
}

void RBF::step(const double &dt, joint &Q) {
  t_ += dt;
  basis_fnc(t_, Q);
  Q.q += intercepts_;
}

void RBF::get_init_pos(vec7 &pos) const {
  joint Q;
  basis_fnc(0.0, Q);
  pos = Q.q + intercepts_;
}

void RBF::set_init_pos(const vec7 &pos) {
  vec7 pos_pre;
  get_init_pos(pos_pre);
  intercepts_ += pos - pos_pre;
}

void RBF::get_goal_pos(vec7 &pos) const {
  joint Q;
  basis_fnc(Tmax_, Q);
  pos = Q.q + intercepts_;
}

CRBF::CRBF()
    : t_(0.0), Tmax_(1.0), intercepts_(zeros<vec>(7)), widths_(zeros<mat>(1, 7)),
      centers_(zeros<mat>(1, 7)), params_(zeros<vec>(1)) {}

CRBF::CRBF(const std::string &param_file)
    : t_(0.0), Tmax_(1.0), intercepts_(zeros<vec>(7)), widths_(zeros<mat>(1, 7)),
      centers_(zeros<mat>(1, 7)), params_(zeros<vec>(1)) {

  std::ifstream stream(param_file);
  json jobs;
  stream >> jobs;
  params_ = json2vec(jobs.at("params"));
  centers_ = zeros<mat>(params_.n_elem, 7);
  widths_ = zeros<mat>(params_.n_elem, 7);

  int idx = 0;
  for (auto elem : jobs.at("joints")) {
    idx = elem.at("ID");
    intercepts_(idx) = elem.at("intercept");
    centers_.col(idx) = json2vec(elem.at("centers"));
    widths_.col(idx) = json2vec(elem.at("widths"));
  }
}

void CRBF::basis_fnc(const double T, joint &Q) const {

  for (int i = 0; i < NDOF; i++) {
    vec x = exp(-pow(T - centers_.col(i), 2) / widths_.col(i));
    Q.q(i) = dot(params_, x);
    Q.qd(i) = dot(params_, (-2 * (T - centers_.col(i)) / widths_.col(i)) % x);
    Q.qdd(i) = dot(params_, ((-2 / widths_.col(i) +
                             4 * pow((T - centers_.col(i)) / widths_.col(i), 2)) %
                            x));
  }
}

void CRBF::step(const double &dt, joint &Q) {
  t_ += dt;
  basis_fnc(t_, Q);
  Q.q += intercepts_;
}

void CRBF::get_init_pos(vec7 &pos) const {
  joint Q;
  basis_fnc(0.0, Q);
  pos = Q.q + intercepts_;
}

void CRBF::set_init_pos(const vec7 &pos) {
  vec7 pos_pre;
  get_init_pos(pos_pre);
  intercepts_ += pos - pos_pre;
}

void CRBF::get_goal_pos(vec7 &pos) const {
  joint Q;
  basis_fnc(Tmax_, Q);
  pos = Q.q + intercepts_;
}

} // namespace serve
