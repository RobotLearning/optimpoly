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

RBF::RBF(const std::string &param_file) {

  std::ifstream stream(param_file);
  json jobs;
  stream >> jobs;
  centers = json2vec(jobs.at("centers"));
  widths = json2vec(jobs.at("widths"));
  intercepts = json2vec(jobs.at("intercepts"));
  int idx = 0;
  params = zeros<mat>(centers.n_elem, 7);
  for (auto elem : jobs.at("joints")) {
    idx = elem.at("ID");
    params.col(idx) = json2vec(elem.at("params"));
  }
}

void RBF::basis_fnc(const double T, joint &Q) const {
  vec xt = exp(-pow(T - centers, 2) / widths);
  Q.q = params.t() * xt;
  Q.qd = params.t() * ((-2 * (T - centers) / widths) % xt);
  Q.qdd =
      params.t() * ((-2 / widths + 4 * pow((T - centers) / widths, 2)) % xt);
}

void RBF::step(const double &dt, joint &Q) {
  t += dt;
  basis_fnc(t, Q);
  Q.q += intercepts;
}

void RBF::get_init_pos(vec7 &pos) const {
  joint Q;
  basis_fnc(0.0, Q);
  pos = Q.q + intercepts;
}

void RBF::set_init_pos(const vec7 &pos) {
  vec7 pos_pre;
  get_init_pos(pos_pre);
  intercepts += pos - pos_pre;
}

void RBF::get_goal_pos(vec7 &pos) const {
  joint Q;
  basis_fnc(Tmax, Q);
  pos = Q.q + intercepts;
}

CRBF::CRBF(const std::string &param_file) {

  std::ifstream stream(param_file);
  json jobs;
  stream >> jobs;
  params = json2vec(jobs.at("params"));
  centers = zeros<mat>(params.n_elem, 7);
  widths = zeros<mat>(params.n_elem, 7);

  int idx = 0;
  for (auto elem : jobs.at("joints")) {
    idx = elem.at("ID");
    intercepts(idx) = elem.at("intercept");
    centers.col(idx) = json2vec(elem.at("centers"));
    widths.col(idx) = json2vec(elem.at("widths"));
  }
}

void CRBF::basis_fnc(const double T, joint &Q) const {

  for (int i = 0; i < NDOF; i++) {
    vec x = exp(-pow(T - centers.col(i), 2) / widths.col(i));
    Q.q(i) = dot(params, x);
    Q.qd(i) = dot(params, (-2 * (T - centers.col(i)) / widths.col(i)) % x);
    Q.qdd(i) = dot(params, ((-2 / widths.col(i) +
                             4 * pow((T - centers.col(i)) / widths.col(i), 2)) %
                            x));
  }
}

void CRBF::step(const double &dt, joint &Q) {
  t += dt;
  basis_fnc(t, Q);
  Q.q += intercepts;
}

void CRBF::get_init_pos(vec7 &pos) const {
  joint Q;
  basis_fnc(0.0, Q);
  pos = Q.q + intercepts;
}

void CRBF::set_init_pos(const vec7 &pos) {
  vec7 pos_pre;
  get_init_pos(pos_pre);
  intercepts += pos - pos_pre;
}

void CRBF::get_goal_pos(vec7 &pos) const {
  joint Q;
  basis_fnc(Tmax, Q);
  pos = Q.q + intercepts;
}

} // namespace serve
