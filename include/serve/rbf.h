#pragma once
#include "optim.h"
#include <armadillo>

using arma::mat;
using arma::vec;
using arma::vec7;
using arma::zeros;

namespace serve {

class RBF { // rbf for each dof independent

private:
  double t_;
  double Tmax_;
  vec intercepts_;
  vec widths_;
  vec centers_;
  mat params_;

  void basis_fnc(const double T, optim::joint &) const;

public:
  RBF();
  RBF(const std::string &param_file);

  void reset() { t_ = 0.0; };
  void step(const double &dt, optim::joint &);
  void get_init_pos(vec7 &pos) const;
  void set_init_pos(const vec7 &pos);
  void get_goal_pos(vec7 &pos) const;
};

class CRBF { // coupled RBF where the dofs are dependent

private:
  double t_;
  double Tmax_;
  vec intercepts_;
  mat widths_;
  mat centers_;
  vec params_;

  // stacked basis function
  void basis_fnc(const double T, optim::joint &) const;

public:
  CRBF();
  CRBF(const std::string &param_file);

  void reset() { t_ = 0.0; };
  void step(const double &dt, optim::joint &);
  void get_init_pos(vec7 &pos) const;
  void set_init_pos(const vec7 &pos);
  void get_goal_pos(vec7 &pos) const;
};

} // namespace serve
