/*
 * optim.cpp
 *
 *  Created on: Apr 9, 2018
 *      Author: okoc
 */

#include "optim.h"
#include "constants.h"
#include "kinematics.h"
#include "lookup.h"
#include "math.h"
#include "stdlib.h"
#include "tabletennis.h"
#include "utils.h"
#include <armadillo>

namespace optim {

/*
 * Give info about the optimization after termination
 */
static bool check_optim_result(const int res);

Optim::~Optim() { nlopt_destroy(opt_); }

Optim::Optim()
    : lookup_(false), verbose_(true), moving_(false), update_(false),
      running_(false), detach_(false),
      lookup_table_(zeros<mat>(1, OPTIM_DIM_ + 2 * NCART)), opt_(nullptr),
      qf_{0.0}, qfdot_{0.0}, T_(1.0), param_des_(nullptr), lb_{0.0}, ub_{0.0},
      qrest_{0.0}, q0_{0.0}, q0dot_{0.0}, time2return_(1.0) {}

void Optim::update_init_state(const joint &qact) {
  for (int i = 0; i < NDOF; i++) {
    q0_[i] = qact.q(i);
    q0dot_[i] = qact.qd(i);
  }
}

bool Optim::check_running() { return running_; }

bool Optim::check_update() { return update_; }

void Optim::set_moving(bool flag_move) { moving_ = flag_move; }

void Optim::set_detach(bool flag_detach) { detach_ = flag_detach; }

void Optim::set_return_time(const double &ret_time) { time2return_ = ret_time; }

void Optim::set_verbose(bool flag_verbose) { verbose_ = flag_verbose; }

bool Optim::get_params(const joint &qact, spline_params &p) {

  bool flag = false;
  if (update_ && !running_) {
    vec7 qf, qfdot, qrest;
    for (int i = 0; i < NDOF; i++) {
      qf(i) = qf_[i];
      qfdot(i) = qfdot_[i];
      qrest(i) = qrest_[i];
    }
    vec7 qnow = qact.q;
    vec7 qdnow = qact.qd;
    p.a.col(0) = 2.0 * (qnow - qf) / pow(T_, 3) + (qfdot + qdnow) / pow(T_, 2);
    p.a.col(1) = 3.0 * (qf - qnow) / pow(T_, 2) - (qfdot + 2.0 * qdnow) / T_;
    p.a.col(2) = qdnow;
    p.a.col(3) = qnow;
    // cout << "A = \n" << p.a << endl;
    p.b.col(0) = 2.0 * (qf - qrest) / pow(time2return_, 3) +
                 (qfdot) / pow(time2return_, 2);
    p.b.col(1) = 3.0 * (qrest - qf) / pow(time2return_, 2) -
                 (2.0 * qfdot) / time2return_;
    p.b.col(2) = qfdot;
    p.b.col(3) = qf;
    p.time2hit = T_;
    // cout << "B = \n" << p.b << endl;
    flag = true;
    update_ = false;
  }
  return flag;
}

void Optim::update_rest_state(const vec7 &q_rest_new) {

  for (int i = 0; i < NDOF; i++)
    qrest_[i] = q_rest_new(i);
}

void Optim::set_des_params(optim_des *params_) { param_des_ = params_; }

void Optim::init_lookup_soln(double *x) {

  vec::fixed<OPTIM_DIM_> robot_params;
  vec6 ball_params;
  for (int i = 0; i < NCART; i++) {
    ball_params(i) = param_des_->ball_pos(i, 0);
    ball_params(i + NCART) = param_des_->ball_vel(i, 0);
  }
  // cout << "Init ball est:" << ball_params << endl;
  player::predict_till_net(ball_params);
  // cout << "Net ball est:" << ball_params << endl;
  // k = 5 nearest neighbour regression
  player::knn(lookup_table_, ball_params, 5, robot_params);
  for (int i = 0; i < OPTIM_DIM_; i++) {
    x[i] = robot_params(i);
    //  printf("x[%d] = %f\n", i, x[i]);
  }
}

void Optim::run() {

  std::thread t = std::thread(&Optim::optim, this);
  if (detach_) {
    t.detach();
  } else {
    t.join();
  }
}

void Optim::optim() {

  update_ = false;
  running_ = true;
  double x[OPTIM_DIM_];

  if (moving_) {
    init_last_soln(x);
  } else {
    if (lookup_) {
      if (verbose_) {
        std::cout
            << "Looking up good initial parameters with k = 5\n"; // kNN
                                                                  // parameter
                                                                  // k = 5
      }
      init_lookup_soln(x);
    } else {
      init_rest_soln(x);
    }
  }

  double init_time = get_time();
  double past_time = 0.0;
  double minf; // the minimum objective value, upon return //
  int res;     // error code

  if ((res = nlopt_optimize(opt_, x, &minf)) < 0) {
    past_time = (get_time() - init_time) / 1e3;
    if (verbose_) {
      printf("NLOPT failed with exit code %d!\n", res);
      printf("NLOPT took %f ms\n", past_time);
    }
  } else {
    past_time = (get_time() - init_time) / 1e3;
    if (verbose_) {
      printf("NLOPT success with exit code %d!\n", res);
      printf("NLOPT took %f ms\n", past_time);
      printf("Found minimum at f = %0.10g\n", minf);
    }
    if (test_soln(x) < 1e-2)
      finalize_soln(x, past_time);
  }
  if (verbose_)
    check_optim_result(res);
  running_ = false;
}

static bool check_optim_result(const int res) {

  bool flag = false;
  switch (res) {
  case NLOPT_SUCCESS:
    printf("Success!\n");
    flag = true;
    break;
  case NLOPT_STOPVAL_REACHED:
    printf("Optimization stopped because stopval (above) was reached.\n");
    flag = true;
    break;
  case NLOPT_FTOL_REACHED:
    printf("Optimization stopped because ftol_rel "
           "or ftol_abs (above) was reached.\n");
    flag = true;
    break;
  case NLOPT_XTOL_REACHED:
    flag = true;
    printf("Optimization stopped because xtol_rel or xtol_abs (above) was "
           "reached.\n");
    break;
  case NLOPT_MAXEVAL_REACHED:
    flag = true;
    printf("Optimization stopped because maxeval (above) was reached.\n");
    break;
  case NLOPT_MAXTIME_REACHED:
    flag = true;
    printf("Optimization stopped because maxtime (above) was reached.\n");
    break;
  case NLOPT_FAILURE:
    printf("Epic fail!\n");
    break;
  case NLOPT_INVALID_ARGS:
    printf("Invalid arguments (e.g. lower bounds are bigger than "
           "upper bounds, an unknown algorithm was specified, etcetera).\n");
    break;
  case NLOPT_OUT_OF_MEMORY:
    printf("Ran out of memory!\n");
    break;
  case NLOPT_ROUNDOFF_LIMITED:
    printf("Halted because roundoff errors limited progress."
           "(In this case, the optimization still typically returns a useful "
           "result.\n");
    break;
  case NLOPT_FORCED_STOP:
    printf("Halted because of a forced termination: "
           "the user called nlopt_force_stop(opt_)"
           "on the optimization’s nlopt_opt object "
           "opt_ from the user’s objective function or constraints.\n");
    break;
  }
  return flag;
}

} // namespace optim
