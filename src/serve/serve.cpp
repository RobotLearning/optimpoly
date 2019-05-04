/*
 * serve.cpp
 *
 *  Created on: Jul 22, 2018
 *      Author: okoc
 */
#include "serve.h"

#include <algorithm>
#include <armadillo>
#include <string>

#include "ball_interface.h"
#include "dmp.h"
#include "kalman.h"
#include "kinematics.hpp"
#include "player.hpp"
#include "rbf.h"
#include "tabletennis.h"

using namespace arma;
using namespace const_tt;
using namespace player;
using namespace optim;

namespace serve {

using dmps = Joint_DMPs;
using vec_str = std::vector<std::string>;

dmps init_dmps(std::string &file) {
  const std::string home = std::getenv("HOME");
  vec_str files = get_files(home + "/projects/table-tennis/json/", "dmp");
  arma_rng::set_seed_random();
  int val = (randi(1, distr_param(0, files.size() - 1)).at(0));
  file = files[val];
  std::cout << "Loading DMP " << file << std::endl;
  std::string full_file = home + "/projects/table-tennis/json/" + file;
  return dmps(full_file);
}

template <typename T>
ServeBall<T>::ServeBall(const serve_flags &sflags, const vec7 &qinit)
    : ran_optim_(false), init_filter_(false), idx_balls_obs_filter_(0),
      Tmax_(1.0), t_clock_(0.0), sflags_(sflags), mp_(T()), opt_(nullptr),
      filter_(player::init_ball_filter(0.3, 0.001, false)) {

  std::string home = std::getenv("HOME");
  std::string file = home + "/projects/table-tennis/json/" + sflags_.json_file;
  mp_ = T(file);

  if (sflags_.start_from_act_state) {
    mp_.set_init_pos(qinit);
  }

  double Tmax_ = 2.0;
  double lb[2 * NDOF + 1], ub[2 * NDOF + 1];
  set_bounds(lb, ub, 0.01, Tmax_);
  mp_.get_init_pos(q_rest_des_);
  opt_ = new FocusedOptim(q_rest_des_, lb, ub);
  timer_.tic();
}

template <typename T> ServeBall<T>::~ServeBall() { delete opt_; }

template <typename T> void ServeBall<T>::set_flags(const serve_flags &sflags) {
  sflags_ = sflags;
}

template <typename T>
void ServeBall<T>::serve(const ball_obs &obs, const joint &qact, joint &qdes) {

  estimate_ball_state(obs);

  if (opt_->get_params(qact, p_)) {
    t_clock_ = 0.0;
    ran_optim_ = true;
    opt_->set_moving(true);
  }

  if (ran_optim_) {
    update_next_state(p_, q_rest_des_, 1.0, t_clock_, qdes);
  } else {
    double dt = sflags_.speedup * DT;
    mp_.step(dt, qdes);
    t_clock_ += dt;
  }

  if (sflags_.mpc && (!ran_optim_ || timer_.toc() > (1.0 / sflags_.freq_mpc))) {
    double t_pred = (Tmax_ - t_clock_) / DT;
    if (!predict_ball_hit(t_pred)) {
      // launch optimizer
      cout << "Ball won't be hit! Launching optimization to correct traj...\n";
      correct_with_optim(qact);
      timer_.tic();
    }
  } else {
    // cout << "Predicting ball hit..." << endl;
  }
}

template <typename T> void ServeBall<T>::correct_with_optim(const joint &qact) {

  double Tpred = 2.0;
  mat balls_pred = zeros<mat>(6, Tpred / DT);
  predict_ball(Tpred, balls_pred, filter_);

  // lookup_soln(filter_.get_mean(),1,qact);
  double time_land_des = sflags_.time_land_des;
  vec2 ball_land_des;
  ball_land_des(X) = sflags_.ball_land_des_offset[0];
  ball_land_des(Y) =
      sflags_.ball_land_des_offset[1] + dist_to_table - 1 * table_length / 4;
  optim_des pred_params;
  pred_params.Nmax = 1000;
  calc_racket_strategy(balls_pred, ball_land_des, time_land_des, pred_params);
  FocusedOptim *fp = static_cast<FocusedOptim *>(opt_);
  fp->set_return_time(1.0);
  fp->set_detach(sflags_.detach);
  fp->set_des_params(&pred_params);
  fp->set_verbose(sflags_.verbose);
  fp->update_init_state(qact);
  fp->run();
}

template <typename T>
void ServeBall<T>::estimate_ball_state(const ball_obs &obs) {

  const int min_ball_to_start_filter = 5;
  static mat init_obs = zeros<mat>(3, min_ball_to_start_filter);
  static vec init_times = zeros<vec>(min_ball_to_start_filter);

  // FILTERING
  if (!init_filter_) {
    if (idx_balls_obs_filter_ == min_ball_to_start_filter && obs.status) {
      vec6 init_filter_state;
      estimate_ball_linear(init_obs, init_times, false, init_filter_state);
      mat P;
      P.eye(6, 6);
      filter_.set_prior(init_filter_state, P);
    } else if (idx_balls_obs_filter_ < min_ball_to_start_filter && obs.status) {
      init_obs.col(idx_balls_obs_filter_) = obs.pos;
      init_times(idx_balls_obs_filter_) = idx_balls_obs_filter_ * DT;
      idx_balls_obs_filter_++;
    }
  } else {
    filter_.predict(DT, true);
    if (obs.status)
      filter_.update(obs.pos);
  }
}

template <typename T>
bool ServeBall<T>::predict_ball_hit(const double &t_pred) {

  // PREDICT IF BALL WILL BE HIT, OTHERWISE LAUNCH TRAJ OPTIM
  try {
    TableTennis tt_pred = TableTennis(filter_.get_mean(), false, false);
    racket robot_pred_racket;
    optim::joint Q_pred;
    unsigned N_pred = t_pred / DT;

    if (!ran_optim_) {
      T pred_mp = mp_;
      for (unsigned j = 0; j < N_pred; j++) {
        pred_mp.step(DT, Q_pred);
        calc_racket_state(Q_pred, robot_pred_racket);
        tt_pred.integrate_ball_state(robot_pred_racket, DT);
      }
    } else {
      double t_poly_pred = t_pred;
      for (unsigned j = 0; j < N_pred; j++) {
        update_next_state(p_, q_rest_des_, 1.0, t_poly_pred, Q_pred);
        calc_racket_state(Q_pred, robot_pred_racket);
        tt_pred.integrate_ball_state(robot_pred_racket, DT);
      }
    }
    return tt_pred.was_legally_served();
  } catch (const std::exception &not_init_error) {
    // filter_ not yet initialized
    return true;
  }
}

vec_str get_files(std::string folder_name, std::string prefix) {

  using std::string;
  vec_str files;
  FILE *stream;
  const int max_buffer = 256;
  char buffer[max_buffer];
  string cmd = "ls " + folder_name;
  // cmd.append(" 2>&1");

  stream = popen(cmd.c_str(), "r");
  if (stream) {
    while (!feof(stream)) {
      if (fgets(buffer, max_buffer, stream) != NULL) {
        std::string str(buffer);
        if (str.find(prefix) == 0)
          files.push_back(str);
      }
    }
    pclose(stream);
  }
  // remove newline from each file
  for (std::string &file : files) {
    file.erase(std::remove(file.begin(), file.end(), '\n'), file.end());
  }
  return files;
}

template class ServeBall<dmps>;
template class ServeBall<RBF>;
template class ServeBall<CRBF>;

} // namespace serve
