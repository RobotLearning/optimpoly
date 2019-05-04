/**
 * @file player.cpp
 *
 * @brief Table Tennis player class and the functions it uses are stored here.
 *
 * Player class launches 3 different optimization algorithms
 * for striking trajectory generation.
 *
 *  Created on: Feb 8, 2017
 *  Author: okoc
 */

#include "player.hpp"
#include "stdlib.h"
#include <armadillo>
#include <string>
#include <thread>

#include "ball_interface.h"
#include "constants.h"
#include "kalman.h"
#include "kinematics.h"
#include "kinematics.hpp"
#include "lookup.h"
#include "optim.h"
#include "tabletennis.h"

using namespace arma;
using namespace optim;

namespace player {

Player::Player(const Player &player)
    : filter_(player.filter_), pflags_(player.pflags_),
      ball_land_des_(player.ball_land_des_), q_rest_des_(player.q_rest_des_) {
  observations_ = zeros<mat>(3, pflags_.min_obs);
  times_ = zeros<vec>(pflags_.min_obs);
  set_optim();
}

Player::Player(const vec7 &q0, player_flags &flags) : pflags_(flags) {

  ball_land_des_(X) += pflags_.ball_land_des_offset[X];
  ball_land_des_(Y) =
      dist_to_table - 3 * table_length / 4 + pflags_.ball_land_des_offset[Y];
  q_rest_des_ = q0;
  observations_ = zeros<mat>(3, pflags_.min_obs); // for initializing filter_
  times_ = zeros<vec>(pflags_.min_obs);           // for initializing filter_
  filter_ = init_ball_filter(0.3, 0.001, pflags_.spin_based_pred);
  // load_lookup_table(lookup_table_);
  set_optim();
}

Player::~Player() { delete opt_; }

void Player::set_optim() {

  double lb[2 * NDOF + 1];
  double ub[2 * NDOF + 1];
  double SLACK = 0.02;
  double Tmax = 1.0;
  set_bounds(lb, ub, SLACK, Tmax);

  switch (pflags_.alg) {
  case FOCUS: {
    opt_ = new FocusedOptim(q_rest_des_, lb, ub);
    pred_params_.Nmax = 1000;
    break;
  }
  case VHP: {
    opt_ = new HittingPlane(q_rest_des_, lb, ub);
    break;
  }
  case DP: {
    opt_ = new DefensiveOptim(q_rest_des_, lb, ub, true, true); // use lookup
    DefensiveOptim *dp = static_cast<DefensiveOptim *>(opt_);
    dp->set_weights(pflags_.weights);
    dp->set_velocity_multipliers(pflags_.mult_vel);
    dp->set_penalty_loc(pflags_.penalty_loc);
    pred_params_.Nmax = 1000;
    break;
  }
  default:
    throw("Algorithm is not recognized!\n");
  }
  opt_->set_return_time(pflags_.time2return);
  opt_->set_verbose(pflags_.verbosity > 1);
  opt_->set_detach(pflags_.detach);
}

Player &Player::operator=(const Player &player) {

  if (this != &player) {
    filter_ = player.filter_;
    pflags_ = player.pflags_;
    ball_land_des_ = player.ball_land_des_;
    q_rest_des_ = player.q_rest_des_;
    observations_ = zeros<mat>(3, pflags_.min_obs);
    times_ = zeros<vec>(pflags_.min_obs);
    set_optim();
  }
  return *this;
}

bool Player::filter_is_initialized() const { return init_ball_state_; }

void Player::estimate_ball_state(const ball_obs &obs, const double &dt) {

  using std::ref;
  using std::thread;
  int verb = pflags_.verbosity;
  valid_obs_ = false;

  if (check_reset_filter(obs.status, verb, pflags_.t_reset_thresh)) {
    filter_ = init_ball_filter(pflags_.var_model, pflags_.var_noise, pflags_.spin_based_pred,
                              pflags_.out_reject_mult);
    num_obs_ = 0;
    init_ball_state_ = false;
    game_state_ = AWAITING;
    t_obs_ = 0.0; // t_cumulative
  }

  if (num_obs_ < pflags_.min_obs && obs.status) {
    times_(num_obs_) = t_obs_;
    observations_.col(num_obs_) = obs.pos;
    num_obs_++;
    if (num_obs_ == pflags_.min_obs) {
      if (verb >= 1)
        cout << "Estimating initial ball state\n";
      thread t =
          thread(estimate_prior, ref(observations_), ref(times_),
                 ref(pflags_.verbosity), ref(init_ball_state_), ref(filter_));
      if (pflags_.detach) {
        t.detach();
      }
      else {
        t.join();
      }
    }

  } else if (init_ball_state_) { // comes here if there are enough balls to start
                                // filter_
    filter_.predict(dt, true);   // true);
    if (obs.status) {
      valid_obs_ = true;
      if (pflags_.outlier_detection)
        valid_obs_ = !filter_.check_outlier(obs.pos, verb > 2);
    }
    if (valid_obs_) {
      filter_.update(obs.pos);
    }
  }
  t_obs_ += dt;
}

vec6 Player::filt_ball_state(const vec3 &obs) {

  ball_obs obs_str;
  obs_str.status = true;
  obs_str.pos = obs;
  estimate_ball_state(obs_str);
  try {
    return filter_.get_mean();
  } catch (const std::exception &exception) {
    return join_vert(obs, zeros<vec>(3));
  }
}

void Player::play(const ball_obs &obs, const joint &qact, joint &qdes) {

  estimate_ball_state(obs, const_tt::DT);

  switch (pflags_.alg) {
  case FOCUS:
    optim_fp_param(qact);
    break;
  case VHP:
    optim_vhp_param(qact);
    break;
  case DP:
    optim_dp_param(qact);
    break;
  default:
    throw("Algorithm is not recognized!\n");
  }

  // generate movement or calculate next desired step
  calc_next_state(qact, qdes);
}

void Player::cheat(const joint &qact, const vec6 &ballstate, joint &qdes) {

  // resetting legal ball detecting to AWAITING state
  if (ballstate(Y) < (dist_to_table - table_length) && ballstate(DY) > 2.0)
    game_state_ = AWAITING;
  filter_.set_prior(ballstate, 0.01 * eye<mat>(6, 6));

  switch (pflags_.alg) {
  case FOCUS:
    optim_fp_param(qact);
    break;
  case VHP:
    optim_vhp_param(qact);
    break;
  case DP:
    optim_dp_param(qact);
    break;
  default:
    throw("Algorithm is not recognized!\n");
  }

  // generate movement or calculate next desired step
  calc_next_state(qact, qdes);
}

void Player::optim_vhp_param(const joint &qact) {

  vec6 ball_pred;
  double time_pred;

  // if ball is fast enough and robot is not moving consider optimization
  if (check_update(qact)) {
    if (predict_hitting_point(pflags_.VHPY, pflags_.check_bounce, ball_pred,
                              time_pred, filter_,
                              game_state_)) { // ball is legal and reaches VHP
      calc_racket_strategy(ball_pred, ball_land_des_, pflags_.time_land_des,
                           pred_params_);
      HittingPlane *vhp = static_cast<HittingPlane *>(opt_);
      vhp->set_des_params(&pred_params_);
      vhp->fix_hitting_time(time_pred);
      vhp->update_init_state(qact);
      vhp->run();
    }
  }
}

void Player::optim_fp_param(const joint &qact) {

  mat balls_pred;

  // if ball is fast enough and robot is not moving consider optimization
  if (check_update(qact)) {
    predict_ball(2.0, balls_pred, filter_);
    if (!pflags_.check_bounce || check_legal_ball(filter_.get_mean(), balls_pred,
                                                 game_state_)) { // ball is legal
      // lookup_soln(filter_.get_mean(),1,qact);
      calc_racket_strategy(balls_pred, ball_land_des_, pflags_.time_land_des,
                           pred_params_);
      FocusedOptim *fp = static_cast<FocusedOptim *>(opt_);
      fp->set_des_params(&pred_params_);
      fp->update_init_state(qact);
      fp->run();
    } else {
      // cout << "Ball is not legal!\n";
    }
  }
}

void Player::optim_dp_param(const joint &qact) {

  mat balls_pred;

  // if ball is fast enough and robot is not moving consider optimization
  if (check_update(qact)) {
    predict_ball(2.0, balls_pred, filter_);
    if (!pflags_.check_bounce || check_legal_ball(filter_.get_mean(), balls_pred,
                                                 game_state_)) { // ball is legal
      // lookup_soln(filter_.get_mean(),1,qact);
      // calc_racket_strategy(balls_pred,ball_land_des_,time_land_des,pred_params_);
      pred_params_.ball_pos = balls_pred.rows(X, Z);
      pred_params_.ball_vel = balls_pred.rows(DX, DZ);
      pred_params_.Nmax = balls_pred.n_cols;
      DefensiveOptim *dp = static_cast<DefensiveOptim *>(opt_);
      dp->set_des_params(&pred_params_);
      dp->update_init_state(qact);
      dp->run();
    }
  }
}

bool Player::check_update(const joint &qact) const {

  static wall_clock timer;
  vec6 state_est;
  bool update = false;
  racket robot_racket;
  bool activate, passed_lim, feasible = false;

  if (!init_ball_state_) {
    timer.tic();
  }

  try {
    state_est = filter_.get_mean();
    feasible =
        (state_est(DY) > 0.5) &&
        (state_est(Y) > (dist_to_table - table_length + pflags_.optim_offset));
    update = !opt_->check_update() && !opt_->check_running();

    // ball is incoming
    if (pflags_.mpc) { // && t_poly_ > 0.0) {
      calc_racket_state(qact, robot_racket);
      activate = timer.toc() > (1.0 / pflags_.freq_mpc);
      passed_lim = state_est(Y) > robot_racket.pos(Y);
      update = update && valid_obs_ && activate && feasible && !passed_lim;
    } else {
      update = update && (t_poly_ == 0.0) && feasible; // only once
    }
    if (update) {
      // cout << "Consider reacting to ball: " << state_est.t() << endl;
      timer.tic();
    }
  } catch (const std::exception &not_init_error) {
    update = false;
  }

  return update;
}

void Player::calc_next_state(const joint &qact, joint &qdes) {

  using std::ref;
  using std::thread;
  // this should be only for MPC?
  if (opt_->get_params(qact, poly_)) {
    if (pflags_.verbosity) {
      std::cout << "Launching/updating strike" << std::endl;
    }
    t_poly_ = const_tt::DT;
    opt_->set_moving(true);
  }

  // make sure we update after optim finished
  if (t_poly_ > 0.0) {
    if (!update_next_state(poly_, q_rest_des_, pflags_.time2return, t_poly_,
                           qdes)) {
      opt_->set_moving(false);
      // optimize to find a better resting state close to predicted balls
      if (pflags_.optim_rest_posture) {
        opt_->run_qrest_optim(q_rest_des_);
      }
    }
  }
}

void Player::reset_filter(double var_model, double var_noise) {

  filter_ = init_ball_filter(var_model, var_noise, pflags_.spin_based_pred,
                            pflags_.out_reject_mult);
  init_ball_state_ = false;
  num_obs_ = 0;
  game_state_ = AWAITING;
  t_obs_ = 0.0;
}

void Player::lookup_soln(const vec6 &ball_state, const int k,
                         const joint &qact) {

  double time2return = pflags_.time2return;
  if (t_poly_ == 0.0) {
    if (pflags_.verbosity) {
      cout << "Starting movement based on lookup, k = 5\n"; // kNN parameter k =
                                                            // 5
    }
    vec::fixed<15> robot_params;
    vec6 ball_est = ball_state;
    // cout << "Init ball est:" << ball_params << endl;
    predict_till_net(ball_est);
    // cout << "Net ball est:" << ball_params << endl;
    knn(lookup_table_, ball_est, k, robot_params);
    vec7 qf, qfdot;
    for (int i = 0; i < NDOF; i++) {
      qf(i) = robot_params(i);
      qfdot(i) = robot_params(i + NDOF);
    }
    double T = robot_params(2 * NDOF);
    vec7 qnow = qact.q;
    vec7 qdnow = qact.qd;
    poly_.a.col(0) = 2.0 * (qnow - qf) / pow(T, 3) + (qfdot + qdnow) / pow(T, 2);
    poly_.a.col(1) = 3.0 * (qf - qnow) / pow(T, 2) - (qfdot + 2.0 * qdnow) / T;
    poly_.a.col(2) = qdnow;
    poly_.a.col(3) = qnow;
    // cout << "A = \n" << p.a << endl;
    poly_.b.col(0) = 2.0 * (qf - q_rest_des_) / pow(time2return, 3) +
                    (qfdot) / pow(time2return, 2);
    poly_.b.col(1) = 3.0 * (q_rest_des_ - qf) / pow(time2return, 2) -
                    (2.0 * qfdot) / time2return;
    poly_.b.col(2) = qfdot;
    poly_.b.col(3) = qf;
    poly_.time2hit = T;
    t_poly_ = const_tt::DT;
  }
}

bool predict_hitting_point(const double &vhpy, const bool &check_bounce,
                           vec6 &ball_pred, double &time_pred, EKF &filter_,
                           game &game_state_) {

  const double time_min = 0.05;
  mat balls_path;
  bool valid_hp = false;
  predict_ball(2.0, balls_path, filter_);
  uvec vhp_index;
  unsigned idx;

  if (!check_bounce || check_legal_ball(filter_.get_mean(), balls_path,
                                        game_state_)) { // ball is legal
    vhp_index = find(balls_path.row(Y) >= vhpy, 1);
    if (vhp_index.n_elem == 1) {
      idx = as_scalar(vhp_index);
      ball_pred = balls_path.col(idx);
      time_pred = const_tt::DT * (idx + 1);
      if (time_pred > time_min)
        valid_hp = true;
    }
  }

  return valid_hp;
}

void predict_ball(const double &time_pred, mat &balls_pred, const EKF &filter_) {

  // static wall_clock timer;
  // timer.tic();
  int N = (int)(time_pred / const_tt::DT);
  balls_pred = filter_.predict_path(const_tt::DT, N);
  // cout << "Pred. ball time: " << 1000 * timer.toc() << " ms." << endl;
}

void check_legal_bounce(const vec6 &ball_est, game &game_state_) {

  static double last_y_pos = 0.0;
  static double last_z_vel = 0.0;
  bool incoming = (ball_est(DY) > 0.0);
  bool on_opp_court = (ball_est(Y) < (dist_to_table - (table_length / 2.0)));
  bool bounce = (last_z_vel < 0.0 && ball_est(DZ) > 0.0) &&
                (fabs(ball_est(Y) - last_y_pos) < 0.1);

  if (bounce && incoming) {
    // incoming ball has bounced
    if (game_state_ == LEGAL) {
      cout << "Ball bounced twice\n";
      game_state_ = ILLEGAL;
    } else if (game_state_ == AWAITING && on_opp_court) {
      cout << "Ball bounced on opponents court!\n";
      game_state_ = ILLEGAL;
    } else {
      cout << "Legal bounce occurred!" << endl;
      game_state_ = LEGAL;
    }
  }
  last_y_pos = ball_est(Y);
  last_z_vel = ball_est(DZ);
}

bool check_legal_ball(const vec6 &ball_est, const mat &balls_predicted,
                      game &game_state_) {

  int num_bounces = 0;
  int N = balls_predicted.n_cols;

  check_legal_bounce(ball_est, game_state_);

  // if sign of z-velocity changes then the ball bounces
  for (int i = 0; i < N - 1; i++) {
    if (balls_predicted(DZ, i) < 0.0 && balls_predicted(DZ, i + 1) > 0.0) {
      num_bounces++;
    }
  }

  // one bounce is predicted before an actual bounce has happened
  if (game_state_ == AWAITING && num_bounces == 1) {
    return true;
  }
  // no bounce is predicted
  if (game_state_ == LEGAL && num_bounces == 0) {
    return true;
  }

  return false;
}

void Player::get_strategy(vec2 &pos_land_des, double &time_land_des) {

  switch (pflags_.alg) {
  case FOCUS:
    pos_land_des = this->ball_land_des_;
    time_land_des = this->pflags_.time_land_des;
    break;
  case VHP:
    pos_land_des = this->ball_land_des_;
    time_land_des = this->pflags_.time_land_des;
    break;
  case DP:
    std::cout << "DP does not have a fixed return strategy!\n";
    pos_land_des = zeros<vec>(2);
    time_land_des = 0.0;
    break;
  default:
    throw("Algorithm is not recognized!\n");
  }
}

} // namespace player
