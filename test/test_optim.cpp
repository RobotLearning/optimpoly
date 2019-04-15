/*
 * test_optim.cpp
 *
 * Unit Tests for polynomial optimization
 *
 *  Created on: Feb 17, 2017
 *      Author: okoc
 */
#include "kinematics.h"
#include "kinematics.hpp"
#include "lookup.h"
#include "optim.h"
#include "player.hpp"
#include "tabletennis.h"
#include "utils.h"
#include "gtest/gtest.h"
#include <armadillo>
#include <iostream>
#include <thread>

using namespace arma;
using namespace const_tt;
using namespace optim;
using namespace player;

struct des_ball_data { // Initial ball data time stamps and position
                       // observations
  vec ball_incoming = zeros<vec>(3);
  vec ball_land_des = zeros<vec>(3);
  double time_land_des = 0.5;
  double topspin = 0.0;
};

static double calc_landing_res(unsigned n, const double *x, double *grad,
                               void *data);
static void optim_spin_outgoing_ball_vel(const des_ball_data &data,
                                         const bool verbose,
                                         vec3 &est); // spin based optimization
static void init_right_posture(vec7 &q0);

TEST(OptimTests, TestVirtualHittingPlaneOptim) {

  const double VHPY = -0.3;
  static double lb[2 * NDOF + 1];
  static double ub[2 * NDOF + 1];
  double SLACK = 0.01;
  double Tmax = 1.0;
  joint qact;
  spline_params poly;

  // update initial parameters from lookup table
  arma_rng::set_seed(1); // looking for a random ball entry
  // arma_rng::set_seed_random();
  vec strike_params = zeros<vec>(15);
  vec ball_state = zeros<vec>(6);
  lookup_random_entry(ball_state, strike_params);
  init_right_posture(qact.q);
  set_bounds(lb, ub, SLACK, Tmax);

  EKF filter = init_ball_filter();
  mat66 P;
  P.eye();
  filter.set_prior(ball_state, P);

  double time_pred = 0.0;
  vec6 ball_pred;
  game game_state = AWAITING;
  vec2 ball_land_des = {0.0, dist_to_table - 3 * table_length / 4};
  double time_land_des = 0.8;
  EXPECT_TRUE(predict_hitting_point(VHPY, true, ball_pred, time_pred, filter,
                                    game_state));
  // cout << ball_pred << endl;
  optim_des racket_params;
  calc_racket_strategy(ball_pred, ball_land_des, time_land_des, racket_params);

  vec3 normal_example = racket_params.racket_normal(span(X, Z), 0);
  EXPECT_EQ(arma::norm(normal_example), 1.0); // tol = 0.01

  HittingPlane opt = HittingPlane(qact.q.memptr(), lb, ub);
  opt.set_des_params(&racket_params);
  opt.fix_hitting_time(time_pred);
  opt.update_init_state(qact);
  opt.run();
  bool update = opt.get_params(qact, poly);

  EXPECT_TRUE(update);
}

TEST(OptimTests, TestFocusedPlayerOptim) {

  double lb[2 * NDOF + 1], ub[2 * NDOF + 1];
  double SLACK = 0.01;
  double Tmax = 1.0;
  joint qact;
  spline_params poly;

  // update initial parameters from lookup table
  arma_rng::set_seed(1); // looking up a random ball entry
  // arma_rng::set_seed_random();
  vec::fixed<15> strike_params;
  vec6 ball_state;
  lookup_random_entry(ball_state, strike_params);
  init_right_posture(qact.q);
  set_bounds(lb, ub, SLACK, Tmax);
  optim_des racket_params;
  int N = 1000;
  racket_params.Nmax = 1000;

  EKF filter = init_ball_filter();
  mat66 P;
  P.eye();
  filter.set_prior(ball_state, P);

  vec6 ball_pred;
  double time_land_des = 0.8;
  mat balls_pred = filter.predict_path(DT, N);
  vec2 ball_land_des = {0.0, dist_to_table - 3 * table_length / 4};
  racket_params = calc_racket_strategy(balls_pred, ball_land_des, time_land_des,
                                       racket_params);

  FocusedOptim opt = FocusedOptim(qact.q.memptr(), lb, ub);
  opt.set_des_params(&racket_params);
  opt.update_init_state(qact);
  opt.run();
  bool update = opt.get_params(qact, poly);

  EXPECT_TRUE(update);
}

TEST(OptimTests, TestDefensivePlayerOptim) {

  double lb[2 * NDOF + 1], ub[2 * NDOF + 1];
  double SLACK = 0.01;
  double Tmax = 1.0;
  joint qact;
  spline_params poly;

  // update initial parameters from lookup table
  arma_rng::set_seed(1); // looking up a random ball entry
  vec::fixed<15> strike_params;
  vec6 ball_state;
  lookup_random_entry(ball_state, strike_params);
  init_right_posture(qact.q);
  set_bounds(lb, ub, SLACK, Tmax);

  int N = 1000;
  EKF filter = init_ball_filter();
  mat66 P;
  P.eye();
  filter.set_prior(ball_state, P);
  mat balls_pred = filter.predict_path(DT, N);
  optim_des ball_params;
  ball_params.ball_pos = balls_pred.rows(X, Z);
  ball_params.ball_vel = balls_pred.rows(DX, DZ);
  ball_params.Nmax = N;
  bool land = true;
  bool lookup = true;
  DefensiveOptim opt = DefensiveOptim(qact.q.memptr(), lb, ub, land,
                                      lookup); // only touch the ball if false!
  opt.set_des_params(&ball_params);
  opt.update_init_state(qact);
  opt.set_verbose(true);
  opt.run();
  bool update = opt.get_params(qact, poly);

  EXPECT_TRUE(update);
}

TEST(OptimTests, CheckRestingPostureOptim) {

  // Find a qf and t such that J(qf) has minimal Frobenius norm
  // while being close to predicted ball state b(t) and being close to q_hit
  double lb[2 * NDOF + 1], ub[2 * NDOF + 1];
  double SLACK = 0.01;
  double Tmax = 1.0;
  vec7 q_rest_des = zeros<vec>(7);
  joint qact;
  // update initial parameters from lookup table
  arma_rng::set_seed_random(); // looking up a random ball entry
  vec::fixed<15> strike_params;
  vec6 ball_state;
  lookup_random_entry(ball_state, strike_params);
  vec7 q_hit = strike_params.head(NDOF);

  set_bounds(lb, ub, SLACK, Tmax);
  int N = 1000;
  EKF filter = init_ball_filter();
  mat66 P;
  P.eye();
  filter.set_prior(ball_state, P);
  mat balls_pred = filter.predict_path(DT, N);
  optim_des ball_params;
  ball_params.ball_pos = balls_pred.rows(X, Z);
  ball_params.ball_vel = balls_pred.rows(DX, DZ);
  ball_params.Nmax = N;
  DefensiveOptim opt = DefensiveOptim(q_rest_des, lb, ub, true, true);
  opt.set_verbose(false);
  opt.set_des_params(&ball_params);
  opt.update_init_state(qact);
  opt.run();
  opt.run_qrest_optim(q_rest_des);
  // nothing to test here
}

TEST(OptimTests, CheckAccuracyOfRacketCalculationsBVPUsingSpinBallModel) {

  // For the Focused Player with a spin model and spin estimation:
  // comparing with the ballistic model

  // Solving BVP for one particular ball...
  // Solving for desired outgoing velocity using a spin model");
  // Optimization should produce accurate inversion...
  arma_rng::set_seed_random();
  double topspin = -50.0;
  TableTennis tt = TableTennis(true, false);
  tt.set_topspin(topspin);
  int ball_launch_side = (randi(1, distr_param(0, 2)).at(0));
  tt.set_ball_gun(0.05, ball_launch_side);
  int N = 50;
  double dt = 0.02;
  // mat balls_pred = zeros<mat>(6,N);
  for (int i = 0; i < N; i++) {
    tt.integrate_ball_state(dt);
    // balls_pred.col(i) = tt.get_ball_state();
  }
  vec3 ball_land_des = {0.0, dist_to_table - 3 * table_length / 4.0,
                        floor_level - table_height + ball_radius};
  double time_land_des = 0.6;
  vec3 ball_vel_out; // estimated outgoing ball velocity with ballistic model
  vec6 ballin = tt.get_ball_state();
  tt.calc_des_ball_out_vel(ball_land_des.head(2), time_land_des, true, ballin,
                           ball_vel_out); // Incoming ball: ballin,

  int N_horizon = time_land_des / dt;
  tt.set_ball_state(join_vert(ballin.head(3), ball_vel_out));
  for (int i = 0; i < N_horizon; i++) {
    tt.integrate_ball_state(dt);
  }
  // Resulting landing position: tt.get_ball_position().t()
  double err_norm_no_spin =
      norm(tt.get_ball_position() - ball_land_des, 2); // error norm
  // initialize optimization with adjusted values (e.g. to reduce optim time)
  ball_vel_out(X) *= 1.1;
  ball_vel_out(X) *= 1.2;
  ball_vel_out(X) *= 0.9;

  des_ball_data data;
  data.ball_land_des = ball_land_des;
  data.time_land_des = time_land_des;
  data.topspin = topspin;
  data.ball_incoming = ballin.head(3);
  optim_spin_outgoing_ball_vel(data, true,
                               ball_vel_out); // soln. to BVP with spin

  tt.set_ball_state(join_vert(ballin.head(3), ball_vel_out));
  for (int i = 0; i < N_horizon; i++) {
    tt.integrate_ball_state(dt);
  } // resulting land. pos: tt.get_ball_position()
  double err_norm_spin =
      norm(tt.get_ball_position() - ball_land_des, 2); // norm of error

  EXPECT_LT(err_norm_spin, err_norm_no_spin);
}

/**
 * @brief Solve BVP for a particular predicted spinning ball's outgoing desired
 * ball velocity
 *
 * BVP is solved using optimization
 */
static void optim_spin_outgoing_ball_vel(const des_ball_data &data,
                                         const bool verbose, vec3 &est) {

  static double x[3]; /* some initial guess */
  static double minf; /* the minimum objective value, upon return */
  static double init_time;
  static int res; // error code
  static nlopt_opt opt;

  opt = nlopt_create(NLOPT_LD_MMA, 3);
  nlopt_set_min_objective(opt, calc_landing_res, (void *)&data);
  nlopt_set_xtol_rel(opt, 1e-2);

  for (int i = 0; i < 3; i++) {
    x[i] = est(i);
  }

  init_time = get_time();
  if ((res = nlopt_optimize(opt, x, &minf)) < 0) {
    if (verbose)
      printf("NLOPT failed!\n");
  } else {
    if (verbose) {
      printf("NLOPT took %f ms\n", (get_time() - init_time) / 1e3);
      printf("Found minimum at f = %0.10g\n", minf);
    }
    for (int i = 0; i < 3; i++) {
      est(i) = x[i];
    }
  }
  // nlopt_destroy(opt);
}

/**
 * Cost function for computing the residual (norm squared)
 * of the outgoing ball landing error
 * Calculates also the gradient if grad is TRUE
 *
 */
static double calc_landing_res(unsigned n, const double *x, double *grad,
                               void *data) {

  static double dt = 0.02;
  static TableTennis tt =
      TableTennis(true, false, false); // no contact checking
  static vec3 vel_out;
  static vec6 init_state;
  static vec3 out_pos;

  des_ball_data *mydata = (des_ball_data *)data;
  tt.set_topspin(mydata->topspin);
  vel_out(X) = x[0];
  vel_out(Y) = x[1];
  vel_out(Z) = x[2];
  assert(n == 3);
  init_state = join_vert(mydata->ball_incoming, vel_out);
  tt.set_ball_state(init_state);
  for (int i = 0; i < mydata->time_land_des / dt; i++)
    tt.integrate_ball_state(dt);
  // for (int i = 0; i < 5; i++)
  //  tt.symplectic_int_fourth(mydata->time_land_des/5.0);

  out_pos = tt.get_ball_position();

  if (grad) {

    grad[0] = mydata->time_land_des * (out_pos(X) - mydata->ball_land_des(X));
    grad[1] = mydata->time_land_des * (out_pos(Y) - mydata->ball_land_des(Y));
    grad[2] = mydata->time_land_des * (out_pos(Z) - mydata->ball_land_des(Z));
    // Finite difference
    /*static double h = 1e-6;
    static double val_plus, val_minus;
    static double xx[3];
    for (unsigned i = 0; i < n; i++)
        xx[i] = x[i];
    for (unsigned i = 0; i < n; i++) {
        xx[i] += h;
        val_plus = calc_landing_res(n, xx, NULL, data);
        xx[i] -= 2*h;
        val_minus = calc_landing_res(n, xx, NULL, data);
        grad[i] = (val_plus - val_minus) / (2*h);
        xx[i] += h;
    }*/
  }

  return pow(norm(out_pos - mydata->ball_land_des), 2);
}

/*
 * Initialize robot posture on the right size of the robot
 */
static void init_right_posture(vec7 &q0) {

  q0(0) = 1.0;
  q0(1) = -0.2;
  q0(2) = -0.1;
  q0(3) = 1.8;
  q0(4) = -1.57;
  q0(5) = 0.1;
  q0(6) = 0.3;
}
