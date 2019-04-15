#include "constants.h"
#include "dmp.h"
#include "kinematics.hpp"
#include "optim.h"
#include "rbf.h"
#include "serve.h"
#include "tabletennis.h"
#include "gtest/gtest.h"
#include <armadillo>
#include <iostream>
#include <thread>

using namespace arma;
using namespace const_tt;
using namespace player;
using namespace serve;
using dmps = Joint_DMPs;
using vec_str = std::vector<std::string>;

static void init_posture(vec7 &q0, int posture, bool verbose);

vec6 init_ball_vertical(const double &T, const double &std_noise,
                        const vec7 &goal_state);

TEST(ServeTests, CheckLegalLandingForServeStartingWithDMPAndSwitchingToOptim) {

  // start evolving dmp
  // have a ball coming down to goal state
  // track the ball with a filter
  // if we predict miss then run an optimizer to correct
  // switch to optimizer

  const double T = 1.0;
  std::string dmp_file;
  dmps multi_dmp = init_dmps(dmp_file);
  optim::joint qdes, qact;
  TableTennis tt = TableTennis(false, true);
  double std_init_noise = 0.0;
  vec7 goal_state;
  multi_dmp.get_goal_pos(goal_state);
  vec6 init_ball_state = init_ball_vertical(T, std_init_noise, goal_state);
  tt.set_ball_state(init_ball_state);
  racket racket_state;
  serve_flags flags;
  flags.json_file = dmp_file;
  ServeBall<dmps> server = ServeBall<dmps>(flags);
  // sflags flags;
  // server.set_flags(flags);

  // check interaction and make sure ball is served correctly
  while (!tt.touched_ground() && !tt.was_legally_served()) {

    tt.integrate_ball_state(racket_state, DT);
    ball_obs obs;
    obs.status = true;
    obs.pos = tt.get_ball_position();
    server.serve(obs, qact, qdes);
    qact = qdes;
    calc_racket_state(qact, racket_state);
  }
  EXPECT_TRUE(tt.was_legally_served());
}

TEST(ServeTests,
     CheckLegalLandingForServeStartingWithRBFBasedDemoAndSwitchingToOpt) {
  const double T = 1.0;
  optim::joint qdes, qact;
  TableTennis tt = TableTennis(false, true);
  double std_init_noise = 0.0;
  vec7 goal_state;
  const std::string home = std::getenv("HOME");
  const std::string date = "16.11.18";
  std::string rbf_file = "rbf_3_" + date + ".json";
  std::string file_path = home + "/projects/table-tennis/json/" + rbf_file;
  CRBF rbf = CRBF(file_path);
  rbf.get_goal_pos(goal_state);
  vec6 init_ball_state = init_ball_vertical(T, std_init_noise, goal_state);
  tt.set_ball_state(init_ball_state);
  racket racket_state;
  serve_flags flags;
  flags.ball_land_des_y_offset -= 0.2;
  flags.json_file = rbf_file;
  ServeBall<CRBF> server = ServeBall<CRBF>(flags);

  // check interaction and make sure ball is served correctly
  while (!tt.touched_ground() && !tt.was_legally_served()) {
    tt.integrate_ball_state(racket_state, DT);
    ball_obs obs;
    obs.status = true;
    obs.pos = tt.get_ball_position();
    server.serve(obs, qact, qdes);
    qact = qdes;
    calc_racket_state(qact, racket_state);
  }
  EXPECT_TRUE(tt.was_legally_served());
}

/*
 * Initialize ball vertically such that assuming only gravity acts on the
 * ball (i.e. no airdrag) the ball will be around the MP goal cartesian position
 * T seconds from the start. [around = some noise added]
 */
vec6 init_ball_vertical(const double &T, const double &std_noise,
                        const vec7 &goal_state) {
  ball_params params;
  optim::joint goal_joint_state;
  goal_joint_state.q = goal_state;
  racket goal_racket_state;
  calc_racket_state(goal_joint_state, goal_racket_state);
  vec3 init_ball_pos = goal_racket_state.pos + std_noise * randn(3, 1);
  init_ball_pos(Z) -= 0.83 * params.gravity * T * T / 2.0; // 0.83 - dmp works!
  vec6 init_ball_state = zeros<vec>(6);
  init_ball_state.head(3) = init_ball_pos;

  return init_ball_state;
}

TEST(ServeTests, CheckDMPIsStableInReachingGoalPosition) {

  // create a dmp and evolve
  const double T = 1.0;
  const std::string home = std::getenv("HOME");
  const std::string date = "16-Nov-2018";
  const std::string file =
      home + "/projects/table-tennis/json/dmp_1_" + date + ".json";
  dmps multi_dmp = dmps(file);
  mat M = multi_dmp.evolve(T);

  // reinitialize in a different init. state and evolve
  vec7 q0;
  init_posture(q0, 2, 1);
  multi_dmp.set_init_pos(q0);
  // evolve it
  mat M2 = multi_dmp.evolve(T);

  // change goal state slightly and evolve
  vec goal = zeros<vec>(NDOF);
  multi_dmp.get_goal_pos(goal);
  vec goal_new = goal + 0.01 * randn(NDOF, 1);
  multi_dmp.set_goal_pos(goal_new);
  // evolve it
  mat M3 = multi_dmp.evolve(T);

  EXPECT_TRUE(approx_equal(goal, M.tail_cols(1), "absdiff", 0.05));
  EXPECT_TRUE(approx_equal(goal, M2.tail_cols(1), "absdiff", 0.05));
  EXPECT_TRUE(approx_equal(goal_new, M3.tail_cols(1), "absdiff", 0.05));
}

TEST(ServeTests, DISABLED_CheckDMPDoesNotExceedMaxAccelerations) {

  // for all the json files
  // check max acc with goal/init pos adjusted
  const double MAX_ACC_ALLOWED = 50.0;
  const std::string home = std::getenv("HOME");
  vec_str files = get_files(home + "/projects/table-tennis/json/", "dmp");
  dmps multi_dmp;

  for (std::string &file : files) {
    // file.erase(std::remove(file.begin(),file.end(),'\n'),file.end());
    std::string full_file = home + "/projects/table-tennis/json/" + file;
    // cout << full_file << endl;
    multi_dmp = dmps(full_file);
    // multi_dmp.turn_on_safe_acc();
    mat Q, Qd, Qdd;
    multi_dmp.evolve(1.0, Q, Qd, Qdd);
    double max_acc = max(max(abs(Qdd)));
    EXPECT_LT(max_acc, MAX_ACC_ALLOWED);

    vec goal = zeros<vec>(NDOF);
    multi_dmp.get_goal_pos(goal);
    vec goal_new = goal + 0.01 * randn(NDOF, 1);
    multi_dmp.set_goal_pos(goal_new);
    multi_dmp.evolve(1.0, Q, Qd, Qdd);
    max_acc = max(max(abs(Qdd)));
    EXPECT_LT(max_acc, MAX_ACC_ALLOWED);

    vec init = zeros<vec>(NDOF);
    multi_dmp.get_init_pos(init);
    vec init_new = init + 0.01 * randn(NDOF, 1);
    multi_dmp.set_init_pos(init_new);
    multi_dmp.evolve(1.0, Q, Qd, Qdd);
    max_acc = max(max(abs(Qdd)));
    EXPECT_LT(max_acc, MAX_ACC_ALLOWED);

    vec7 q0;
    init_posture(q0, 2, 0);
    multi_dmp.set_init_pos(q0);
    multi_dmp.evolve(1.0, Q, Qd, Qdd);
    max_acc = max(max(abs(Qdd)));
    EXPECT_LT(max_acc, MAX_ACC_ALLOWED);
  }
}

TEST(ServeTests, CompareSpedUpDMPToSubsampledNormalDMP) {

  const double T = 1.0;
  const std::string home = std::getenv("HOME");
  const std::string date = "16-Nov-2018";
  const std::string file =
      home + "/projects/table-tennis/json/dmp_1_" + date + ".json";
  dmps multi_dmp = dmps(file);
  double tau = 2.0;
  multi_dmp.set_time_constant(tau);
  mat M_fast = multi_dmp.evolve(T / tau);

  multi_dmp.set_time_constant(tau / 2);
  mat M_slow = multi_dmp.evolve(T);
  unsigned int N = T / DT;
  uvec idx_sub = linspace<uvec>(1, N - 1, N / 2);
  mat M_slow_subsamp = M_slow.cols(idx_sub);

  EXPECT_TRUE(approx_equal(M_fast, M_slow_subsamp, "absdiff", 0.01));
}

/*
 * Initialize robot posture
 */
static void init_posture(vec7 &q0, int posture, bool verbose) {

  rowvec qinit;
  switch (posture) {
  case 2: // right
    if (verbose)
      cout << "Initializing robot on the right side.\n";
    qinit << 1.0 << -0.2 << -0.1 << 1.8 << -1.57 << 0.1 << 0.3 << endr;
    break;
  case 1: // center
    if (verbose)
      cout << "Initializing robot on the center.\n";
    qinit << 0.0 << 0.0 << 0.0 << 1.5 << -1.75 << 0.0 << 0.0 << endr;
    break;
  case 0: // left
    if (verbose)
      cout << "Initializing robot on the left side\n";
    qinit << -1.0 << 0.0 << 0.0 << 1.5 << -1.57 << 0.1 << 0.3 << endr;
    break;
  default: // default is the right side
    qinit << 1.0 << -0.2 << -0.1 << 1.8 << -1.57 << 0.1 << 0.3 << endr;
    break;
  }
  q0 = qinit.t();
}
