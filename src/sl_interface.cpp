/**
 *
 * @file sl_interface.cpp
 *
 * @brief Interface of the Player class to the SL real-time simulator and to
 * the robot.
 *
 */

#include "sl_interface.h"

#include <armadillo>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <string>
#include <sys/time.h>

#include "ball_interface.h"
#include "dmp.h"
#include "json.hpp"
#include "kalman.h"
#include "player.hpp"
#include "rbf.h"
#include "serve.h"
#include "tabletennis.h"
#include "yaml-cpp/yaml.h"

using namespace arma;
using namespace player;
using namespace serve;
using namespace const_tt;

player_flags pflags; //!< global structure for setting Player task options
serve_flags sflags;  //!< global structure for setting Serve task options

/*
 *  Set algorithm to initialize Player with.
 *  alg_num selects between three algorithms: VHP/FOCUSED/DEFENSIVE.
 */
static void set_player_algorithm(const int alg_num) {

  switch (alg_num) {
  case 0:
    std::cout << "Setting to FOCUSED player..." << std::endl;
    pflags.alg = FOCUS;
    break;
  case 1:
    std::cout << "Setting to DEFENSIVE player..." << std::endl;
    pflags.alg = DP;
    break;
  case 2:
    std::cout << "Setting to VHP player..." << std::endl;
    pflags.alg = VHP;
    break;
  default:
    pflags.alg = FOCUS;
  }
}

void *load_player_options() {

  std::string home = std::getenv("HOME");
  std::string config_file = home + "/projects/table-tennis/config/" + "player.yaml";

  YAML::Node config = YAML::LoadFile(config_file);
  YAML::Node player_conf = config["player"];
  set_player_algorithm(player_conf["algorithm"].as<int>());
  pflags.mpc = player_conf["mpc"]["enable"].as<bool>();
  pflags.freq_mpc = player_conf["mpc"]["frequency"].as<int>();
  pflags.verbosity = player_conf["verbose"].as<int>();
  pflags.save = player_conf["save_data"].as<bool>();
  pflags.check_bounce = player_conf["check_bounce"].as<bool>();
  pflags.optim_rest_posture = player_conf["optim_rest_posture"].as<bool>();
  YAML::Node filter_conf = config["filter"];
  pflags.spin_based_pred = filter_conf["model"]["spin_based_pred"].as<bool>();
  pflags.var_noise = filter_conf["model"]["var_obs_noise"].as<double>();
  pflags.var_model = filter_conf["model"]["var_process_noise"].as<double>();
  pflags.min_obs = filter_conf["reset"]["min_obs"].as<int>();
  pflags.t_reset_thresh = filter_conf["reset"]["time_to_reset"].as<double>();
  pflags.outlier_detection =
      filter_conf["outlier_detection"]["enable"].as<bool>();
  pflags.out_reject_mult =
      filter_conf["outlier_detection"]["multiplier"].as<double>();
  YAML::Node optim_conf = config["optim"];
  pflags.optim_offset = optim_conf["start_optim_offset"].as<double>();
  pflags.time2return = optim_conf["time_to_return"].as<double>();
  pflags.ball_land_des_offset =
      optim_conf["ball_land_des_offset"].as<std::vector<double>>();
  pflags.time_land_des = optim_conf["time_land_des"].as<double>();
  // defensive player options are extra
  pflags.weights =
      optim_conf["defensive_player"]["weights"].as<std::vector<double>>();
  pflags.penalty_loc =
      optim_conf["defensive_player"]["penalty_loc"].as<std::vector<double>>();
  pflags.mult_vel = optim_conf["defensive_player"]["mult_outgoing_ball_vel"]
                        .as<std::vector<double>>();
  pflags.VHPY = optim_conf["virtual_hitting_player"]["y_location"].as<double>();
  YAML::Node network_conf = config["network"];
  pflags.zmq_url = network_conf["url"].as<std::string>();
  pflags.debug_vision = network_conf["vision"]["debug"].as<bool>();
  pflags.listen_2d = network_conf["vision"]["listen_2d"].as<bool>();
  pflags.triangulation = network_conf["vision"]["triangulation"].as<std::string>();

  pflags.detach = true; // always detached in SL/REAL ROBOT!
  pflags.reset = true;
  return reinterpret_cast<void *>(&pflags);
}

void play(const SL_Jstate joint_state[NDOF + 1],
          SL_DJstate joint_des_state[NDOF + 1]) {

  // connect to ZMQ server
  // acquire ball info from ZMQ server
  // if new ball add status true else false
  // call play function
  static Listener *listener;

  // since old code support multiple blobs
  static ball_obs blob;
  static vec7 q0;
  static optim::joint qact;
  static optim::joint qdes;
  static Player *robot = nullptr;

  if (pflags.reset) {
    delete listener;
    listener =
        new Listener(pflags.zmq_url, pflags.listen_2d, pflags.debug_vision, pflags.triangulation);
    for (int i = 0; i < NDOF; i++) {
      qdes.q(i) = q0(i) = joint_state[i + 1].th;
      qdes.qd(i) = 0.0;
      qdes.qdd(i) = 0.0;
    }
    delete robot;
    robot = new Player(q0, pflags);
    pflags.reset = false;
  } else {
    // update ball info
    listener->fetch(blob);

    for (int i = 0; i < NDOF; i++) {
      qact.q(i) = joint_state[i + 1].th;
      qact.qd(i) = joint_state[i + 1].thd;
      qact.qdd(i) = joint_state[i + 1].thdd;
    }
    robot->play(blob, qact, qdes);
  }

  // update desired joint state
  for (int i = 0; i < NDOF; i++) {
    joint_des_state[i + 1].th = qdes.q(i);
    joint_des_state[i + 1].thd = qdes.qd(i);
    joint_des_state[i + 1].thdd = qdes.qdd(i);
  }
}

void cheat(const SL_Jstate joint_state[NDOF + 1],
           const SL_Cstate sim_ball_state,
           SL_DJstate joint_des_state[NDOF + 1]) {

  static vec7 q0;
  static vec6 ball_state;
  static optim::joint qact;
  static optim::joint qdes;
  static Player *cp = nullptr; // centered player

  if (pflags.reset) {
    for (int i = 0; i < NDOF; i++) {
      qdes.q(i) = q0(i) = joint_state[i + 1].th;
      qdes.qd(i) = 0.0;
      qdes.qdd(i) = 0.0;
    }
    cp = new Player(q0, pflags);
    pflags.reset = false;
  } else {
    for (int i = 0; i < NDOF; i++) {
      qact.q(i) = joint_state[i + 1].th;
      qact.qd(i) = joint_state[i + 1].thd;
      qact.qdd(i) = joint_state[i + 1].thdd;
    }
    for (int i = 0; i < NCART; i++) {
      ball_state(i) = sim_ball_state.x[i + 1];
      ball_state(i + NCART) = sim_ball_state.xd[i + 1];
    }
    cp->cheat(qact, ball_state, qdes);
  }

  // update desired joint state
  for (int i = 0; i < NDOF; i++) {
    joint_des_state[i + 1].th = qdes.q(i);
    joint_des_state[i + 1].thd = qdes.qd(i);
    joint_des_state[i + 1].thdd = qdes.qdd(i);
  }
}

void save_joint_data(const SL_Jstate joint_state[NDOF + 1],
                     const SL_DJstate joint_des_state[NDOF + 1],
                     const int save_qdes, const int reset) {

  using namespace std::chrono;
  const std::string home = std::getenv("HOME");
  const std::string joint_file = home + "/table-tennis/joints.txt";

  static rowvec q = zeros<rowvec>(NDOF);
  static rowvec qdes = zeros<rowvec>(NDOF);
  static std::ofstream stream_joints;

  if (reset) {
    stream_joints.close();
    stream_joints.open(joint_file, std::ofstream::out);
  } else {
    milliseconds ms =
        duration_cast<milliseconds>(system_clock::now().time_since_epoch());
    double time = ms.count();
    for (int i = 1; i <= NDOF; i++) {
      q(i - 1) = joint_state[i].th;
      if (save_qdes) {
        qdes(i - 1) = joint_des_state[i].th;
      }
    }

    if (stream_joints.is_open()) {
      if (save_qdes) {
        stream_joints << std::fixed << time << join_horiz(q, qdes);
      } else {
        stream_joints << std::fixed << time << q;
      }
    }
    // stream_joints.close();
  }
}

void save_ball_data(const char *url_string, const int listen_2d,
                    const int debug_vision, const int reset) {

  using namespace std::chrono;
  static Listener *listener;
  static ball_obs obs;
  static std::ofstream stream_balls;
  static std::string home = std::getenv("HOME");
  static std::string ball_file = home + "/projects/table-tennis/balls.txt";
  static std::string triangulation_method = "DLT";

  if (reset) {
    delete listener;
    listener = new Listener(url_string, listen_2d, (bool)debug_vision, triangulation_method);
    stream_balls.close();
    stream_balls.open(ball_file, std::ofstream::out);
  } else {
    milliseconds ms =
        duration_cast<milliseconds>(system_clock::now().time_since_epoch());
    double time = ms.count();
    // update blobs structure with new ball data
    listener->fetch(obs);
    if (obs.status) {
      if (stream_balls.is_open()) {
        stream_balls << std::fixed << time << obs.pos.t();
      }
    }
  }
}

void *load_serve_options(double custom_pose[], serve_task_options *options) {

  const std::string home = std::getenv("HOME");
  std::string config_file = home + "/projects/table-tennis/config/serve.yaml";

  YAML::Node config = YAML::LoadFile(config_file);
  double perturb_init = config["perturb_init"].as<double>();
  options->use_inv_dyn_fb = config["use_inv_dyn_fb"].as<bool>();

  YAML::Node serve_conf = config["serve"];
  sflags.alg = serve_conf["algorithm"].as<int>();
  sflags.speedup = serve_conf["speed_up"].as<double>();
  sflags.start_from_act_state = serve_conf["start_from_act_state"].as<bool>();
  sflags.json_file = serve_conf["json_file"].as<std::string>();
  // save settings
  sflags.save_joint_act_data = serve_conf["save"]["act_joint"].as<bool>();
  sflags.save_joint_des_data = serve_conf["save"]["des_joint"].as<bool>();
  sflags.save_ball_data = serve_conf["save"]["ball_data"].as<bool>();

  YAML::Node network_conf = config["network"];
  sflags.zmq_url = network_conf["url"].as<std::string>();
  sflags.debug_vision = network_conf["vision"]["debug"].as<bool>();
  sflags.listen_2d = network_conf["vision"]["listen_2d"].as<bool>();
  sflags.triangulation = network_conf["vision"]["triangulation"].as<std::string>();

  YAML::Node optim_conf = config["optim"];
  sflags.detach = config["optim"]["detach"].as<bool>();
  sflags.mpc = config["mpc"]["enable"].as<bool>();
  sflags.freq_mpc = config["mpc"]["frequency"].as<int>();
  sflags.time_land_des = config["optim"]["time_land_des"].as<double>();
  sflags.ball_land_des_offset = config["optim"]["ball_land_des_offset"].as<std::vector<double>>();

  std::string json_file = home + "/projects/table-tennis/json/" + sflags.json_file;
  vec7 pose;
  if (sflags.alg == 1) {
    RBF rbfs = RBF(json_file);
    rbfs.get_init_pos(pose);
  } else if (sflags.alg == 2) {
    CRBF crbfs = CRBF(json_file);
    crbfs.get_init_pos(pose);
  } else {
    using dmps = Joint_DMPs;
    dmps multi_dmp = dmps(json_file);
    multi_dmp.get_init_pos(pose);
  }
  pose += perturb_init * randn(7, 1);

  for (int i = 0; i < NDOF; i++) {
    custom_pose[i] = pose(i);
  }
  sflags.reset = true;
  return reinterpret_cast<void *>(&sflags);
}

void serve_ball(const SL_Jstate joint_state[], SL_DJstate joint_des_state[]) {

  static Listener *listener;
  static ball_obs blob;
  static optim::joint qact;
  static optim::joint qdes;
  static ServeBall<dmps> *robot_dmp = nullptr; // class for serving ball
  static ServeBall<RBF> *robot_rbf = nullptr;
  static ServeBall<CRBF> *robot_crbf = nullptr;

  if (sflags.reset) {
    delete listener;
    listener =
        new Listener(sflags.zmq_url, sflags.listen_2d, sflags.debug_vision, sflags.triangulation);
    vec7 qinit;
    if (sflags.start_from_act_state) {
      for (int i = 0; i < NDOF; i++) {
        qinit(i) = joint_state[i + 1].th;
      }
    }
    for (int i = 0; i < NDOF; i++) {
      qdes.q(i) = joint_state[i + 1].th;
      qdes.qd(i) = 0.0;
      qdes.qdd(i) = 0.0;
    }
    if (sflags.alg == 1) {
      robot_rbf = new ServeBall<RBF>(sflags, qinit);
    } else if (sflags.alg == 2) {
      robot_crbf = new ServeBall<CRBF>(sflags, qinit);
    } else {
      robot_dmp = new ServeBall<dmps>(sflags, qinit);
    }
    sflags.reset = false;
  } else {
    // update blobs structure with new ball data
    listener->fetch(blob);
    for (int i = 0; i < NDOF; i++) {
      qact.q(i) = joint_state[i + 1].th;
      qact.qd(i) = joint_state[i + 1].thd;
      qact.qdd(i) = joint_state[i + 1].thdd;
    }
    if (sflags.alg == 1) {
      robot_rbf->serve(blob, qact, qdes);
    } else if (sflags.alg == 2) {
      robot_crbf->serve(blob, qact, qdes);
    } else {
      robot_dmp->serve(blob, qact, qdes);
    }
  }

  // update desired joint state
  for (int i = 0; i < NDOF; i++) {
    joint_des_state[i + 1].th = qdes.q(i);
    joint_des_state[i + 1].thd = qdes.qd(i);
    joint_des_state[i + 1].thdd = qdes.qdd(i);
  }
}
