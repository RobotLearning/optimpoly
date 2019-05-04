#include "ball_interface.h"
#include "json.hpp"
#include "kalman.h"
#include "kinematics.hpp"
#include "optim.h"
#include "sl_interface.h"
#include "tabletennis.h"
#include "zmqpp/zmqpp.hpp"
#include "gtest/gtest.h"
#include <armadillo>
#include <chrono>
#include <stdio.h>
#include <thread>
#include <unistd.h>

using namespace player;
using namespace optim;
using namespace const_tt;
using namespace arma;
using json = nlohmann::json;

static mat create_ball_path(const unsigned &num_balls);

static void pub_zmq(const int &num_sent) {
  const std::string sendpoint = "tcp://*:7660";

  // initialize the 0MQ context
  zmqpp::context context;

  // generate a push socket
  zmqpp::socket_type type = zmqpp::socket_type::pub;
  zmqpp::socket socket(context, type);

  // open the connection
  std::cout << "Binding to " << sendpoint << "...\n";
  socket.bind(sendpoint);

  // pause to connect
  std::this_thread::sleep_for(std::chrono::milliseconds(1000));

  // send a message
  std::cout << "Sending increasing numbers...\n";

  for (int i = 0; i < num_sent; i++) {
    zmqpp::message message;
    // compose a message from a string and a number
    message << i;
    socket.send(message);
    std::this_thread::sleep_for(std::chrono::microseconds(1000));
  }

  std::cout << "Sent all messages.";
  // socket.disconnect(sendpoint);
}

static void sub_zmq(const int &num_max, int &num_received) {
  const std::string url = "tcp://localhost:7660";

  // initialize the 0MQ context
  zmqpp::context context;

  // generate a subscribe socket
  zmqpp::socket_type type = zmqpp::socket_type::sub;
  zmqpp::socket socket(context, type);
  socket.subscribe(""); // subscribe to the default channel

  // connect to the socket
  std::cout << "Connecting to " << url << "...\n";
  socket.connect(url);

  // receive the message
  std::cout << "Receiving messages...\n";

  while (num_received < num_max) {
    zmqpp::message message;
    // decompose the message
    socket.receive(message);
    int number;
    message >> number;
    num_received++;
  }
  std::cout << "Received all messages...\n";
  // socket.disconnect(url);
}

static void pub_ball(const unsigned &num_balls) {
  // create zmq server
  zmqpp::context ctx;
  auto socket_type = zmqpp::socket_type::pub;
  zmqpp::socket position_push = zmqpp::socket(ctx, socket_type);
  const std::string &url = "tcp://*:4242";
  position_push.bind(url);

  // pause to connect
  std::this_thread::sleep_for(std::chrono::milliseconds(1000));

  // create table tennis ball and create path
  mat ball_path = create_ball_path(num_balls);

  // send message to listener object
  for (unsigned int i = 0; i < num_balls; i++) {
    zmqpp::message mes;
    vec3 obs3 = ball_path.col(i).rows(0, 2);
    json obs = {obs3(0), obs3(1), obs3(2)};
    json jobs{{"num", i}, {"time", 0.002 * i}, {"obs", obs}};
    mes << jobs.dump();
    position_push.send(mes);
  }
}

static mat create_ball_path(const unsigned &num_balls) {
  EKF filter = init_ball_filter();
  TableTennis tt = TableTennis();
  tt.set_ball_gun(0.01);
  filter.set_prior(tt.get_ball_state(), 0.01 * eye<mat>(6, 6));
  return filter.predict_path(0.002, num_balls);
}

TEST(TestZMQ, TestZMQPPInSubPubMode) {

  using std::thread;
  std::cout << "\nTesting ZMQPP in sub/pub mode...\n";

  const int num_sent = 1000;
  int num_received = 0;
  thread t_sub(sub_zmq, std::ref(num_sent), std::ref(num_received));
  t_sub.detach();
  thread t_pub(pub_zmq, std::ref(num_sent));
  t_pub.join();
  std::cout << "Finished.\n";
  EXPECT_EQ(num_sent, num_received);
}

TEST(TestZMQ, TestZMQPForEqualityOfSentAndReceivedBalls) {

  // start listener server
  Listener listener("tcp://localhost:4242", false, false, "DLT");
  int num_sent_balls = 500;
  pub_ball(num_sent_balls);
  // stop listener
  std::this_thread::sleep_for(std::chrono::milliseconds(1000));

  // print data
  int num_received_balls = listener.give_info();
  listener.stop();
  EXPECT_EQ(num_received_balls, num_sent_balls);
}

TEST(TestBallListener, Test3DInterfaceForBallListener) {

  // create zmq server
  zmqpp::context ctx;
  auto socket_type = zmqpp::socket_type::pub;
  zmqpp::socket position_pub = zmqpp::socket(ctx, socket_type);
  const std::string &url = "tcp://*:4243";
  position_pub.bind(url);

  // pause to connect
  std::this_thread::sleep_for(std::chrono::milliseconds(1000));

  // create table tennis ball
  TableTennis tt = TableTennis(false, true);
  tt.set_ball_gun(0.001);
  int N = 2000;
  SL_Jstate joint_state[NDOF + 1];
  SL_DJstate joint_des_state[NDOF + 1];
  rowvec qinit;
  qinit << 0.0 << 0.0 << 0.0 << 0.0 << 1.5 << -1.75 << 0.0 << 0.0
        << endr; // center pose
  for (int i = 0; i <= NDOF; i++) {
    joint_state[i].th = qinit(i);
    joint_state[i].thd = 0.0;
    joint_state[i].thdd = 0.0;
    joint_state[i].u = 0.0;
    joint_state[i].ufb = 0.0;
    joint_state[i].load = 0.0;
    joint_des_state[i].th = 0.0;
    joint_des_state[i].thd = 0.0;
    joint_des_state[i].thdd = 0.0;
    joint_des_state[i].uex = 0.0;
    joint_des_state[i].uff = 0.0;
  }
  joint qdes, qact;
  qdes.q = qinit.tail(NDOF).t();
  racket robot_racket;

  void *vflags = load_player_options(); // load settings from config file
  player_flags *pflags = reinterpret_cast<player_flags *>(vflags);
  // 3d ball info sent and received
  pflags->listen_2d = false;
  pflags->zmq_url = "tcp://localhost:4243";
  // pflags->debug_vision = true;

  // robot.reset_filter(std_model,std_noise);
  for (int i = 0; i < N; i++) { // one trial
    vec3 obs3 = tt.get_ball_position();
    zmqpp::message mes;
    json obs = {obs3(0), obs3(1), obs3(2)};
    json jobs{{"num", i}, {"time", DT * i}, {"obs", obs}};
    mes << jobs.dump();
    position_pub.send(mes);
    std::this_thread::sleep_for(std::chrono::milliseconds(2));
    play(joint_state, joint_des_state);
    for (int i = 0; i < NDOF; i++) {
      qdes.q(i) = joint_des_state[i + 1].th;
      qdes.qd(i) = joint_des_state[i + 1].thd;
      qdes.qdd(i) = joint_des_state[i + 1].thdd;
    }
    calc_racket_state(qdes, robot_racket);
    tt.integrate_ball_state(robot_racket, DT);
    for (int i = 0; i < NDOF; i++) {
      joint_state[i + 1].th = qdes.q(i);
      joint_state[i + 1].thd = qdes.qd(i);
      joint_state[i + 1].thdd = qdes.qdd(i);
    }
  }
  EXPECT_TRUE(tt.has_legally_landed());
}

TEST(TestBallListener, Test2DInterfaceForBallListener) {

  const int NUM_CAM = 2;
  // create zmq server
  zmqpp::context ctx;
  auto socket_type = zmqpp::socket_type::pub;
  zmqpp::socket position_pub = zmqpp::socket(ctx, socket_type);
  const std::string &url = "tcp://*:4244";
  position_pub.bind(url);

  // pause to connect
  std::this_thread::sleep_for(std::chrono::milliseconds(1000));

  // create table tennis ball
  TableTennis tt = TableTennis(false, true);
  tt.set_ball_gun(0.001);
  int N = 2000;
  SL_Jstate joint_state[NDOF + 1];
  SL_DJstate joint_des_state[NDOF + 1];
  rowvec qinit;
  qinit << 0.0 << 0.0 << 0.0 << 0.0 << 1.5 << -1.75 << 0.0 << 0.0
        << endr; // center pose
  for (int i = 0; i <= NDOF; i++) {
    joint_state[i].th = qinit(i);
    joint_state[i].thd = 0.0;
    joint_state[i].thdd = 0.0;
    joint_state[i].u = 0.0;
    joint_state[i].ufb = 0.0;
    joint_state[i].load = 0.0;
    joint_des_state[i].th = 0.0;
    joint_des_state[i].thd = 0.0;
    joint_des_state[i].thdd = 0.0;
    joint_des_state[i].uex = 0.0;
    joint_des_state[i].uff = 0.0;
  }
  joint qdes, qact;
  qdes.q = qinit.tail(NDOF).t();
  racket robot_racket;

  void *vflags = load_player_options();
  player_flags *pflags = reinterpret_cast<player_flags *>(vflags);
  // 3d ball info sent and received
  pflags->listen_2d = true;
  pflags->zmq_url = "tcp://localhost:4244";
  pflags->debug_vision = false;
  pflags->time_land_des += 0.1;
  pflags->ball_land_des_offset[Y] -= 0.5;

  std::map<unsigned, mat34> proj_mats; // Calibration matrices
  proj_mats = load_proj_mats("server_3d_conf_ping_okan.json");

  // robot.reset_filter(std_model,std_noise);
  for (int i = 0; i < N; i++) { // one trial
    vec3 obs3 = tt.get_ball_position();
    zmqpp::message mes;
    for (int j = 0; j < NUM_CAM; j++) {
      vec4 obs4;
      obs4.head(3) = obs3;
      obs4(3) = 1.0;
      vec3 proj_obs = proj_mats[j + 2] * obs4; // in hom. coordinates
      json obs = {proj_obs(0), proj_obs(1)};
      json jobs{{"num", i}, {"time", DT * i}, {"obs", obs}, {"cam_id", j + 2}};
      mes << jobs.dump();
      position_pub.send(mes);
    }
    // cout << "Num: " << i << ", Projected: " << obs3.t();
    std::this_thread::sleep_for(std::chrono::milliseconds(2));
    play(joint_state, joint_des_state);
    for (int i = 0; i < NDOF; i++) {
      qdes.q(i) = joint_des_state[i + 1].th;
      qdes.qd(i) = joint_des_state[i + 1].thd;
      qdes.qdd(i) = joint_des_state[i + 1].thdd;
    }
    calc_racket_state(qdes, robot_racket);
    tt.integrate_ball_state(robot_racket, DT);
    for (int i = 0; i < NDOF; i++) {
      joint_state[i + 1].th = qdes.q(i);
      joint_state[i + 1].thd = qdes.qd(i);
      joint_state[i + 1].thdd = qdes.qdd(i);
    }
  }
  EXPECT_TRUE(tt.has_legally_landed());
}

TEST(TestBallListener, TestTriangulationFor3DEstimateFrom2DBallObservations) {

  // create a ball motion
  // pick some balls and project
  // triangulate to 3d and check for error
  const int num_balls = 10;
  const int num_balls_in_path = 100;
  mat ball_path = create_ball_path(num_balls_in_path);
  uvec idx_rand = randi<uvec>(num_balls, distr_param(0, num_balls_in_path));
  mat ball_states = ball_path.cols(idx_rand);
  mat balls = ball_states.rows(X, Z);
  mat hom_balls = join_vert(balls, ones<rowvec>(num_balls));

  std::map<unsigned, mat34> proj_mats =
      load_proj_mats("server_3d_conf_ping_okan.json");
  mat proj_balls_2 = proj_mats[2] * hom_balls;
  mat proj_balls_3 = proj_mats[3] * hom_balls;

  mat balls_pred = zeros<mat>(3, num_balls);
  std::vector<pixels> v_proj_balls;
  for (int i = 0; i < num_balls; i++) {
    pixels p;
    p.cam_id = 2;
    p.vals[0] = proj_balls_2(0, i);
    p.vals[1] = proj_balls_2(1, i);
    v_proj_balls.push_back(p);
    pixels p2;
    p2.cam_id = 3;
    p2.vals[0] = proj_balls_3(0, i);
    p2.vals[1] = proj_balls_3(1, i);
    v_proj_balls.push_back(p2);
    vec3 obs3d;
    bool success = triangulate(proj_mats, v_proj_balls, "invert", obs3d);
    // TODO: DLT method won't work for triangulation here!
    if (success)
      balls_pred.col(i) = obs3d;
  }
  EXPECT_TRUE(approx_equal(balls, balls_pred, "absdiff", 0.001));
  // to pass requires the second triangulation method, i.e. inv. of P1 and P2
  // 1-2 columns stacked
}
