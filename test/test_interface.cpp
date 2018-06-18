/*
 * Test cases for the SL/ZMQ interface.
 */
#include <armadillo>
#include <thread>
#include <chrono>
#include <boost/test/unit_test.hpp>
#include "sl_interface.h"
#include "zmqpp/zmqpp.hpp"
#include "json.hpp"
#include "tabletennis.h"
#include "kalman.h"
#include "optim.h"
#include "kinematics.hpp"
#include <stdio.h>
#include <unistd.h>

using namespace player;
using namespace arma;
using json = nlohmann::json;

static void pub_zmq(const int & num_sent) {
    const std::string sendpoint = "tcp://*:7660";

    // initialize the 0MQ context
    zmqpp::context context;

    // generate a push socket
    zmqpp::socket_type type = zmqpp::socket_type::pub;
    zmqpp::socket socket (context, type);

    // open the connection
    BOOST_TEST_MESSAGE("Binding to " << sendpoint << "...");
    socket.bind(sendpoint);

    // pause to connect
    std::this_thread::sleep_for(std::chrono::milliseconds(1000));

    // send a message
    BOOST_TEST_MESSAGE("Sending increasing numbers...");

    for (int i = 0; i < num_sent; i++) {
        zmqpp::message message;
        // compose a message from a string and a number
        message << i;
        socket.send(message);
        std::this_thread::sleep_for(std::chrono::microseconds(1000));
    }

    BOOST_TEST_MESSAGE("Sent all messages.");
    //socket.disconnect(sendpoint);
}

static void sub_zmq(const int & num_max, int & num_received) {
    const std::string url = "tcp://localhost:7660";

    // initialize the 0MQ context
    zmqpp::context context;

    // generate a subscribe socket
    zmqpp::socket_type type = zmqpp::socket_type::sub;
    zmqpp::socket socket (context, type);
    socket.subscribe(""); // subscribe to the default channel

    // connect to the socket
    BOOST_TEST_MESSAGE("Connecting to " << url << "...");
    socket.connect(url);

    // receive the message
    BOOST_TEST_MESSAGE("Receiving messages...");

    while (num_received < num_max) {
        zmqpp::message message;
        // decompose the message
        socket.receive(message);
        int number;
        message >> number;
        //BOOST_TEST_MESSAGE("Received " << number);
        num_received++;
    }
    BOOST_TEST_MESSAGE("Received all messages...");
    //socket.disconnect(url);
}

static void pub_ball(const unsigned & num_balls) {
    // create zmq server
    zmqpp::context ctx;
    auto socket_type = zmqpp::socket_type::pub;
    zmqpp::socket position_push = zmqpp::socket(ctx, socket_type);
    const std::string & url = "tcp://*:4242";
    position_push.bind(url);

    // pause to connect
    std::this_thread::sleep_for(std::chrono::milliseconds(1000));

    // create table tennis ball and create path
    EKF filter = init_filter();
    TableTennis tt = TableTennis();
    tt.set_ball_gun(0.01);
    filter.set_prior(tt.get_ball_state(),0.01*eye<mat>(6,6));
    mat ball_path = filter.predict_path(0.002,num_balls);

    //send message to listener object
    for (unsigned int i = 0; i < num_balls; i++) {
        zmqpp::message mes;
        vec3 obs3 = ball_path.col(i).rows(0,2);
        json obs = {obs3(0),obs3(1),obs3(2)};
        json jobs {
            {"num", i},
            {"time", 0.002*i},
            {"obs", obs}
        };
        mes << jobs.dump();
        position_push.send(mes);
    }
}

void test_zmqpp() {

    using std::thread;
    BOOST_TEST_MESSAGE("\nTesting ZMQPP in sub/pub mode...");

    const int num_sent = 1000;
    int num_received = 0;
    thread t_sub(sub_zmq,std::ref(num_sent),std::ref(num_received));
    t_sub.detach();
    thread t_pub(pub_zmq,std::ref(num_sent));
    t_pub.join();
    BOOST_TEST_MESSAGE("Finished.");
    BOOST_TEST(num_sent == num_received);
}

void test_zmq_listener() {

    BOOST_TEST_MESSAGE("\nTesting ZMQ connection for ball position info...");

    //start listener server
    Listener listener("tcp://localhost:4242");
    int num_balls = 500;
    pub_ball(num_balls);
    // stop listener
    sleep(0.1);
    listener.stop();

    //print data
    int num_received = listener.give_info();
    //check length of data
    BOOST_TEST_MESSAGE("Checking equality of num. sent balls and received balls...");
    BOOST_TEST(num_received == num_balls);
}

void test_new_interface() {

    BOOST_TEST_MESSAGE("\nTesting new SL interface...");

    using namespace player;
    using namespace optim;
    using namespace const_tt;

    // create zmq server
    BOOST_TEST_MESSAGE("Creating ZMQ server...");
    zmqpp::context ctx;
    auto socket_type = zmqpp::socket_type::pub;
    zmqpp::socket position_pub = zmqpp::socket(ctx, socket_type);
    const std::string & url = "tcp://*:7660";
    position_pub.bind(url);

    // pause to connect
    std::this_thread::sleep_for(std::chrono::milliseconds(1000));

    // create table tennis ball
    BOOST_TEST_MESSAGE("Setting up table tennis setup...");
    TableTennis tt = TableTennis(false,true);
    tt.set_ball_gun(0.001);
    int N = 2000;
    SL_Jstate joint_state[NDOF+1];
    SL_DJstate joint_des_state[NDOF+1];
    rowvec qinit;
    qinit << 0.0 << 0.0 << 0.0 << 0.0 << 1.5 << -1.75 << 0.0 << 0.0 << endr; // center pose
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

    BOOST_TEST_MESSAGE("Loading settings from config file...");
    load_options();

    //robot.reset_filter(std_model,std_noise);
    for (int i = 0; i < N; i++) { // one trial
        vec3 obs3 = tt.get_ball_position();
        zmqpp::message mes;
        json obs = {obs3(0),obs3(1),obs3(2)};
        json jobs {
            {"num", i},
            {"time", DT*i},
            {"obs", obs}
        };
        mes << jobs.dump();
        position_pub.send(mes);
        std::this_thread::sleep_for(std::chrono::milliseconds(2));
        play_new(joint_state, joint_des_state);
        for (int i = 0; i < NDOF; i++) {
            qdes.q(i) = joint_des_state[i+1].th;
            qdes.qd(i) = joint_des_state[i+1].thd;
            qdes.qdd(i) = joint_des_state[i+1].thdd;
        }
        calc_racket_state(qdes,robot_racket);
        tt.integrate_ball_state(robot_racket,DT);
        for (int i = 0; i < NDOF; i++) {
            joint_state[i+1].th = qdes.q(i);
            joint_state[i+1].thd = qdes.qd(i);
            joint_state[i+1].thdd = qdes.qdd(i);
        }
    }
    BOOST_TEST(tt.has_legally_landed());
}
