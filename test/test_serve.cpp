#include <boost/test/unit_test.hpp>
#include <iostream>
#include <armadillo>
#include <thread>
#include "kinematics.hpp"
#include "tabletennis.h"
#include "optim.h"
#include "constants.h"
#include "dmp.h"
#include "serve.h"

using namespace arma;
using namespace const_tt;
using namespace player;
using namespace serve;
using dmps = Joint_DMPs;
using vec_str = std::vector<std::string>;

static void init_posture(vec7 & q0, int posture, bool verbose);

vec6 init_ball_vertical(const double & T,
						const double & std_noise,
                       const dmps & multi_dmp);

void test_serve() {

    // start evolving dmp
    // have a ball coming down to goal state
    // track the ball with a filter
    // if we predict miss then run an optimizer to correct
    // switch to optimizer

    BOOST_TEST_MESSAGE("\nTesting OPTIMIZATION after DMP...");
    const double T = 1.0;
    dmps multi_dmp = init_dmps();
    optim::joint qdes, qact;
    TableTennis tt = TableTennis(false,true);
    double std_init_noise = 0.0;
    vec6 init_ball_state = init_ball_vertical(T,std_init_noise,multi_dmp);
    tt.set_ball_state(init_ball_state);
    racket racket_state;
    ServeBall server = ServeBall(multi_dmp);
    //sflags flags;
    //server.set_flags(flags);

    // check interaction and make sure ball is served correctly
    while (!tt.touched_ground() && !tt.was_legally_served()) {

        tt.integrate_ball_state(racket_state,DT);
        ball_obs obs;
        obs.status = true;
        obs.pos = tt.get_ball_position();
        server.serve(obs,qact,qdes);
        qact = qdes;
        calc_racket_state(qact,racket_state);
        //double dist = norm(obs.pos - racket_state.pos);
        //BOOST_TEST_MESSAGE("Dist between racket and ball: " << dist);

    }
    BOOST_TEST(tt.was_legally_served());

}

/*
 * Initialize ball vertically such that assuming only gravity acts on the
 * ball (i.e. no airdrag) the ball will be around the DMP goal cartesian position
 * T seconds from the start. [around = some noise added]
 */
vec6 init_ball_vertical(const double & T,
						const double & std_noise,
                        const dmps & multi_dmp) {
    ball_params params;
    vec7 goal_state;
    optim::joint goal_joint_state;
    multi_dmp.get_goal_pos(goal_state);
    goal_joint_state.q = goal_state;
    racket goal_racket_state;
    calc_racket_state(goal_joint_state,goal_racket_state);
    vec3 init_ball_pos = goal_racket_state.pos + std_noise * randn(3,1);
    init_ball_pos(Z) -= 0.83 * params.gravity * T*T/2.0; // 0.83 - dmp works!
    vec6 init_ball_state = zeros<vec>(6);
    init_ball_state.head(3) = init_ball_pos;

    return init_ball_state;
}

void test_evolve_dmp() {

    // create a dmp
    const double T = 1.0;
    const std::string home = std::getenv("HOME");
    const std::string file = home + "/table-tennis/json/dmp2.json";
    dmps multi_dmp = dmps(file);
    mat M = multi_dmp.evolve(T);
    BOOST_TEST_MESSAGE("\nEvolving DMP to goal: " << M.tail_cols(1).t());

    // reinitialize in a different init. state
    vec7 q0;
    init_posture(q0,2,1);
    multi_dmp.set_init_pos(q0);
    // evolve it
    mat M2 = multi_dmp.evolve(T);
    BOOST_TEST_MESSAGE("Evolving DMP from a different posture: " << M2.tail_cols(1).t());

    // change goal state slightly
    vec goal = zeros<vec>(NDOF);
    multi_dmp.get_goal_pos(goal);
    BOOST_TEST_MESSAGE("Goal positions: " << goal.t());
    vec goal_new = goal + 0.01*randn(NDOF,1);
    multi_dmp.set_goal_pos(goal_new);
    // evolve it
    mat M3 = multi_dmp.evolve(T);
    BOOST_TEST_MESSAGE("Evolving DMP to a slightly different goal: " << M3.tail_cols(1).t());

    BOOST_TEST(approx_equal(goal,M.tail_cols(1),"absdiff",0.01));
    BOOST_TEST(approx_equal(goal,M2.tail_cols(1),"absdiff",0.01));
    BOOST_TEST(approx_equal(goal_new,M3.tail_cols(1),"absdiff",0.01));

}

void test_dmp_acc() {

    // for all the json files
    // check max acc with goal/init pos adjusted
    const double MAX_ACC_ALLOWED = 50.0;
    const std::string home = std::getenv("HOME");
    vec_str files = get_files(home + "/table-tennis/json/", "dmp");
    dmps multi_dmp;

    for (std::string & file : files) {
        //file.erase(std::remove(file.begin(),file.end(),'\n'),file.end());
        BOOST_TEST_MESSAGE("\nTesting maximum accelerations for " << file);
        std::string full_file = home + "/table-tennis/json/" + file;
        //cout << full_file << endl;
        multi_dmp = dmps(full_file);
        //multi_dmp.turn_on_safe_acc();
        mat Q, Qd, Qdd;
        multi_dmp.evolve(1.0,Q,Qd,Qdd);
        double max_acc = max(max(abs(Qdd)));
        //BOOST_TEST_MESSAGE("Maximum acc: " << max_acc);
        BOOST_TEST(max_acc < MAX_ACC_ALLOWED);

        //BOOST_TEST_MESSAGE("Evolving DMP to a slightly different goal");
        vec goal = zeros<vec>(NDOF);
        multi_dmp.get_goal_pos(goal);
        vec goal_new = goal + 0.01*randn(NDOF,1);
        multi_dmp.set_goal_pos(goal_new);
        multi_dmp.evolve(1.0,Q,Qd,Qdd);
        max_acc = max(max(abs(Qdd)));
        //BOOST_TEST_MESSAGE("Maximum acc: " << max_acc);
        BOOST_TEST(max_acc < MAX_ACC_ALLOWED);

        //BOOST_TEST_MESSAGE("Evolving DMP from a slightly different init posture");
        vec init = zeros<vec>(NDOF);
        multi_dmp.get_init_pos(init);
        vec init_new = init + 0.01*randn(NDOF,1);
        multi_dmp.set_init_pos(init_new);
        multi_dmp.evolve(1.0,Q,Qd,Qdd);
        max_acc = max(max(abs(Qdd)));
        //BOOST_TEST_MESSAGE("Maximum acc: " << max_acc);
        BOOST_TEST(max_acc < MAX_ACC_ALLOWED);

        //BOOST_TEST_MESSAGE("Evolving DMP from a VERY different init posture");
        vec7 q0;
        init_posture(q0,2,0);
        multi_dmp.set_init_pos(q0);
        multi_dmp.evolve(1.0,Q,Qd,Qdd);
        max_acc = max(max(abs(Qdd)));
        //BOOST_TEST_MESSAGE("Maximum acc: " << max_acc);
        BOOST_TEST(max_acc < MAX_ACC_ALLOWED);
    }

}

void test_speedup_dmp() {

    // evolve dmp fast
    // evolve dmp slowly
    // subsample to check if equal
    BOOST_TEST_MESSAGE("\nComparing sped-up dmp evolution to subsampled normal dmp...");

    const double T = 1.0;
    const std::string home = std::getenv("HOME");
    const std::string file = home + "/table-tennis/json/dmp4.json";
    dmps multi_dmp = dmps(file);
    double tau = 2.0;
    multi_dmp.set_time_constant(tau);
    mat M_fast = multi_dmp.evolve(T/tau);

    multi_dmp.set_time_constant(tau/2);
    mat M_slow = multi_dmp.evolve(T);
    unsigned int N = T/DT;
    uvec idx_sub = linspace<uvec>(1,N-1,N/2);
    mat M_slow_subsamp = M_slow.cols(idx_sub);

    BOOST_TEST(approx_equal(M_fast,M_slow_subsamp,"absdiff",0.01));
}

/*
 * Initialize robot posture
 */
static void init_posture(vec7 & q0, int posture, bool verbose) {

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
