#include <boost/test/unit_test.hpp>
#include <iostream>
#include <armadillo>
#include "constants.h"
#include "dmp.h"

using namespace arma;
using namespace const_tt;
using namespace DMP;
using dmps = Joint_DMPs;
using vec_str = std::vector<std::string>;

static void init_posture(vec7 & q0, int posture, bool verbose);
static vec_str get_files(std::string folder_name);

void test_evolve_dmp() {

    // create a dmp
    const double T = 1.0;
    const std::string home = std::getenv("HOME");
    const std::string file = home + "/table-tennis/json/dmp2.json";
    cout << file << endl;
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
    const std::string home = std::getenv("HOME");
    vec_str files = get_files("json/");
    dmps multi_dmp;

    for (const std::string & file : files) {
        BOOST_TEST_MESSAGE("\nTesting maximum accelerations for " << file);
        std::string full_file = home + "/table-tennis/json/" + file;
        multi_dmp = dmps(full_file);
        /*mat Q, Qd, Qdd;
        multi_dmp.evolve(1.0,Q,Qd,Qdd);
        double max_acc = max(max(abs(Qdd)));
        BOOST_TEST_MESSAGE("Maximum acc: " << max_acc);

        BOOST_TEST_MESSAGE("Evolving DMP to a slightly different goal");
        vec goal = zeros<vec>(NDOF);
        multi_dmp.get_goal_pos(goal);
        vec goal_new = goal + 0.01*randn(NDOF,1);
        multi_dmp.set_goal_pos(goal_new);
        multi_dmp.evolve(1.0,Q,Qd,Qdd);
        double max_acc_new = max(max(abs(Qdd)));
        BOOST_TEST_MESSAGE("Maximum acc: " << max_acc_new);

        BOOST_TEST_MESSAGE("Evolving DMP from a slightly different init posture");
        vec init = zeros<vec>(NDOF);
        multi_dmp.get_init_pos(init);
        vec init_new = init + 0.01*randn(NDOF,1);
        multi_dmp.set_init_pos(init_new);
        multi_dmp.evolve(1.0,Q,Qd,Qdd);
        max_acc_new = max(max(abs(Qdd)));
        BOOST_TEST_MESSAGE("Maximum acc: " << max_acc_new);

        BOOST_TEST_MESSAGE("Evolving DMP from a VERY different init posture");
        vec7 q0;
        init_posture(q0,2,0);
        multi_dmp.set_init_pos(q0);
        multi_dmp.evolve(1.0,Q,Qd,Qdd);
        max_acc_new = max(max(abs(Qdd)));
        BOOST_TEST_MESSAGE("Maximum acc: " << max_acc_new);*/
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

static vec_str get_files(std::string folder_name) {

    using std::string;
    vec_str files;
    FILE * stream;
    const int max_buffer = 256;
    char buffer[max_buffer];
    string cmd = "ls " + folder_name;
    //cmd.append(" 2>&1");

    stream = popen(cmd.c_str(), "r");
    if (stream) {
        while (!feof(stream)) {
            if (fgets(buffer, max_buffer, stream) != NULL) {
                files.push_back(buffer);
            }
        }
        pclose(stream);
    }
    return files;
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
