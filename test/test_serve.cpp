#include <boost/test/unit_test.hpp>
#include <iostream>
#include <armadillo>
#include "constants.h"
#include "dmp.h"

using namespace arma;
using namespace const_tt;
using namespace DMP;
using dmps = Joint_DMPs;

static void init_posture(vec7 & q0, int posture, bool verbose);

void test_evolve_dmp() {

    // create a dmp
    const double T = 1.0;
    std::string file = "json/dmp2.json";
    dmps multi_dmp = dmps(file);
    mat M = multi_dmp.evolve(T);
    BOOST_TEST_MESSAGE("\nEvolving DMP: " << M.tail_cols(1).t());

    // reinitialize in a different init. state
    vec7 q0;
    init_posture(q0,2,0);
    multi_dmp.set_init_pos(q0);
    // evolve it
    mat M2 = multi_dmp.evolve(T);
    BOOST_TEST_MESSAGE("Evolving DMP from a different posture: " << M2.tail_cols(1).t());

    // change goal state slightly
    vec goal = zeros<vec>(NDOF);
    multi_dmp.get_goal_pos(goal);
    BOOST_TEST_MESSAGE("Goal positions: " << goal.t());
    vec goal_new = goal + 0.01*randn(NDOF);
    multi_dmp.set_goal_pos(goal_new);
    // evolve it
    mat M3 = multi_dmp.evolve(T);
    BOOST_TEST_MESSAGE("Evolving DMP to a slightly different goal: " << M3.tail_cols(1).t());

    BOOST_TEST(approx_equal(goal,M.tail_cols(1),"absdiff",0.01));
    BOOST_TEST(approx_equal(goal,M2.tail_cols(1),"absdiff",0.01));
    BOOST_TEST(approx_equal(goal_new,M3.tail_cols(1),"absdiff",0.01));

}

void test_speedup_dmp() {

    // evolve dmp fast
    // evolve dmp slowly
    // subsample to check if equal
    BOOST_TEST_MESSAGE("\nComparing sped-up dmp evolution to subsampled normal dmp...");

    const double T = 1.0;
    std::string file = "json/dmp4.json";
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
