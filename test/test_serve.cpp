#include <boost/test/unit_test.hpp>
#include <iostream>
#include <armadillo>
#include <algorithm>
#include <string>
#include <thread>
#include "kinematics.hpp"
#include "tabletennis.h"
#include "optim.h"
#include "constants.h"
#include "dmp.h"

using namespace arma;
using namespace const_tt;
using namespace player;
using namespace optim;
using namespace DMP;
using dmps = Joint_DMPs;
using vec_str = std::vector<std::string>;

static void init_posture(vec7 & q0, int posture, bool verbose);
static vec_str get_files(std::string folder_name);
void correct_with_optim(const joint & Q,
                        const EKF & filter,
                        FocusedOptim & fp,
                        bool & ran_optim,
                        bool & read_optim);
bool predict_ball_hit(const EKF & filter,
                      const vec7 & q_rest_des,
                      const dmps & multi_dmp,
                      const spline_params & p,
                      const double & t_poly,
                      const int N,
                      const int i,
                      const bool & ran_optim);
void estimate_ball_state(const vec3 & obs, EKF & filter);
void serve(const vec7 & q_rest_des,
           const bool & ran_optim,
           const dmps & multi_dmp,
           const EKF & filter,
           const int N,
           const int i,
           FocusedOptim & fp,
           joint & Q,
           double & t_poly,
           spline_params & p,
           bool & read_optim);

void test_optim_with_dmp() {

    // start evolving dmp
    // have a ball coming down to goal state
    // track the ball with a filter
    // if we predict miss then run an optimizer to correct
    // switch to optimizer

    const double T = 1.0;
    const std::string home = std::getenv("HOME");
    vec_str files = get_files(home + "/table-tennis/json/");
    arma_rng::set_seed_random();
    int val = (randi(1,distr_param(0,files.size()-1)).at(0));
    std::string file = files[val];
    BOOST_TEST_MESSAGE("\nTesting OPTIMIZATION after DMP " << file);
    std::string full_file = home + "/table-tennis/json/" + file;
    dmps multi_dmp = dmps(full_file);
    unsigned int N = T/DT;
    //multi_dmp.reset();
    vec7 goal_state;
    vec7 q_rest_des;
    multi_dmp.get_init_pos(q_rest_des);
    optim::joint goal_joint_state;
    multi_dmp.get_goal_pos(goal_state);
    goal_joint_state.q = goal_state;
    racket goal_racket_state;
    calc_racket_state(goal_joint_state,goal_racket_state);
    optim::joint Q;
    TableTennis tt = TableTennis(false,true);
    ball_params params;
    vec3 init_ball_pos = goal_racket_state.pos;
    init_ball_pos(Z) -= 0.85 * params.gravity * T*T/2.0; // 0.83 - dmp works!
    vec6 init_ball_state = zeros<vec>(6);
    init_ball_state.head(3) = init_ball_pos;
    tt.set_ball_state(init_ball_state);
    racket racket_state;
    EKF filter = init_ball_filter(0.03,0.0001);

    bool ran_optim = false;
    bool read_optim = false;
    int i = 0;
    double Tmax = 2.0;
    double lb[2*NDOF+1], ub[2*NDOF+1];
    set_bounds(lb,ub,0.01,Tmax);
    FocusedOptim fp = FocusedOptim(q_rest_des,lb,ub);
    double t_poly = 0.0;
    spline_params p;

    // check interaction and make sure ball is served correctly
    while (!tt.touched_ground() && !tt.was_served()) {

        serve(q_rest_des,ran_optim,multi_dmp,filter,N,i,fp,Q,t_poly,p,read_optim);
        calc_racket_state(Q,racket_state);
        tt.integrate_ball_state(racket_state,DT);
        vec3 obs = tt.get_ball_position();
        estimate_ball_state(obs,filter);
        //double dist = norm(obs - racket_state.pos);
        //BOOST_TEST_MESSAGE("Dist between racket and ball: " << dist);
        i++;
    }
    BOOST_TEST(tt.was_served());

}

void serve(const vec7 & q_rest_des,
           const bool & ran_optim,
           const dmps & multi_dmp,
           const EKF & filter,
           const int N,
           const int i,
           FocusedOptim & fp,
           joint & Q,
           double & t_poly,
           spline_params & p,
           bool & read_optim) {
    if (ran_optim) {
        if (!read_optim) {
            fp.get_params(Q,p);
            t_poly = 0.0;
            read_optim = true;
        }

        update_next_state(p,q_rest_des,1.0,t_poly,Q);
    }
    else {
        multi_dmp.step(DT,Q);
    }

    if (!predict_ball_hit(filter,q_rest_des,multi_dmp,p,t_poly,N,i,ran_optim)) {
        // launch optimizer
        BOOST_TEST_MESSAGE("Ball won't be hit! Launching optimization to correct traj...");
        correct_with_optim(Q,filter,fp,ran_optim,read_optim);
    }
    else {
        //BOOST_TEST_MESSAGE("Predicting ball hit...");
    }
}

void estimate_ball_state(const vec3 & obs, EKF & filter) {

    const int min_ball_to_start_filter = 5;
    static mat init_obs = zeros<mat>(3,min_ball_to_start_filter);
    static vec init_times = zeros<vec>(min_ball_to_start_filter);
    static int idx = 0;

    // FILTERING
    if (idx == min_ball_to_start_filter) {
        vec6 init_filter_state;
        estimate_ball_linear(init_obs,init_times,false,init_filter_state);
        mat P;
        P.eye(6,6);
        filter.set_prior(init_filter_state,P);
    }
    else if (idx < min_ball_to_start_filter) {
        init_obs.col(idx) = obs;
        init_times(idx) = idx*DT;
    }
    else {
        filter.predict(DT,true);
        filter.update(obs);
    }
}

bool predict_ball_hit(const EKF & filter,
                      const vec7 & q_rest_des,
                      const dmps & multi_dmp,
                      const spline_params & p,
                      const double & t_poly,
                      const int N,
                      const int i,
                      const bool & ran_optim) {

    // PREDICT IF BALL WILL BE HIT, OTHERWISE LAUNCH TRAJ OPTIM

    try {
        TableTennis tt_pred = TableTennis(filter.get_mean(),false,false);
        racket robot_pred_racket;
        optim::joint Q_pred;
        unsigned N_pred;

        if (!ran_optim) {
            dmps multi_pred_dmp = multi_dmp;
            N_pred = N-i;
            for (unsigned j = 0; j < N_pred; j++) {
                multi_pred_dmp.step(DT,Q_pred);
                calc_racket_state(Q_pred,robot_pred_racket);
                tt_pred.integrate_ball_state(robot_pred_racket,DT);
            }
        }
        else {
            double t_poly_pred = t_poly;
            N_pred = (p.time2hit - t_poly)/DT;
            for (unsigned j = 0; j < N_pred; j++) {
                update_next_state(p,q_rest_des,1.0,t_poly_pred,Q_pred);
                calc_racket_state(Q_pred,robot_pred_racket);
                tt_pred.integrate_ball_state(robot_pred_racket,DT);
            }
        }
        return tt_pred.was_served();
    }
    catch (const std::exception & not_init_error) {
        // filter not yet initialized
        return true;
    }
}

void correct_with_optim(const joint & Q,
                        const EKF & filter,
                        FocusedOptim & fp,
                        bool & ran_optim,
                        bool & read_optim) {

    double Tpred = 2.0;
    mat balls_pred = zeros<mat>(6,Tpred/DT);
    predict_ball(Tpred,balls_pred,filter);

    //lookup_soln(filter.get_mean(),1,qact);
    double time_land_des = 0.6;
    vec2 ball_land_des;
    ball_land_des(X) = 0.0;
    ball_land_des(Y) = dist_to_table - 1*table_length/4;
    optim_des pred_params;
    pred_params.Nmax = 1000;
    calc_racket_strategy(balls_pred,ball_land_des,time_land_des,pred_params);
    fp.set_return_time(1.0);
    fp.set_detach(false);
    fp.set_des_params(&pred_params);
    fp.set_verbose(false);
    fp.update_init_state(Q);
    fp.run();

    if (fp.check_update()) {
        ran_optim = true;
        read_optim = false;
    }
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
    vec_str files = get_files(home + "/table-tennis/json/");
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
    // remove newline from each file
    for (std::string & file : files) {
        file.erase(std::remove(file.begin(),file.end(),'\n'),file.end());
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
