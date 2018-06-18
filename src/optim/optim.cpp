/*
 * optim.cpp
 *
 *  Created on: Apr 9, 2018
 *      Author: okoc
 */

#include <armadillo>
#include "constants.h"
#include "utils.h"
#include "stdlib.h"
#include "math.h"
#include "kinematics.h"
#include "optim.h"
#include "tabletennis.h"
#include "lookup.h"

namespace optim {

/*
 * Give info about the optimization after termination
 */
static bool check_optim_result(const int res);

Optim::~Optim() {

    nlopt_destroy(opt);
}

void Optim::update_init_state(const joint & qact) {
    for (int i = 0; i < NDOF; i++) {
        q0[i] = qact.q(i);
        q0dot[i] = qact.qd(i);
    }
}

bool Optim::check_running() {
    return running;
}

bool Optim::check_update() {
    return update;
}

void Optim::set_moving(bool flag_move) {
    moving = flag_move;
}

void Optim::set_detach(bool flag_detach) {
    detach = flag_detach;
}

void Optim::set_return_time(const double & ret_time) {
    time2return = ret_time;
}

void Optim::set_verbose(bool flag_verbose) {
    verbose = flag_verbose;
}

bool Optim::get_params(const joint & qact, spline_params & p) {

    bool flag = false;
    if (update && !running) {
        vec7 qf_, qfdot_, qrest_;
        for (int i = 0; i < NDOF; i++) {
            qf_(i) = qf[i];
            qfdot_(i) = qfdot[i];
            qrest_(i) = qrest[i];
        }
        vec7 qnow = qact.q;
        vec7 qdnow = qact.qd;
        p.a.col(0) = 2.0 * (qnow - qf_) / pow(T,3) + (qfdot_ + qdnow) / pow(T,2);
        p.a.col(1) = 3.0 * (qf_ - qnow) / pow(T,2) - (qfdot_ + 2.0*qdnow) / T;
        p.a.col(2) = qdnow;
        p.a.col(3) = qnow;
        //cout << "A = \n" << p.a << endl;
        p.b.col(0) = 2.0 * (qf_ - qrest_) / pow(time2return,3) + (qfdot_) / pow(time2return,2);
        p.b.col(1) = 3.0 * (qrest_ - qf_) / pow(time2return,2) - (2.0*qfdot_) / time2return;
        p.b.col(2) = qfdot_;
        p.b.col(3) = qf_;
        p.time2hit = T;
        //cout << "B = \n" << p.b << endl;
        flag = true;
        update = false;
    }
    return flag;
}

void Optim::update_rest_state(const vec7 & q_rest_new) {

    for (int i = 0; i < NDOF; i++)
        qrest[i] = q_rest_new(i);
}

void Optim::set_des_params(optim_des *params_) {
    param_des = params_;
}

void Optim::init_lookup_soln(double *x) {

    vec::fixed<OPTIM_DIM> robot_params;
    vec6 ball_params;
    for (int i = 0; i < NCART; i++) {
        ball_params(i) = param_des->ball_pos(i,0);
        ball_params(i+NCART) = param_des->ball_vel(i,0);
    }
    //cout << "Init ball est:" << ball_params << endl;
    player::predict_till_net(ball_params);
    //cout << "Net ball est:" << ball_params << endl;
    // k = 5 nearest neighbour regression
    player::knn(lookup_table,ball_params,5,robot_params);
    for (int i = 0; i < OPTIM_DIM; i++) {
        x[i] = robot_params(i);
        //  printf("x[%d] = %f\n", i, x[i]);
    }
}

void Optim::run() {

    std::thread t = std::thread(&Optim::optim,this);
    if (detach) {
        t.detach();
    }
    else {
        t.join();
    }
}

void Optim::optim() {

    update = false;
    running = true;
    double x[OPTIM_DIM];

    if (moving) {
        init_last_soln(x);
    }
    else {
        if (lookup) {
            if (verbose) {
                std::cout << "Looking up good initial parameters with k = 5\n"; // kNN parameter k = 5
            }
            init_lookup_soln(x);
        }
        else {
            init_rest_soln(x);
        }
    }

    double init_time = get_time();
    double past_time = 0.0;
    double minf; // the minimum objective value, upon return //
    int res; // error code

    if ((res = nlopt_optimize(opt, x, &minf)) < 0) {
        past_time = (get_time() - init_time)/1e3;
        if (verbose) {
            printf("NLOPT failed with exit code %d!\n", res);
            printf("NLOPT took %f ms\n", past_time);
        }
    }
    else {
        past_time = (get_time() - init_time)/1e3;
        if (verbose) {
            printf("NLOPT success with exit code %d!\n", res);
            printf("NLOPT took %f ms\n", past_time);
            printf("Found minimum at f = %0.10g\n", minf);
        }
        if (test_soln(x) < 1e-2)
            finalize_soln(x,past_time);
    }
    if (verbose)
        check_optim_result(res);
    running = false;
}

static bool check_optim_result(const int res) {

    bool flag = false;
    switch (res) {
    case NLOPT_SUCCESS:
        printf("Success!\n");
        flag = true;
        break;
    case NLOPT_STOPVAL_REACHED:
        printf("Optimization stopped because stopval (above) was reached.\n");
        flag = true;
        break;
    case NLOPT_FTOL_REACHED:
        printf("Optimization stopped because ftol_rel "
                "or ftol_abs (above) was reached.\n");
        flag = true;
        break;
    case NLOPT_XTOL_REACHED:
        flag = true;
        printf("Optimization stopped because xtol_rel or xtol_abs (above) was reached.\n");
        break;
    case NLOPT_MAXEVAL_REACHED:
        flag = true;
        printf("Optimization stopped because maxeval (above) was reached.\n");
        break;
    case NLOPT_MAXTIME_REACHED:
        flag = true;
        printf("Optimization stopped because maxtime (above) was reached.\n");
        break;
    case NLOPT_FAILURE:
        printf("Epic fail!\n");
        break;
    case NLOPT_INVALID_ARGS:
        printf("Invalid arguments (e.g. lower bounds are bigger than "
                "upper bounds, an unknown algorithm was specified, etcetera).\n");
        break;
    case NLOPT_OUT_OF_MEMORY:
        printf("Ran out of memory!\n");
        break;
    case NLOPT_ROUNDOFF_LIMITED:
        printf("Halted because roundoff errors limited progress."
            "(In this case, the optimization still typically returns a useful result.\n");
        break;
    case NLOPT_FORCED_STOP:
        printf("Halted because of a forced termination: "
                "the user called nlopt_force_stop(opt)"
                "on the optimization’s nlopt_opt object "
                "opt from the user’s objective function or constraints.\n");
        break;

    }
    return flag;
}

}
