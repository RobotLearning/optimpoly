/*
 * Dynamic Movement Primitives
 *
 * @author Okan Koc
 * @update Fri, Jun 8, 2018
 */

#include <armadillo>
#include "dmp.h"
#include "json.hpp"
#include "constants.h"

using namespace arma;
using json = nlohmann::json;

void Canonical::reset() {
    phase = 1.0;
}

void Canonical::step(const double & dt) {
    phase -= dt * a * tau * phase;
}

DMP::DMP(std::vector<double> weights_,
        std::vector<double> centers_,
        std::vector<double> heights_,
        double alpha_,
        double beta_,
        double goal_,
        double init_pos) {

    set_weights(weights_,centers_,heights_);
    set_goal_state(goal_);
    alpha = alpha_;
    beta = beta_;
    set_init_state(init_pos);
    reset();
}

void DMP::reset() {
    x = x0;
}

void DMP::set_weights(const std::vector<double> & w_stdvec,
                     const std::vector<double> & c_stdvec,
                     const std::vector<double> & h_stdvec) {
    w = conv_to<vec>::from(w_stdvec);
    c = conv_to<vec>::from(c_stdvec);
    h = conv_to<vec>::from(h_stdvec);
}

void DMP::set_goal_state(const double & goal) {
    g = goal;
}

void DMP::set_init_state(const double & init_pos) {
    x0(0) = init_pos;
    x0(1) = 0.0;
    x0(2) = 0.0;
}

void DMP::get_goal_state(double & goal) const {
    goal = g;
}

void DMP::get_init_state(double & init_pos) const {
    init_pos = x0(0);
}

vec3 DMP::step(const Canonical & can, const double & dt) {

    const double amp = 1.0;
    double f = forcing(can.phase);
    mat22 A = {{0,1}, {-alpha*beta*can.tau*can.tau, -alpha*can.tau}};
    vec2 b = {0,can.tau*can.tau*(alpha*beta*g + amp*f)};
    vec2 xdot = A * x.head(2) + b;
    x(0) += dt * xdot(0);
    x(1) += dt * xdot(1);
    x(2) = xdot(1);
    return x;
}


vec DMP::basis(const double & x) const {
    return exp(-h % ((x-c) % (x-c)));
}

double DMP::forcing(const double & phase) const {
    double f = 0.0;
    double scale = 0.0;
    vec psi = basis(phase);
    for (unsigned int i = 0; i < w.n_elem; i++) {
        f += psi(i) * w(i) * phase;
        scale += psi(i);
    }
    f /= scale + w.n_elem * 1e-10;
    return f;
}

void Joint_DMPs::reset() {
    can.reset();
    for (unsigned int i = 0; i < dmps.size(); i++) {
        dmps[i].reset();
    }
}

Joint_DMPs::Joint_DMPs(const std::string & param_file) {

    using std::vector;
    // load from json file
    std::ifstream stream(param_file);
    json jobs;
    stream >> jobs;
    can.tau = jobs.at("tau");
    for (auto elem : jobs.at("joints")) {
        DMP dmp = DMP(elem.at("weights"),
                    elem.at("centers"),
                    elem.at("heights"),
                    jobs.at("alpha"),
                    jobs.at("beta"),
                    elem.at("goal"),
                    elem.at("init_pos"));
        dmps.push_back(dmp);

    }
}

vec Joint_DMPs::step(const double & dt) {

    vec q_next = zeros<vec>(dmps.size());
    for (unsigned int i = 0; i < dmps.size(); i++) {
        vec x = dmps[i].step(can,dt);
        q_next(i) = x(0);
    }
    can.step(dt);
    return q_next;
}

mat Joint_DMPs::evolve(const double & T) {

    using namespace const_tt;
    unsigned int N = T/DT;
    reset();
    mat joints = zeros<mat>(NDOF,N);
    for (unsigned int i = 0; i < N; i++) {
        joints.col(i) = step(DT);
    }
    reset();
    return joints;

}

void Joint_DMPs::get_init_pos(vec & pos) const {

    using namespace std;
    try {
        for (unsigned int i = 0; i < dmps.size(); i++) {
            dmps[i].get_init_state(pos(i));
        }
    }
    catch (std::exception & ex) {
        cerr << "Array length incorrect: " << ex.what() << endl;
    }
}

void Joint_DMPs::set_init_pos(const vec & pos) {

    using namespace std;
    try {
        for (unsigned int i = 0; i < dmps.size(); i++) {
            dmps[i].set_init_state(pos(i));
        }
    }
    catch (std::exception & ex) {
        cerr << "Array length incorrect: " << ex.what() << endl;
    }
}

void Joint_DMPs::set_goal_pos(const vec & pos) {

    using namespace std;
    try {
        for (unsigned int i = 0; i < dmps.size(); i++) {
            dmps[i].set_goal_state(pos(i));
        }
    }
    catch (std::exception & ex) {
        cerr << "Array length incorrect: " << ex.what() << endl;
    }
}

void Joint_DMPs::get_goal_pos(vec & pos) const {

    using namespace std;
    try {
        for (unsigned int i = 0; i < dmps.size(); i++) {
            dmps[i].get_goal_state(pos(i));
        }
    }
    catch (std::exception & ex) {
        cerr << "Array length incorrect: " << ex.what() << endl;
    }
}

double Joint_DMPs::get_time_constant() const {
    return can.tau;
}
