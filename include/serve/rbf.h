/*
 * rbf.h
 *
 *  Created on: Sep 1, 2018
 *      Author: robolab
 */

#ifndef INCLUDE_SERVE_RBF_H_
#define INCLUDE_SERVE_RBF_H_

using arma::vec;
using arma::zeros;
using arma::mat;
using arma::vec7;

namespace serve {

class RBF { // rbf for each dof independent

private:

    double t = 0.0;
    double Tmax = 1.0;
    vec intercepts = zeros<vec>(7);
    vec widths = zeros<vec>(10);
    vec centers = zeros<vec>(10);
    mat params = zeros<mat>(10,7);

    void basis_fnc(const double T, optim::joint &) const;

public:

    RBF() {};
    RBF(const std::string & param_file);

    void reset() {t = 0.0; };
    void step(const double & dt, optim::joint &);
    void get_init_pos(vec7 & pos) const;
    void set_init_pos(const vec7 & pos);
    void get_goal_pos(vec7 & pos) const;

};

class CRBF { // coupled RBF where the dofs are dependent

private:

    double t = 0.0;
    double Tmax = 1.0;
    vec intercepts = zeros<vec>(7);
    mat widths = zeros<mat>(1,7);
    mat centers = zeros<mat>(1,7);
    vec params = zeros<vec>(1);

    // stacked basis function
    void basis_fnc(const double T, optim::joint &) const;

public:

    CRBF() {};
    CRBF(const std::string & param_file);

    void reset() {t = 0.0; };
    void step(const double & dt, optim::joint &);
    void get_init_pos(vec7 & pos) const;
    void set_init_pos(const vec7 & pos);
    void get_goal_pos(vec7 & pos) const;
};

}

#endif /* INCLUDE_SERVE_RBF_H_ */
