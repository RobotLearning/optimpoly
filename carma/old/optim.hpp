/*
 * optim.hpp
 *
 *  Created on: Feb 12, 2017
 *      Author: okoc
 */

#ifndef OPTIM_HPP_
#define OPTIM_HPP_

#include <armadillo>
#include "nlopt.hpp"

static const bool PRINT_VERBOSE = true;
static const bool DEBUG_SOLN = true; // seperate source of verbosity for testing
static const int EQ_CONSTR_DIM = 3*NCART; // dimension of equality constraints
static const int INEQ_CONSTR_DIM = 2*NDOF + 2*NDOF; // strike and return, min and max
static const int OPTIM_DIM = 2*NDOF + 1; // number of parameters in the optimization
static const double MAX_VEL = 200.0; // maximum allowed velocity
static const double SLACK = 0.10; // slack on the joint limits
static const double TOL_CONSTR = 0.0001; // relative constraint tolerances
static const double TOL_REL_OPTIM = 0.01; // relative optim parameter tolerances

using namespace std;
using namespace arma;
using namespace nlopt;

typedef std::vector<double> stdvec;

typedef struct {
	vec7 a3; // first poly coeffs
	vec7 a2; // second poly coeffs
	vector<double> x_last; // opt params to compare to
} poly;

class PolyOptim { // public static class for running optimizations

public:

	static optim init;
	static optim rest;
	static optim max;
	static optim min;
	static optim* strike; // this will be updated

	static racket racket_des; // Player class provides racket info

	static poly coeff_str; // striking polynomial coeffs
	static poly coeff_ret; // returning polynomial coeffs

	opt minim; // optimizer

	int count; // number of optimizations so far
	wall_clock timer;

	// equality constraints
	static void kinematics_eq_constr(unsigned m, double *result, unsigned n,
			            const double *x, double *grad, void *my_func_params);
	static void interp(const double &T, vec3 & racket_pos,
			    vec3 & racket_vel, vec3 & racket_n);
	void finalize_soln(const vector<double> & x);
	void nlopt_optim_poly_run();
	void init_soln(vector<double> & x) const;

	static double const_costfunc(const vector<double> &x,
			vector<double> & grad, void *my_func_params);
	static double costfunc(const vector<double> &x,
			vector<double> &grad, void* my_func_params);

	// inequality constraint functions
	static void joint_lim_ineq_constr(unsigned m, double *result, unsigned n,
			            const double *x, double *grad, void *my_func_params);
	static void calc_strike_poly_coeff(const vector<double> & x);
	static void calc_return_poly_coeff(const vector<double> & x);
	static void calc_strike_ext_cand(const double T,
			vec7 & joint_max_cand, vec7 & joint_min_cand);
	static void calc_return_ext_cand(const vector<double> &x,
			vec7 & joint_max_cand, vec7 & joint_min_cand);

	// Constructor - simple
	PolyOptim(const vec7 & q0, const racket & des_racket,
			             const double & time2return, optim & params0,
						 opt & optimizer);
	// MPC constructor
	PolyOptim(const optim & initial, const optim & resting,
			 const racket & des_racket, optim & params0, opt & optimizer);

	void set_bounds();
	void setup();
	void operator()();

};

bool check_optim_result(const result & res);
void print_optim_vec(const vector<double> & x);


#endif /* OPTIM_HPP_ */
