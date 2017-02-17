/*
 * optim.cpp
 *
 *  Created on: Feb 10, 2017
 *      Author: okoc
 */

#include <thread>
#include <armadillo>
#include "constants.h"
#include "nlopt.hpp"
#include "kinematics.hpp"
#include "player.hpp"
#include "optim.hpp"

using namespace std;
using namespace arma;
using namespace nlopt;

// static class members are initialized here
// these are used as globals in cost and constraint functions!
optim PolyOptim::init = {zeros<vec>(NDOF), zeros<vec>(NDOF), 0.0, false};
optim PolyOptim::rest = {zeros<vec>(NDOF), zeros<vec>(NDOF), 0.0, false};
optim PolyOptim::max = {zeros<vec>(NDOF), zeros<vec>(NDOF), 0.0, false};
optim PolyOptim::min = {zeros<vec>(NDOF), zeros<vec>(NDOF), 0.0, false};
optim* PolyOptim::strike = 0;
racket PolyOptim::racket_des = {zeros<mat>(1,1),zeros<mat>(1,1),zeros<mat>(1,1)};
poly PolyOptim::coeff_str = {zeros<vec>(NDOF), zeros<vec>(NDOF), vector<double>(1,1)};
poly PolyOptim::coeff_ret = {zeros<vec>(NDOF), zeros<vec>(NDOF), vector<double>(1,1)};

/*
 * Simple constructor
 * that is not suitable for MPC.
 *
 * It assumes that the resting and initial states are the same.
 */
PolyOptim::PolyOptim(const vec7 & q0, const racket & des_racket,
		             const double & time2return, optim & params0,
					 opt & optimizer) : minim(optimizer) {

	optim initial = {q0, zeros<vec>(NDOF), 0.0, false};
	optim resting = {q0, zeros<vec>(NDOF), time2return, false};
	init = initial;
	rest = resting;
	racket_des = des_racket;
	strike = &params0;
	count = 0;
}

/*
 * Constructor for the polynomial optimization. Useful for MPC.
 *
 * Sets up the parameters and the optimization settings, optimizer and the routine
 *
 *
 */
PolyOptim::PolyOptim(const optim & initial, const optim & resting,
		             const racket & des_racket, optim & params,
					 opt & optimizer) : minim(optimizer) {

	init = initial;
	rest = resting;
	strike = &params;
	racket_des = des_racket;
	count = 0;
}

/*
 * Initialize optimizer settings for NLOPT library
 *
 */
void PolyOptim::setup() {

	//set up the optimization
	stdvec tol_eq(EQ_CONSTR_DIM,TOL_CONSTR);
	stdvec tol_ineq(INEQ_CONSTR_DIM,TOL_CONSTR);
	set_bounds();
	minim.set_min_objective(PolyOptim::costfunc,NULL);
	minim.add_inequality_mconstraint(PolyOptim::joint_lim_ineq_constr,NULL,tol_ineq);
	minim.add_equality_mconstraint(PolyOptim::kinematics_eq_constr,NULL,tol_eq);
	minim.set_xtol_rel(TOL_REL_OPTIM);
}

/*
 * Caller thread executes this function
 */
void PolyOptim::operator ()() {

	if (PRINT_VERBOSE) {
		cout << "==========================================" << endl;
		cout << "Running NLOPT" << endl;
	}

	timer.tic();
	nlopt_optim_poly_run();
	double time_elapsed = timer.toc();
	strike->t -= time_elapsed; // correct optimal hitting time for time elapsed computing

	if (PRINT_VERBOSE) {
		cout << "Optim count: " << ++count << endl;
		cout << "Optim time: " << time_elapsed << endl;
		cout << "==========================================" << endl;
	}
}


/*
 * NLOPT optimization routine for table tennis centred player (CP)
 *
 * If constraints are violated, it will not modify the lookup values (input)
 * TODO: give more info about failure!
 */
void PolyOptim::nlopt_optim_poly_run() {

	vector<double> x(OPTIM_DIM);
	//parameters are the initial joint positions q0
	init_soln(x);
	double minf; // the minimum objective value, upon return
	result res;
	try {
		res = minim.optimize(x, minf);
		//printf("Found minimum at f = %0.10g\n", minf);
		//if (test_optim(x))
		if (DEBUG_SOLN) {
			print_optim_vec(x);
		}
		if (check_optim_result(res))
			finalize_soln(x);
	}
	catch (const runtime_error & err) {
		cout << "NLOPT failed!\n";
		//freeze();
	}

}

/*
 * Finalize the solution and update target SL structure and hitTime value
 */
void PolyOptim::finalize_soln(const vector<double> & x) {

	vec p = conv_to<colvec>::from(x);
	strike->q = p(span(0,NDOF-1));
	strike->qdot = p(span(NDOF,2*NDOF-1));
	strike->t = p(2*NDOF);
	strike->update = true;
}

/*
 * Set the initial solution for NLOPT 2*dof + 1 dimensional problem
 *
 * (usually to initial posture with zero velocity)
 *
 *
 */
void PolyOptim::init_soln(vector<double> & x) const {

	// initialize first dof entries to q0
	for (int i = 0; i < NDOF; i++) {
		x[i] = rest.q(i);
		x[i + NDOF] = rest.qdot(i);
	}
	x[2*NDOF] = rest.t;
}

/*
 * Constant cost function.
 * Used to project lookup table solutions to equality constraints.
 *
 * Note: seems to be much better than costfunc when the optim is restricted to 1ms.
 */
double PolyOptim::const_costfunc(const vector<double> &x,
		vector<double> & grad, void *my_func_params) {

	int i;
	if (!grad.empty()) {
		std::fill(grad.begin(),grad.end(),0);
	}

	return 1;
}

/*
 * Calculates the cost function for table tennis centred player
 * to find spline (3rd order strike+return) polynomials
 */
double PolyOptim::costfunc(const vector<double> &x,
		vector<double> &grad, void* my_func_params) {

	double T = x[2*NDOF];
	// calculate the polynomial coeffs which are used in the cost calculation
	calc_strike_poly_coeff(x);

	return T * (3*T*T*dot(coeff_str.a3,coeff_str.a3) +
			3*T*dot(coeff_str.a3,coeff_str.a2) +
			dot(coeff_str.a2,coeff_str.a2));
}



/*
 * Set upper and lower bounds on the optimization
 *
 * Sets as well the max and min fields (to be used in optimization)
 */
void PolyOptim::set_bounds() {

	vector<double> ub(OPTIM_DIM);
	vector<double> lb(OPTIM_DIM);
	double upp[NDOF] = {0.0};
	double low[NDOF] = {0.0};
	read_joint_limits(low,upp);
	// lower bounds and upper bounds for qf are the joint limits
	for (int i = 0; i < NDOF; i++) {
		max.q(i) = ub[i] = upp[i] - SLACK;
		min.q(i) = lb[i] = low[i] - SLACK;
		max.qdot(i) = ub[i + NDOF] = MAX_VEL;
		min.qdot(i) = lb[i + NDOF] = -MAX_VEL;
	}
	// constraints on final time
	max.t = ub[2*NDOF] = 2.00;
	min.t = lb[2*NDOF] = 0.00;
	minim.set_lower_bounds(lb);
	minim.set_upper_bounds(ub);
}

/*
 * This is the inequality constraint that makes sure we never exceed the
 * joint limits during the striking and returning motion
 *
 */
void PolyOptim::joint_lim_ineq_constr(unsigned m, double *result, unsigned n,
		const double *x, double *grad, void *my_func_params) {

	vec7 str_max_cands;
	vec7 str_min_cands;
	vec7 ret_max_cands;
	vec7 ret_min_cands;

	// calculate the polynomial coeffs which are used for checking joint limits
	vector<double> x_(x,x + OPTIM_DIM);
	calc_strike_poly_coeff(x_);
	calc_return_poly_coeff(x_);
	// calculate the candidate extrema both for strike and return
	calc_strike_ext_cand(x[2*NDOF],str_max_cands,str_min_cands);
	calc_return_ext_cand(x_,ret_max_cands,ret_min_cands);

	vec res = join_vert(join_vert(str_max_cands - max.q,
			                     -(str_min_cands - min.q)),
			            join_vert(ret_max_cands - max.q,
					             -(ret_min_cands - min.q)));

	result = res.memptr();


}

/*
 * Calculate the polynomial coefficients from the optimized variables qf,qfdot,T
 * p(t) = a1*t^3 + a2*t^2 + a3*t + a4
 *
 *
 */
void PolyOptim::calc_strike_poly_coeff(const vector<double> & x) {

	if (x == coeff_str.x_last) {
		// dont do anything
	}
	else {
		vec xarma = conv_to<colvec>::from(x);
		// variables to be optimized
		double T = x[2*NDOF];
		vec7 qf = xarma(span(0,NDOF-1));
		vec7 qfdot = xarma(span(NDOF,2*NDOF-1));

		coeff_str.a3 = 2/pow(T,3) * (init.q - qf) + 1/pow(T,2) * (init.qdot + qfdot);
		coeff_str.a2 = 3/pow(T,2) * (qf - init.q) - (1/T) * (qfdot + 2*init.qdot);
		coeff_str.x_last = x;
	}

}

/*
 * Calculate the returning polynomial coefficients from the optimized variables qf,qfdot
 * and time2return variable defined in the header
 * p(t) = a1*t^3 + a2*t^2 + a3*t + a4
 */
void PolyOptim::calc_return_poly_coeff(const vector<double> & x) {

	if (x == coeff_ret.x_last) {
		// dont update anything!
	}
	else {
		vec xarma = conv_to<colvec>::from(x);
		// variables to be optimized
		double T = rest.t;
		vec7 qf = xarma(span(0,NDOF-1));
		vec7 qfdot = xarma(span(NDOF,2*NDOF-1));

		coeff_ret.a3 = 2/pow(T,3) * (qf - rest.q) + 1/pow(T,2) * (rest.qdot + qfdot);
		coeff_ret.a2 = 3/pow(T,2) * (rest.q - qf) - (1/T) * (rest.qdot + 2*qfdot);
		coeff_ret.x_last = x;
	}
}

/*
 * Calculate the extrema candidates for each joint (2*7 candidates in total)
 * For the striking polynomial
 * Clamp to [0,T]
 *
 */
void PolyOptim::calc_strike_ext_cand(const double T,
		vec7 & joint_max_cand, vec7 & joint_min_cand) {

	vec7 a3 = coeff_str.a3;
	vec7 a2 = coeff_str.a2;
	// find the extrema time candidates
	vec7 cand1 = clamp(-a2 + sqrt(a2 % a2 - 3 * a3 % init.qdot) / (3*a3), 0, T);
	vec7 cand2 = clamp(-a2 - sqrt(a2 % a2 - 3 * a3 % init.qdot) / (3*a3), 0, T);
	// find the joint extrema values at candidate times
	cand1 = a3 % pow(cand1,3) + a2 % pow(cand1,2) + init.qdot % cand1 + init.q;
	cand2 = a3 % pow(cand2,3) + a2 % pow(cand2,2) + init.qdot % cand2 + init.q;
	joint_max_cand = arma::max(cand1,cand2);
	joint_min_cand = arma::min(cand1,cand2);
}

/*
 * Calculate the extrema candidates for each joint (2*7 candidates in total)
 * For the return polynomial
 * Clamp to [0,TIME2RETURN]
 *
 */
void PolyOptim::calc_return_ext_cand(const vector<double> &x,
		vec7 & joint_max_cand, vec7 & joint_min_cand) {

	vec p = conv_to<colvec>::from(x);
	strike->q = p(span(0,NDOF-1));
	strike->qdot = p(span(NDOF,2*NDOF-1));
	strike->t = p(2*NDOF);
	vec7 a3 = coeff_ret.a3;
	vec7 a2 = coeff_ret.a2;
	double T = rest.t; // time2return;

	// find the extrema time candidates
	vec7 cand1 = clamp(-a2 + sqrt(a2 % a2 - 3 * a3 % strike->qdot) / (3*a3), 0, T);
	vec7 cand2 = clamp(-a2 - sqrt(a2 % a2 - 3 * a3 % strike->qdot) / (3*a3), 0, T);
	// find the joint extrema values at candidate times
	cand1 = a3 % pow(cand1,3) + a2 % pow(cand1,2) + strike->qdot % cand1 + strike->q;
	cand2 = a3 % pow(cand2,3) + a2 % pow(cand2,2) + strike->qdot % cand2 + strike->q;
	joint_max_cand = arma::max(cand1,cand2);
	joint_min_cand = arma::min(cand1,cand2);
}

/*
 * This is the constraint that makes sure we hit the ball and hit it well!
 *
 */
void PolyOptim::kinematics_eq_constr(unsigned m, double *result, unsigned n,
		const double *x, double *grad, void *my_func_params){

	vector<double> x_(x,x + OPTIM_DIM);
	vec p = conv_to<colvec>::from(x_);
	strike->q = p(span(0,NDOF-1));
	strike->qdot = p(span(NDOF,2*NDOF-1));
	strike->t = p(2*NDOF);
	vec3 racket_des_pos, racket_des_vel, racket_des_normal;
	vec3 racket_pos, racket_vel, racket_normal;
	// get desired racket parameters at new time T
	interp(strike->t,racket_des_pos,racket_des_vel,racket_des_normal);
	// compute actual racket pos,vel and normal
	calc_racket_state(strike->q,strike->qdot,racket_pos,racket_vel,racket_normal);
	// get the deviation
	vec res = join_vert(join_vert(racket_pos - racket_des_pos,
			                      racket_vel - racket_des_vel),
								  racket_normal - racket_des_normal);
	/*cout << "vec:" << endl;
	cout << res << endl;
	result = res.memptr();
	cout << "pointer:" << endl;
	for (int i = 0; i < 9; i++)
		cout << res[i] << "  ";
	cout << endl;*/

}

/*
 * First order hold to interpolate linearly at time T between ball prediction matrix ballMat entries
 *
 * TODO: global dt, Tpred and Tmax variables needed?
 *
 */
void PolyOptim::interp(const double & T, vec3 & racket_pos,
		               vec3 & racket_vel, vec3 & racket_n) {

	int i;
	const double dt = 0.002;
	int N = (int) (T/dt);
	double Tdiff = T - N*dt;
	static int iter;
	int Nmax = racket_des.pos.n_cols - 1;

	if (isnan(T)) {
		cout << "Warning: T value is nan!" << endl;
		racket_pos = racket_des.pos.col(0);
		racket_vel = racket_des.vel.col(0);
		racket_n = racket_des.normal.col(0);
	}
	else {
		if (N < Nmax) {
			racket_pos = racket_des.pos.col(N) + (Tdiff/dt) * (racket_des.pos.col(N+1) - racket_des.pos.col(N));
			racket_vel = racket_des.vel.col(N) + (Tdiff/dt) * (racket_des.vel.col(N+1) - racket_des.vel.col(N));
			racket_n = racket_des.normal.col(N) + (Tdiff/dt) * (racket_des.normal.col(N+1) - racket_des.normal.col(N));
		}
		else {
			racket_pos = racket_des.pos.col(N);
			racket_vel = racket_des.vel.col(N);
			racket_n = racket_des.normal.col(N);
		}
	}
}

/*
 * Prints the 2*NDOF + 1 dimensional solution in user-friendly format
 */
void print_optim_vec(const vector<double> & x) {

	int i;
	printf("qf = [");
	for (i = 0; i < NDOF; i++) {
		printf("%.2f  ", x[i]);
	}
	printf("]\n");
	printf("qfdot = [");
	for (i = 0; i < NDOF; i++) {
		printf("%.2f  ", x[i+NDOF]);
	}
	printf("]\n");
	printf("T = %.2f\n", x[2*NDOF]);
}

/*
 * Give info about the optimization
 * TODO: only prinf if DEBUG_SOLN flag is on
 *
 */
bool check_optim_result(const result & res) {

	bool flag = false;
	switch (res) {
	case NLOPT_SUCCESS:
		cout << "Success!" << endl;
		flag = true;
		break;
	case NLOPT_STOPVAL_REACHED:
		cout << "Optimization stopped because stopval (above) was reached." << endl;
		flag = true;
		break;
	case NLOPT_FTOL_REACHED:
		cout << "Optimization stopped because ftol_rel or ftol_abs (above) was reached." << endl;
		flag = true;
		break;
	case NLOPT_XTOL_REACHED:
		flag = true;
		cout << "Optimization stopped because xtol_rel or xtol_abs (above) was reached." << endl;
		break;
	case NLOPT_MAXEVAL_REACHED:
		flag = true;
		cout << "Optimization stopped because maxeval (above) was reached." << endl;
		break;
	case NLOPT_MAXTIME_REACHED:
		flag = true;
		cout << "Optimization stopped because maxtime (above) was reached." << endl;
		break;
	case NLOPT_FAILURE:
		cout << "Epic fail!" << endl;
		break;
	case NLOPT_INVALID_ARGS:
		cout << "Invalid arguments (e.g. lower bounds are bigger than "
				"upper bounds, an unknown algorithm was specified, etcetera)." << endl;
		break;
	case NLOPT_OUT_OF_MEMORY:
		cout << "Ran out of memory!" << endl;
		break;
	case NLOPT_ROUNDOFF_LIMITED:
		cout << "Halted because roundoff errors limited progress."
				"(In this case, the optimization still typically returns a useful result.)" << endl;
		break;
	case NLOPT_FORCED_STOP:
		cout << "Halted because of a forced termination: the user called nlopt_force_stop(opt)"
				"on the optimization’s nlopt_opt object opt from the user’s objective function or constraints."
		     << endl;
		break;

	}
	return flag;
}
