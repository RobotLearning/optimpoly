/*
 * test_optim.cpp
 *
 * Unit Tests for polynomial optimization
 *
 *  Created on: Feb 17, 2017
 *      Author: okoc
 */

#ifndef BOOST_TEST_MODULE
#define BOOST_TEST_MODULE test_table_tennis
#endif

#include <boost/test/unit_test.hpp>
#include <iostream>
#include <armadillo>
#include <thread>
#include "kinematics.h"
#include "utils.h"
#include "optim.h"
#include "lookup.h"
#include "player.hpp"
#include "tabletennis.h"
#include "kinematics.hpp"

using namespace arma;
int randval;

static double cost_fnc(unsigned n, const double *x,
		                   double *grad, void *data);
static void interp_ball(const mat & ballpred,
		                const double T, vec3 & ballpos);

/*
 * Initialize robot posture on the right size of the robot
 */
inline void init_right_posture(vec7 & q0) {

	q0(0) = 1.0;
	q0(1) = -0.2;
	q0(2) = -0.1;
	q0(3) = 1.8;
	q0(4) = -1.57;
	q0(5) = 0.1;
	q0(6) = 0.3;
}

/*
 *
 * Here testing NLOPT optimization for VHP player
 *
 */
BOOST_AUTO_TEST_CASE(test_vhp_optim) {

	cout << "Testing VHP Trajectory Optimizer...\n";
	const double VHPY = -0.3;
	double lb[2*NDOF+1], ub[2*NDOF+1];
	double SLACK = 0.01;
	double Tmax = 1.0;
	joint qact;
	spline_params poly;

	// update initial parameters from lookup table
	std::cout << "Looking up a random ball entry..." << std::endl;
	arma_rng::set_seed_random();
	randval = (randi(1).at(0));
	arma_rng::set_seed(randval);
	vec::fixed<15> strike_params;
	vec6 ball_state;
	lookup_random_entry(ball_state,strike_params);
	//std::cout << ball_state << std::endl;
	init_right_posture(qact.q);
	set_bounds(lb,ub,SLACK,Tmax);

	EKF filter = init_filter();
	mat66 P; P.eye();
	filter.set_prior(ball_state,P);

	double time_pred;
	vec6 ball_pred;
	game game_state = AWAITING;
	vec2 ball_land_des = {0.0, dist_to_table - 3*table_length/4};
	double time_land_des = 0.8;
	BOOST_TEST(predict_hitting_point(VHPY,true,ball_pred,time_pred,filter,game_state));
	//cout << ball_pred << endl;
	optim_des racket_params;
	calc_racket_strategy(ball_pred,ball_land_des,time_land_des,racket_params);

	vec3 normal_example = racket_params.racket_normal(span(X,Z),0);
	BOOST_TEST(arma::norm(normal_example) == 1.0, boost::test_tools::tolerance(0.01));

	HittingPlane opt = HittingPlane(qact.q.memptr(),lb,ub);
	opt.set_des_params(&racket_params);
	opt.fix_hitting_time(time_pred);
	opt.update_init_state(qact);
	opt.run();
	bool update = opt.get_params(qact,poly);

	BOOST_TEST(update);
}

/*
 * Testing Fixed Player (or Focused Player)
 */
BOOST_AUTO_TEST_CASE(test_fp_optim) {

	cout << "Testing FP Trajectory Optimizer...\n";
	double lb[2*NDOF+1], ub[2*NDOF+1];
	double SLACK = 0.01;
	double Tmax = 1.0;
	joint qact;
	spline_params poly;

	// update initial parameters from lookup table
	std::cout << "Looking up a random ball entry..." << std::endl;
	arma_rng::set_seed(randval);
	//arma_rng::set_seed_random();
	vec::fixed<15> strike_params;
	vec6 ball_state;
	lookup_random_entry(ball_state,strike_params);
	init_right_posture(qact.q);
	set_bounds(lb,ub,SLACK,Tmax);
	optim_des racket_params;
	int N = 1000;
	racket_params.Nmax = 1000;

	EKF filter = init_filter();
	mat66 P; P.eye();
	filter.set_prior(ball_state,P);

	vec6 ball_pred;
	double time_land_des = 0.8;
	mat balls_pred = filter.predict_path(DT,N);
	vec2 ball_land_des = {0.0, dist_to_table - 3*table_length/4};
	racket_params = calc_racket_strategy(balls_pred,ball_land_des,time_land_des,racket_params);

	FocusedOptim opt = FocusedOptim(qact.q.memptr(),lb,ub);
	opt.set_des_params(&racket_params);
	opt.update_init_state(qact);
	opt.run();
	bool update = opt.get_params(qact,poly);

	BOOST_TEST(update);
}

/*
 * Testing Lazy Player (or Defensive Player)
 */
BOOST_AUTO_TEST_CASE(test_dp_optim) {

	cout << "Testing LAZY Trajectory Optimizer...\n";
	double lb[2*NDOF+1], ub[2*NDOF+1];
	double SLACK = 0.01;
	double Tmax = 1.0;
	joint qact;
	spline_params poly;

	// update initial parameters from lookup table
	std::cout << "Looking up a random ball entry..." << std::endl;
	arma_rng::set_seed(randval);
	//arma_rng::set_seed_random();
	vec::fixed<15> strike_params;
	vec6 ball_state;
	lookup_random_entry(ball_state,strike_params);
	init_right_posture(qact.q);
	set_bounds(lb,ub,SLACK,Tmax);

	int N = 1000;
	EKF filter = init_filter();
	mat66 P; P.eye();
	filter.set_prior(ball_state,P);
	mat balls_pred = filter.predict_path(DT,N);
	optim_des ball_params;
	ball_params.ball_pos = balls_pred.rows(X,Z);
	ball_params.ball_vel = balls_pred.rows(DX,DZ);
	ball_params.Nmax = N;
	bool land = true;
	bool lookup = true;
	LazyOptim opt = LazyOptim(qact.q.memptr(),lb,ub,land,lookup); //only touch the ball if false!
	opt.set_des_params(&ball_params);
	opt.update_init_state(qact);
	opt.set_verbose(true);
	opt.run();
	bool update = opt.get_params(qact,poly);

	BOOST_TEST(update);
}

struct rest_optim_data {
	mat ball_pred;
	vec7 q_hit;
};

/*
 * Find a qf and t such that J(qf) has minimal Frobenius norm
 * while being close to predicted ball state b(t) and being close to q_hit
 */
BOOST_AUTO_TEST_CASE(find_rest_posture) {

	cout << "\nOptimizing resting posture based on jacobian...\n";

	double lb[NDOF+1], ub[NDOF+1];
	double x[NDOF] = {0.0};
	for(int i = 0; i < NDOF; i++)
		x[i] = as_scalar(randn(1)); // joint estimates
	x[NDOF] = 1.0; // time estimate

	// update initial parameters from lookup table
	std::cout << "Looking up a random ball entry..." << std::endl;
	arma_rng::set_seed_random();
	vec::fixed<15> strike_params;
	vec6 ball_state;
	lookup_random_entry(ball_state,strike_params);
	vec7 q_hit = strike_params.head(NDOF);
	int N = 1000;
	EKF filter = init_filter();
	mat66 P; P.eye();
	filter.set_prior(ball_state,P);
	mat balls_pred = filter.predict_path(DT,N);
	rest_optim_data rest_data;
	rest_data.ball_pred = balls_pred;
	rest_data.q_hit = q_hit;
	read_joint_limits(lb,ub);
	lb[NDOF] = 0.0;
	ub[NDOF] = 2.0;
	nlopt_opt opt = nlopt_create(NLOPT_LN_COBYLA, NDOF+1);
	nlopt_set_xtol_rel(opt, 1e-2);
	nlopt_set_lower_bounds(opt, lb);
	nlopt_set_upper_bounds(opt, ub);
	nlopt_set_min_objective(opt, cost_fnc, (void*)&rest_data);

	double init_time = get_time();
	double past_time = 0.0;
	double minf; // the minimum objective value, upon return //
	int res; // error code

	if ((res = nlopt_optimize(opt, x, &minf)) < 0) {
		printf("NLOPT failed with exit code %d!\n", res);
	    past_time = (get_time() - init_time)/1e3;
		printf("NLOPT took %f ms\n", past_time);
	}
	else {
		past_time = (get_time() - init_time)/1e3;
		printf("NLOPT success with exit code %d!\n", res);
		printf("NLOPT took %f ms\n", past_time);
		printf("Found minimum at f = %0.10g\n", minf);
		for (int i = 0; i < NDOF; i++) {
			printf("q_rest[%d] = %f\n", i, x[i]);
		}
		printf("Time: %f\n", x[NDOF]);
		/*mat::fixed<6,7> jac = zeros<mat>(6,7);
		vec3 pos = get_jacobian(q,jac);
		cout << "cart_pos:\n" << pos;
		cout << "jac:\n" << jac;*/
	}

}

static double cost_fnc(unsigned n, const double *x,
		                   double *grad, void *data) {

	rest_optim_data *rest_data = (rest_optim_data*)data;
	static mat::fixed<6,7> jac = zeros<mat>(6,7);
	vec q_rest(x,NDOF);
	double T = x[NDOF];
	vec3 robot_pos = get_jacobian(q_rest,jac);
	vec3 ball_pos;

	if (grad) {
		static double h = 1e-6;
		static double val_plus, val_minus;
		static double xx[NDOF+1];
		for (unsigned i = 0; i < n; i++)
			xx[i] = x[i];
		for (unsigned i = 0; i < n; i++) {
			xx[i] += h;
			val_plus = cost_fnc(n, xx, NULL, data);
			xx[i] -= 2*h;
			val_minus = cost_fnc(n, xx, NULL, data);
			grad[i] = (val_plus - val_minus) / (2*h);
			xx[i] += h;
		}
	}

	interp_ball(rest_data->ball_pred,T,ball_pos);
	double move_cost = 12 * pow(norm(q_rest - rest_data->q_hit),2);
	return pow(norm(jac,"fro"),2) + pow(norm(robot_pos - ball_pos),2) + move_cost;
}

static void interp_ball(const mat & ballpred, const double T, vec3 & ballpos) {

    const double dt = DT;
	if (std::isnan(T)) {
		printf("Warning: T value is nan!\n");
		ballpos = ballpred.col(0).head(3);
	}
	else {
		unsigned N = (int) (T/dt);
		double Tdiff = T - N*dt;
		ballpos = ballpred.col(N).head(3) +
				(Tdiff/dt) * (ballpred.col(N+1).head(3) - ballpred.col(N).head(3));
	}
}


