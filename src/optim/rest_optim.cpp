/**
 * @file rest_optim.cpp
 *
 * @brief This file contains the resting state optimization routines.
 *
 *  Created on: Jul 31, 2017
 *      Author: okoc
 */

#include <armadillo>
#include "optim.h"
#include "utils.h"
#include "kinematics.hpp"

using namespace arma;
using namespace const_tt;

namespace optim {

/*
 * Interpolate ball positions at a given time point T.
 */
static void interp_ball(const mat & ballpred,
                        const double T,
                        vec3 & ballpos);

/*
 * Make sure that the robot resting posture touches the ball path at some point.
 */
static void intersect_ball_path(unsigned m,
                                double *result,
                                unsigned n,
                                const double *x,
                                double *grad,
                                void *data);

/*
 * Cost function for the resting posture optim
 */
static double cost_fnc(unsigned n,
                        const double *x,
		                double *grad,
		                void *data);

void Optim::run_qrest_optim(vec7 & q_rest_des) {
	std::thread t = std::thread(&Optim::optim_rest_posture,this,std::ref(q_rest_des));
	if (detach) {
		t.detach();
	}
	else {
		t.join();
	}
}

void Optim::optim_rest_posture(vec7 & q_rest_des) {

	double x[NDOF+1];
	double tol_eq[NCART];
	const_vec(NCART,1e-2,tol_eq);
	rest_optim_data rest_data;
	rest_data.ball_pred = param_des->ball_pos;
	double lb_[NDOF+1], ub_[NDOF+1];
	for (int i = 0; i < NDOF; i++) {
		rest_data.q_hit(i) = qf[i];
		lb_[i] = lb[i];
		ub_[i] = ub[i];
		x[i] = qf[i];
	}
	lb_[NDOF] = lb[2*NDOF];
	ub_[NDOF] = ub[2*NDOF];
	x[NDOF] = 1.0; // time estimate

	/*for (int i = 0; i < NDOF; i++) {
		printf("q_rest_old[%d] = %f\n", i, x[i]);
	}*/

	nlopt_opt opt = nlopt_create(NLOPT_LN_COBYLA, NDOF+1);
	nlopt_set_xtol_rel(opt, 1e-2);
	nlopt_set_lower_bounds(opt, lb_);
	nlopt_set_upper_bounds(opt, ub_);
	nlopt_set_min_objective(opt, cost_fnc, (void*)&rest_data);
	nlopt_add_inequality_mconstraint(opt, NCART, intersect_ball_path, (void*)&rest_data, tol_eq);

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
			//printf("Time of ball predicted: %f\n", x[NDOF]);
		}
		for (int i = 0; i < NDOF; i++) {
			//printf("q_rest_new[%d] = %f\n", i, x[i]);
			q_rest_des(i) = x[i];
		}
		update_rest_state(q_rest_des);
		/*mat::fixed<6,7> jac = zeros<mat>(6,7);
		vec3 pos = get_jacobian(q,jac);
		cout << "cart_pos:\n" << pos;
		cout << "jac:\n" << jac;*/
	}

	nlopt_destroy(opt);
}

static double cost_fnc(unsigned n,
                        const double *x,
		                double *grad,
		                void *data) {

	rest_optim_data *rest_data = (rest_optim_data*)data;
	static mat::fixed<6,7> jac = zeros<mat>(6,7);
	vec q_rest(x,NDOF);
	player::get_jacobian(q_rest,jac);

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

	double move_cost = pow(norm(q_rest - rest_data->q_hit),2);
	return pow(norm(jac,"fro"),2) +  + move_cost;
}

static void intersect_ball_path(unsigned m,
                                double *result,
                                unsigned n,
                                const double *x,
                                double *grad,
                                void *data) {

	rest_optim_data *rest_data = (rest_optim_data*)data;
	vec3 ball_pos;
	vec q_rest(x,NDOF);
	double T = x[NDOF];
	static mat::fixed<6,7> jac = zeros<mat>(6,7);
	vec3 robot_pos = player::get_jacobian(q_rest,jac);
	interp_ball(rest_data->ball_pred,T,ball_pos);

	for (int i = 0; i < NCART; i++) {
		result[i] = robot_pos(i) - ball_pos(i);
	}
}

static void interp_ball(const mat & ballpred,
                        const double T,
                        vec3 & ballpos) {

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

}
