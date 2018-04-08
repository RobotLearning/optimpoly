/**
 * @file racket_optim.cpp
 *
 * @brief Boundary Value Problem for computing desired outgoing ball velocities
 * is solved with spin
 *
 * We're mostly interested in the solution speed (time elapsed) using optimization
 * since an argument needs to be added to the paper
 *
 *
 *  Created on: Jan 7, 2018
 *      Author: okoc
 */

#include <stdio.h>
#include <stdlib.h>
#include "sys/time.h"

// optimization and math libraries
#include <math.h>
#include <nlopt.h>
#include <armadillo>
// table tennis prediction functions
#include "tabletennis.h"
#include "player.hpp"
#include "utils.h"

namespace optim {

/**
 * @brief Incoming ball state (+topspin) and desired landing point and time
 */
struct des_ball_data {
	vec3 ball_incoming;
	vec3 ball_land_des;
	double time_land_des;
	double topspin;
};

/*
 * Cost function for computing the residual (norm squared)
 * of the outgoing ball landing error
 * Calculates also the gradient if grad is TRUE
 *
 */
static double calc_landing_res(unsigned n,
                               const double *x,
                               double *grad,
                               void *data);

/*
 * Solve BVP for a particular predicted spinning ball's outgoing desired ball velocity
 *
 * BVP is solved using optimization
 */
static void optim_spin_outgoing_ball_vel(const des_ball_data & data,
                                         const bool verbose,
                                         vec3 & est); // spin based optimization

optim_des calc_spin_racket_strategy(const mat & balls_predicted,
								    const double & topspin,
								    const vec3 & ball_land_des,
								    const double time_land_des,
								    optim_des & racket_params) {

	TableTennis tennis = TableTennis(true,false);
	tennis.set_topspin(topspin);
	des_ball_data data;
	data.ball_land_des = ball_land_des;
	data.time_land_des = time_land_des;
	data.topspin = topspin;

	int N = balls_predicted.n_cols;
	mat balls_out_vel = zeros<mat>(3,N);
	mat racket_des_normal = zeros<mat>(3,N);
	mat racket_des_vel = zeros<mat>(3,N);
	vec3 vel_out;

	// get initial outgoing ball velocity estimates
	tennis.calc_des_ball_out_vel(ball_land_des.head(2),time_land_des,false,balls_predicted,balls_out_vel);
	vel_out = balls_out_vel.col(0);
	// refine the estimates with optimization
	for (int i = 0; i < N; i++) {
		data.ball_incoming = balls_predicted.col(i).head(3);
		optim_spin_outgoing_ball_vel(data,false,vel_out);
		balls_out_vel.col(i) = vel_out;
	}
	tennis.calc_des_racket_normal(balls_predicted.rows(DX,DZ),balls_out_vel,racket_des_normal);
	tennis.calc_des_racket_vel(balls_predicted.rows(DX,DZ),balls_out_vel,racket_des_normal,racket_des_vel);

	// place racket centre on the predicted ball
	racket_params.ball_pos = balls_predicted.rows(X,Z);
	racket_params.ball_vel = balls_predicted.rows(DX,DZ);
	racket_params.racket_pos = balls_predicted.rows(X,Z);
	racket_params.racket_vel = racket_des_vel;
	racket_params.racket_normal = racket_des_normal;
	return racket_params;
}

optim_des calc_racket_strategy(const mat & balls_predicted,
		                       const vec2 & ball_land_des,
		                       const double time_land_des,
							   optim_des & racket_params) {

	//static wall_clock timer;
	//timer.tic();
	bool hack = true; // for modifying outgoing ball velocities
	TableTennis tennis = TableTennis(false,false);

	int N = balls_predicted.n_cols;
	mat balls_out_vel = zeros<mat>(3,N);
	mat racket_des_normal = zeros<mat>(3,N);
	mat racket_des_vel = zeros<mat>(3,N);
	tennis.calc_des_ball_out_vel(ball_land_des,time_land_des,hack,balls_predicted,balls_out_vel);
	tennis.calc_des_racket_normal(balls_predicted.rows(DX,DZ),balls_out_vel,racket_des_normal);
	tennis.calc_des_racket_vel(balls_predicted.rows(DX,DZ),balls_out_vel,racket_des_normal,racket_des_vel);

	// place racket centre on the predicted ball
	racket_params.ball_pos = balls_predicted.rows(X,Z);
	racket_params.ball_vel = balls_predicted.rows(DX,DZ);
	racket_params.racket_pos = balls_predicted.rows(X,Z);
	racket_params.racket_vel = racket_des_vel;
	racket_params.racket_normal = racket_des_normal;

	//cout << "Pred. racket time: " << 1000 * timer.toc() << " ms." << endl;
	return racket_params;
}

static void optim_spin_outgoing_ball_vel(const des_ball_data & data,
                                         const bool verbose,
                                         vec3 & est) {

	static double x[3];  /* some initial guess */
	static double minf; /* the minimum objective value, upon return */
	static double init_time;
	static int res; // error code
	static nlopt_opt opt = nlopt_create(NLOPT_LD_MMA, 3);
	nlopt_set_min_objective(opt, calc_landing_res, (void*)&data);
	nlopt_set_xtol_rel(opt, 1e-2);

	for(int i = 0; i < 3; i++) {
		x[i] = est(i);
	}

	init_time = get_time();
	if ((res = nlopt_optimize(opt, x, &minf)) < 0) {
		if (verbose)
			printf("NLOPT failed!\n");
	}
	else {
		if (verbose) {
			printf("NLOPT took %f ms\n", (get_time() - init_time)/1e3);
			printf("Found minimum at f = %0.10g\n", minf);
		}
		for(int i = 0; i < 3; i++) {
			est(i) = x[i];
		}
	}
	//nlopt_destroy(opt);
}

static double calc_landing_res(unsigned n,
                                const double *x,
                                double *grad,
                                void *data) {

	static double dt = 0.02;
    static TableTennis tt = TableTennis(true,false,false); // no contact checking!
    static vec3 vel_out;
    static vec3 out_pos;

    des_ball_data *mydata = (des_ball_data*) data;
    tt.set_topspin(mydata->topspin);
    vel_out(X) = x[0];
    vel_out(Y) = x[1];
    vel_out(Z) = x[2];
    tt.set_ball_state(join_vert(mydata->ball_incoming,vel_out));
	for (int i = 0; i < mydata->time_land_des/dt; i++)
		tt.integrate_ball_state(dt);

	out_pos = tt.get_ball_position();

    if (grad) {

    	grad[0] = mydata->time_land_des * (out_pos(X) - mydata->ball_land_des(X));
    	grad[1] = mydata->time_land_des * (out_pos(Y) - mydata->ball_land_des(Y));
    	grad[2] = mydata->time_land_des * (out_pos(Z) - mydata->ball_land_des(Z));
    	// Finite difference
		/*static double h = 1e-6;
		static double val_plus, val_minus;
		static double xx[3];
		for (unsigned i = 0; i < n; i++)
			xx[i] = x[i];
		for (unsigned i = 0; i < n; i++) {
			xx[i] += h;
			val_plus = calc_landing_res(n, xx, NULL, data);
			xx[i] -= 2*h;
			val_minus = calc_landing_res(n, xx, NULL, data);
			grad[i] = (val_plus - val_minus) / (2*h);
			xx[i] += h;
		}*/
    }

	return pow(norm(out_pos - mydata->ball_land_des),2);
}

}
