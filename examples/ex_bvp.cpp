/*
 * ex_bvp.cpp
 *
 * Example Boundary Value Problem with spin
 * where we investigate the solution speed (time elapsed) using optimization
 * (TNEWTON in NLOPT)
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

double topspin;
double time_land_des;
vec2 ball_land_des;
vec6 ballin;

void nlopt_optimize(vec3 & est);
double costfunc(unsigned n, const double *x, double *grad, void *data);
long get_time();

int main() {

	// initialize ball with certain position, incoming velocity and spin
	// ignore spin and calculate outgoing velocity with ballistic model
	// (no airdrag)

	// initialize optimization with ballistic solution
	// solve optimization and output the time elapsed

	arma_rng::set_seed_random();
	topspin = -50.0;
	TableTennis tt = TableTennis(true,false);
	tt.set_topspin(topspin);
	int ball_launch_side = (randi(1,distr_param(0,2)).at(0));
	tt.set_ball_gun(0.05,ball_launch_side);
	for (int i = 0; i < 50; i++)
		tt.integrate_ball_state(0.02);

	ball_land_des = {0.0, dist_to_table - 3*table_length/4.0};
	time_land_des = 0.6;
	vec3 ball_vel_out;
	ballin = tt.get_ball_state();
	tt.calc_des_ball_out_vel(ball_land_des,time_land_des,ballin,ball_vel_out);
	cout << "Incoming ball: " << ballin.t();
	cout << "Estimated outgoing ball velocity: " << ball_vel_out.t();
	nlopt_optimize(ball_vel_out);
	cout << "Solution to BVP: " << ball_vel_out.t();
}

void nlopt_optimize(vec3 & est) {

	double x[3];  /* some initial guess */
	double minf; /* the minimum objective value, upon return */
	double init_time;
	int res; // error code
	nlopt_opt opt;
	opt = nlopt_create(NLOPT_LD_MMA, 3);
	nlopt_set_min_objective(opt, costfunc, nullptr);
	nlopt_set_xtol_rel(opt, 1e-2);

	for(int i = 0; i < 3; i++) {
		x[i] = est(i);
	}
	x[0] = 1.1 * x[0];
	x[1] = 1.2 * x[1];
	x[2] = 0.9 * x[2];

	init_time = get_time();
	if ((res = nlopt_optimize(opt, x, &minf)) < 0) {
	    printf("NLOPT failed!\n");
	}
	else {
		printf("NLOPT took %f ms\n", (get_time() - init_time)/1e3);
	    printf("Found minimum at f = %0.10g\n", minf);
		for(int i = 0; i < 3; i++) {
			est(i) = x[i];
		}
	}
	nlopt_destroy(opt);
}

/*
 * Cost function
 * Calculates also the gradient if grad is TRUE
 *
 */
double costfunc(unsigned n, const double *x, double *grad, void *data) {

	static double dt = 0.02;
    static TableTennis tt = TableTennis(true,false);
    static vec3 vel_out;
    static vec6 init_state;
    static vec3 out_pos;
    static vec3 out_des_pos = {ball_land_des(X), ball_land_des(Y), floor_level - table_height + ball_radius};
    static int N = time_land_des/dt;
    tt.set_topspin(topspin);

    if (grad) {
		static double h = 1e-6;
		static double val_plus, val_minus;
		static double xx[3];
		for (unsigned i = 0; i < n; i++)
			xx[i] = x[i];
		for (unsigned i = 0; i < n; i++) {
			xx[i] += h;
			val_plus = costfunc(n, xx, NULL, data);
			xx[i] -= 2*h;
			val_minus = costfunc(n, xx, NULL, data);
			grad[i] = (val_plus - val_minus) / (2*h);
			xx[i] += h;
		}
    }

    vel_out(X) = x[0];
    vel_out(Y) = x[1];
    vel_out(Z) = x[2];
    init_state = join_vert(ballin.head(3),vel_out);
    tt.set_ball_state(init_state);
	for (int i = 0; i < N; i++)
		tt.integrate_ball_state(dt);

	out_pos = tt.get_ball_position();
	return pow(norm(out_pos - out_des_pos),2);

}

/*
 * Return time of day as micro seconds
 */
long get_time() {
	struct timeval tv;
	if (gettimeofday(&tv, (struct timezone *)0) == 0)
		return (tv.tv_sec*1000*1000 + tv.tv_usec);  //us

	return 0.;
}
