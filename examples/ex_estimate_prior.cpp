/*
 * ex_estimate_prior.cpp
 *
 * Estimate initial ball state using
 * nonlinear least squares with NLOPT
 *
 *  Created on: Jul 1, 2017
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

using namespace arma;

#define DIM 7
#define SAMPLES 12
#define DT 0.016

vec times;
mat obs;

void nlopt_optimize(vec6 & est);
void generate_ball_samples(const double topspin,
		vec & times, mat & obs, vec6 & act_state);
void estimate_ball_linear(const int num_samples,
		                  const mat & obs,
						  const vec & t,
						  vec6 & x);
double costfunc(unsigned n, const double *x,
		        double *grad, void *my_func_data);
long get_time();
void actual_ball_example();
void synthetic_ball_example(const double topspin);

int main() {

	synthetic_ball_example(-50);
	//actual_ball_example();
}

void actual_ball_example() {
	// TRIAL 3 BALL ESTIMATION
	times << 0.0 << 0.0140 << 0.0320 << 0.1120 << 0.1280 << endr;
	obs << 0.1379 << 0.1444 << 0.1529 << 0.1752 << 0.1760 << endr
		<< -3.7870 << -3.6840 << -3.6240 << -3.2770 << -3.2270 << endr
		<< -0.4147 << -0.3595 << -0.3327 << -0.2147 << -0.2052 << endr;
	vec6 init_state_est;
	estimate_ball_linear(SAMPLES,obs,times,init_state_est);
	cout << "LS estimate:" << init_state_est.t();
	nlopt_optimize(init_state_est);
	cout << "NLS estimate:" << init_state_est.t();
}

void synthetic_ball_example(const double topspin) {
	times = zeros<vec>(SAMPLES);
	obs = zeros<mat>(3,SAMPLES);
	vec6 init_state_est;
	vec6 actual_state;
	generate_ball_samples(topspin,times,obs,actual_state);
	estimate_ball_linear(SAMPLES,obs,times,init_state_est);
	cout << "Actual state:" << actual_state.t();
	cout << "LS estimate:" << init_state_est.t();
	cout << "Error norm: " << norm(actual_state - init_state_est) << endl;
	nlopt_optimize(init_state_est);
	cout << "NLS estimate:" << init_state_est.t();
	cout << "Error norm: " << norm(actual_state - init_state_est) << endl;
}

void nlopt_optimize(vec6 & est) {

	double x[DIM];  /* some initial guess */
	double minf; /* the minimum objective value, upon return */
	double init_time;
	int res; // error code
	nlopt_opt opt;
	opt = nlopt_create(NLOPT_LD_TNEWTON, DIM);
	nlopt_set_min_objective(opt, costfunc, nullptr);
	nlopt_set_xtol_rel(opt, 1e-2);

	for(int i = 0; i < DIM-1; i++) {
		x[i] = est(i);
	}
	x[DIM-1] = -0.0;

	init_time = get_time();
	if ((res = nlopt_optimize(opt, x, &minf)) < 0) {
	    printf("NLOPT failed!\n");
	}
	else {
		printf("NLOPT took %f ms\n", (get_time() - init_time)/1e3);
	    printf("Found minimum at f = %0.10g\n", minf);
		for(int i = 0; i < DIM-1; i++) {
			est(i) = x[i];
		}
		printf("Topspin est: %f\n", 100*x[DIM-1]);
	}
	nlopt_destroy(opt);
}

/*
 * Generate ball samples using a ball model and random ballgun loc.
 */
void generate_ball_samples(const double topspin,
		        vec & times, mat & obs, vec6 & act_state) {

	arma_rng::set_seed_random();
	TableTennis tt = TableTennis(true,false);
	tt.set_topspin(topspin);
	int ball_launch_side = (randi(1,distr_param(0,2)).at(0));
	tt.set_ball_gun(0.05,ball_launch_side);
	act_state = tt.get_ball_state();
	for (int i = 0; i < SAMPLES; i++) {
		tt.integrate_ball_state(DT);
		times(i) = (i+1)*DT;
		obs.col(i) = tt.get_ball_position();
	}
}

/*
 * Estimate ball state using Least Squares
 */
void estimate_ball_linear(const int num_samples,
		                  const mat & obs,
						  const vec & t,
						  vec6 & x) {

	mat M = zeros<mat>(num_samples,3);
	// and create the data matrix
	for (int i = 0; i < num_samples; i++) {
		M(i,0) = 1.0;
		M(i,1) = t(i);
		M(i,2) = t(i) * t(i);
	}
	// solving for the parameters
	//cout << "Data matrix:" << endl << M << endl;
	mat Beta = solve(M,obs.t());
	//cout << "Parameters:" << endl << Beta << endl;
	x = join_horiz(Beta.row(0),Beta.row(1)).t();
}

/*
 * Cost function
 * Calculates also the gradient if grad is TRUE
 *
 */
double costfunc(unsigned n, const double *x, double *grad, void *data) {

    static TableTennis tt = TableTennis(true,false);
    tt.set_topspin(100*x[DIM-1]);

    if (grad) {
		static double h = 1e-6;
		static double val_plus, val_minus;
		static double xx[DIM];
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

    vec7 state(x);
    vec6 init_state = state.head(DIM-1);
    tt.set_ball_state(init_state);
    tt.integrate_ball_state(times(0));
    double dt;
    double residual = pow(norm(tt.get_ball_position() - obs.col(0)),2);
    for (int i = 1; i < SAMPLES; i++) {
    	dt = times(i) - times(i-1);
    	tt.integrate_ball_state(dt);
    	residual += pow(norm(tt.get_ball_position() - obs.col(i)),2);
    }

    return residual;

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


