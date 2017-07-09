/*
 * estimate_ball.cpp
 *
 *  Created on: Jul 9, 2017
 *      Author: okoc
 */

#include <armadillo>
#include "nlopt.h"
#include "player.hpp"
#include "tabletennis.h"
#include "kalman.h"

using namespace arma;

static void estimate_ball_linear(const mat & observations,
		                  const vec & times,
					      const bool verbose,
					      vec6 & init_est);
static void nlopt_estimate_ball(const mat & obs, const vec & times, const bool verbose, vec6 & est);
static double calc_residual(unsigned n, const double *x, double *grad, void *data);
static long get_time();

struct init_ball_data {
	vec times;
	mat obs;
};

void estimate_prior(const mat & observations,
        			const vec & times,
					const bool verbose,
					EKF & filter) {

	vec6 x;
	vec times_z = times - times(0); // times zeroed
	estimate_ball_linear(observations,times_z,verbose,x);
	nlopt_estimate_ball(observations,times_z,verbose,x);
	mat P;
	P.eye(6,6);
	filter.set_prior(x,P);
	filter.update(observations.col(0));

	double dt;
	for (unsigned i = 1; i < times.n_elem; i++) {
		dt = times_z(i) - times_z(i-1);
		filter.predict(dt,true);
		filter.update(observations.col(i));
	}
}

/*
 * Least squares to estimate prior given
 * matrix of observations Y of column length N, each row
 * corresponding to one
 *
 * If observations arrive as column vectors then we take
 * transpose of it.
 *
 * Velocity estimation is biased, we multiply velocities by 0.8
 * since LSE often overestimates the model with spin.
 *
 *
 */
static void estimate_ball_linear(const mat & observations,
		                  const vec & times,
					      const bool verbose,
					      vec6 & init_est) {

	int num_samples = times.n_elem;
	mat M = zeros<mat>(num_samples,3);

	// and create the data matrix
	for (int i = 0; i < num_samples; i++) {
		M(i,0) = 1.0;
		M(i,1) = times(i);
		M(i,2) = times(i) * times(i);
	}
	// solving for the parameters
	mat Beta = solve(M,observations.t());
	init_est = join_horiz(Beta.row(0),Beta.row(1)).t();

	if (verbose) {
		cout << "Times:" << times.t() << endl;
		cout << "Data:\n" << observations.t() << endl;
		cout << "Initial est:" << init_est.t() << endl;
	}
}

static void nlopt_estimate_ball(const mat & obs, const vec & times, const bool verbose, vec6 & est) {

	double x[7];  /* some initial guess */
	double minf; /* the minimum objective value, upon return */
	double init_time;
	int res; // error code
	nlopt_opt opt;
	init_ball_data data;
	data.obs = obs;
	data.times = times;

	opt = nlopt_create(NLOPT_LD_TNEWTON, 7);
	nlopt_set_min_objective(opt, calc_residual, (void*)&data);
	nlopt_set_xtol_rel(opt, 1e-2);

	for(int i = 0; i < 6; i++) {
		x[i] = est(i);
	}
	x[6] = -0.0;

	init_time = get_time();
	if ((res = nlopt_optimize(opt, x, &minf)) < 0) {
	    printf("NLOPT failed in ball state estimation!\n");
	}
	else {
		for(int i = 0; i < 6; i++) {
			est(i) = x[i];
		}
		if (verbose) {
			printf("NLOPT took %f ms\n", (get_time() - init_time)/1e3);
			printf("Found minimum at f = %0.10g\n", minf);
			cout << "Initial state est:" << est.t() << endl;
			printf("Topspin est: %f\n", 100*x[6]);
		}
	}
	nlopt_destroy(opt);
}

/*
 * Cost function
 * Calculates also the gradient if grad is TRUE
 *
 */
static double calc_residual(unsigned n, const double *x, double *grad, void *void_data) {

    static TableTennis tt = TableTennis(true,false);
    tt.set_topspin(100*x[6]);
    init_ball_data *data = (init_ball_data*) void_data;
    int num_samples = data->times.n_elem;

    if (grad) {
		static double h = 1e-6;
		static double val_plus, val_minus;
		static double xx[7];
		for (unsigned i = 0; i < n; i++)
			xx[i] = x[i];
		for (unsigned i = 0; i < n; i++) {
			xx[i] += h;
			val_plus = calc_residual(n, xx, NULL, data);
			xx[i] -= 2*h;
			val_minus = calc_residual(n, xx, NULL, data);
			grad[i] = (val_plus - val_minus) / (2*h);
			xx[i] += h;
		}
    }

    vec7 state(x);
    vec6 init_state = state.head(6);
    tt.set_ball_state(init_state);
    tt.integrate_ball_state(data->times(0));
    double dt;
    double residual = pow(norm(tt.get_ball_position() - data->obs.col(0)),2);
    for (int i = 1; i < num_samples; i++) {
    	dt = data->times(i) - data->times(i-1);
    	tt.integrate_ball_state(dt);
    	residual += pow(norm(tt.get_ball_position() - data->obs.col(i)),2);
    }

    return residual;

}

/*
 * Checks to see if the observation is new (updated)
 *
 * The blobs need to be at least tol apart from each other in distance
 *
 */
bool check_new_obs(const vec3 & obs, double tol) {

	static vec3 last_obs = zeros<vec>(3);

	if (norm(obs - last_obs) > tol) {
		last_obs = obs;
		return true;
	}
	return false;
}

/*
 * Check to see if we want to reset the filter.
 *
 * Basically if a new ball appears 300 ms later than the last new ball
 * we reset the filter.
 *
 */
bool check_reset_filter(const bool newball, const int verbose, const double threshold) {

	bool reset = false;
	static int reset_cnt = 0;
	static bool firsttime = true;
	static wall_clock timer;

	if (firsttime) {
		firsttime = false;
		timer.tic();
	}

	if (newball) {
		if (timer.toc() > threshold) {
			reset = true;
			if (verbose > 0) {
				std::cout << "Resetting filter! Count: " << ++reset_cnt << std::endl;
			}
		}
		timer.tic();
	}
	return reset;
}

/*
 * Initialize an Extended Kalman Filter
 * useful for passing to Player constructor
 *
 */
EKF init_filter(double std_model, double std_noise, bool spin) {

	mat C = eye<mat>(3,6);
	mat66 Q = std_model * eye<mat>(6,6);
	mat33 R = std_noise * eye<mat>(3,3);
	if (spin)
		return EKF(calc_spin_ball,C,Q,R);
	else
		return EKF(calc_next_ball,C,Q,R);
}

/*
 * Return time of day as micro seconds
 */
static long get_time() {
	struct timeval tv;
	if (gettimeofday(&tv, (struct timezone *)0) == 0)
		return (tv.tv_sec*1000*1000 + tv.tv_usec);  //us

	return 0.;
}
