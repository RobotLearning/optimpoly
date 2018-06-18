/*
 * estimate_ball.cpp
 *
 *  Created on: Jul 9, 2017
 *      Author: okoc
 */

#include <armadillo>
#include "nlopt.h"
#include "optim.h"
#include "tabletennis.h"
#include "kalman.h"

using namespace arma;
static const int fitnorm = 2;

/*
 * Least squares to estimate prior given
 * matrix of observations Y of column length N, each row
 * corresponding to one
 *
 * If observations arrive as column vectors then we take
 * transpose of it.
 *
 */
static void estimate_ball_linear(const mat & observations,
		                          const vec & times,
		                          const bool verbose,
		                          vec6 & init_est);

/*
 * Returns the estimated topspin value
 */
static double nlopt_estimate_ball(const mat & obs,
                                    const vec & times,
                                    const bool verbose,
                                    vec6 & est);

/*
 * Cost function
 * Calculates also the gradient if grad is TRUE
 */
static double calc_residual(unsigned n,
                            const double *x,
                            double *grad,
                            void *data);

/*
 * Return time of day as micro seconds
 */
static long get_time();

namespace optim {

/**
 * @brief Initial ball data time stamps and position observations
 */
struct init_ball_data {
	vec times = zeros<vec>(10);
	mat obs = zeros<mat>(3,10);
};

void estimate_prior(const mat & observations,
        			const mat & times,
					const int & verbose,
					bool & NLOPT_FINISHED,
					player::EKF & filter) {

	NLOPT_FINISHED = false;
	static double topspin = -50;
	vec6 x;
	vec times_z = times - times(0); // times zeroed
	estimate_ball_linear(observations,times_z,verbose > 2,x);
	topspin = nlopt_estimate_ball(observations,times_z,verbose > 2,x);
	mat P;
	P.eye(6,6);
	filter.set_prior(x,P);
	filter.set_fun_params((void*)&topspin);
	filter.update(observations.col(0));

	double dt;
	for (unsigned i = 1; i < times.n_elem; i++) {
		dt = times_z(i) - times_z(i-1);
		filter.predict(dt,true);
		filter.update(observations.col(i));
	}
	//cout << "Detached ball state estimation finished!\n";
	NLOPT_FINISHED = true;
}

}

static long get_time() {
	struct timeval tv;
	if (gettimeofday(&tv, (struct timezone *)0) == 0)
		return (tv.tv_sec*1000*1000 + tv.tv_usec);  //us

	return 0.;
}

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
    //mat Beta = pinv(M,0.01) * observations.t();
    init_est = join_horiz(Beta.row(0),Beta.row(1)).t();

    if (verbose) {
        cout << "Times:" << times.t() << endl;
        cout << "Data:\n" << observations.t() << endl;
        cout << "Initial est:" << init_est.t() << endl;
    }
}

static double nlopt_estimate_ball(const mat & obs,
                                    const vec & times,
                                    const bool verbose,
                                    vec6 & est) {

    using namespace const_tt;
    static const int TS = 6; // topspin
    double lb[7], ub[7];
    double x[7];  /* some initial guess */
    double minf; /* the minimum objective value, upon return */
    double init_time;
    int res; // error code
    nlopt_opt opt;
    optim::init_ball_data data;
    data.obs = obs;
    data.times = times;

    ub[X] = 1.0; ub[Y] = 0.5; ub[Z] = 2.0; ub[DX] = 2.0; ub[DY] = 5.0; ub[DZ] = 4.0; ub[TS] = 1.0;
    lb[X] = -1.0; lb[Y] = -5.0; lb[Z] = -2.0; lb[DX] = -2.0; lb[DY] = -5.0; lb[DZ] = -4.0; lb[TS] = -1.0;

    opt = nlopt_create(NLOPT_LD_TNEWTON, 7);
    nlopt_set_min_objective(opt, calc_residual, (void*)&data);
    nlopt_set_lower_bounds(opt, lb);
    nlopt_set_upper_bounds(opt, ub);
    nlopt_set_xtol_rel(opt, 1e-2);
    //nlopt_set_maxtime(opt, 0.002);
    for(int i = 0; i < 6; i++) {
        x[i] = est(i);
    }
    x[6] = -0.5;

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
            cout << "Initial state est:" << est.t();
            printf("Topspin est: %f\n", 100*x[6]);
        }
    }
    nlopt_destroy(opt);
    return 100*x[6];
}

static double calc_residual(unsigned n,
                            const double *x,
                            double *grad,
                            void *void_data) {

    using namespace player;
    using namespace optim;
    /*const double lambda_spin = 1e-5;
    const double lambda_zvel = 5e-2;
    const double lambda_yvel = 5e-2;*/
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
    double residual = norm(tt.get_ball_position() - data->obs.col(0),fitnorm);
    for (int i = 1; i < num_samples; i++) {
        dt = data->times(i) - data->times(i-1);
        tt.integrate_ball_state(dt);
        residual += norm(tt.get_ball_position() - data->obs.col(i),fitnorm);
    }
    // add regularization to prevent overfitting in spin estimation
    //residual += lambda_spin * pow(x[6] - 0.5,2) + lambda_zvel * pow(x[DZ]-2.5,2) + lambda_yvel * pow(x[DY]-4.5,2);

    return residual;

}
