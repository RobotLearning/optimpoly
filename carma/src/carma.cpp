#include <armadillo>
#include "kalman.h"
#include "tabletennis.h"
#include <cmath>
#include <sys/time.h>
#include "carma.h"

using namespace std;
using namespace arma;

/*
 * Filter the blob information with a kalman Filter and save the result in filtState.
 * KF is only used in simulation mode.
 *
 * Returns the bounce variable.
 * It will be reset to FALSE when filter is re-initialized.
 * Bounce variable is used for legal ball detection
 *
 */
void ekf(double x[], double y[], double racket_pos[], double racket_orient[], int *reset) {

	static wall_clock timer;
	static bool firsttime = true;
	static const int min_obs = 5;
	static int num_obs = 0;
	// observation matrix
	static mat OBS = zeros<mat>(3,min_obs);
	static vec TIMES = zeros<vec>(min_obs);
	static vec3 obs;
	static double t_cum;
	static EKF* filter;

	if (firsttime) {
		firsttime = false;
		timer.tic();
	}

	// get elapsed time
	double dt = timer.toc();
	t_cum += dt;

	if (*reset) {
		num_obs = 0;
		*reset = false;
		t_cum = 0.0; // t_cumulative
	}

	obs << y[0] << endr << y[1] << endr << y[2] << endr;
	bool new_ball = check_new_obs(obs);

	if (num_obs < min_obs) {
		if (new_ball) {
			TIMES(num_obs) = t_cum;
			OBS.col(num_obs) = obs;
			num_obs++;
			x[0] = obs(0); x[1] = obs(1); x[2] = obs(2);
			if (num_obs == min_obs) {
				filter = init_ball_filter(OBS,TIMES);
				pass_mean_estimate(x,filter);
			}
		}
	}
	else { // comes here if there are enough balls to start filter
		filter->predict(dt);
		if (new_ball)
			filter->update(obs);
		pass_mean_estimate(x,filter);
	}
	timer.tic();
}

/*
 * Pass mean estimate in case filter is initialized
 */
void pass_mean_estimate(double x[], EKF * filter) {

	vec6 est = filter->get_mean();
	for (int i = 0; i < 6; i++) {
		x[i] = est(i);
	}
}

/*
 * Initialize filter (only once!)
 */
EKF* init_ball_filter(const mat & obs, const vec & times) {

	double std = 0.01;
	mat C = eye<mat>(3,6);
	mat66 Q = zeros<mat>(6,6);
	mat33 R = std * eye<mat>(3,3);

	vec6 x; mat66 P;
	estimate_prior(obs,times,x,P);
	EKF* filter = new EKF(calc_next_ball,C,Q,R);
	filter->set_prior(x,P);
	filter->update(obs.col(0));

	double dt;
	for (int i = 1; i < times.n_elem; i++) {
		dt = times(i) - times(i-1);
		filter->predict(dt);
		filter->update(obs.col(i));
	}

	return filter;
}

/*
 * Empirical Bayes procedure to estimate prior given
 * matrix of observations Y of column length N, each row
 * corresponding to one
 *
 * If observations arrive as column vectors then we take
 * transpose of it.
 *
 *
 */
void estimate_prior(const mat & observations, const vec & times,
		            vec & x, mat & P) {

	int num_samples = times.n_elem;
	mat M = zeros<mat>(num_samples,3);

	// and create the data matrix
	for (int i = 0; i < num_samples; i++) {
		M(i,0) = 1.0;
		M(i,1) = times(i);
		M(i,2) = times(i) * times(i);
	}
	// solving for the parameters
	cout << "Data matrix:" << endl << M << endl;
	mat Beta = solve(M,observations.t());
	cout << "Parameters:" << endl << Beta << endl;
	x = join_horiz(Beta.row(0),Beta.row(1)).t(); //vectorise(Beta.rows(0,1));
	P.eye(6,6);
	P *= 0.1;
}

/*
 * Checks to see if the observation is new (updated)
 *
 * The blobs need to be at least tol = 1e-3 apart from each other in distance
 *
 */
bool check_new_obs(const vec3 & obs) {

	static vec3 last_obs = zeros<vec>(3);
	static const double tol = 1e-3;

	if (norm(obs - last_obs) > tol) {
		last_obs = obs;
		return true;
	}
	return false;
}

/*
 * Inverts the given SL matrix matc [with indices starting from 1]
 * by initializing ARMADILLO equivalent and inverting it
 *
 * After inversion the matrix contents are overwritten to the SL input matrix (double**)
 * Memory for output matrix needs to be allocated before in SL (e.g. my_matrix()).
 *
 */
void invert_matrix(double** matc, int nrows, double** out) {

	double* bowels = new double[nrows*nrows];

	// fill the bowels
	int i,j;
	for(i = 1; i <= nrows; i++) {
		for(j = 1; j <= nrows; j++) {
			bowels[(j-1)*nrows + (i-1)] = matc[i][j]; // mat from SL
		}
	}

	// initialize ARMA mat
	mat A(bowels,nrows,nrows);
	mat B = inv(A);

	for(i = 1; i <= nrows; i++) {
		for(j = 1; j <= nrows; j++) {
			out[i][j] = B(i-1,j-1); // mat from ARMADILLO
		}
	}

}

/*
 * Taking pseudoinverse with default tolerance using ARMADILLO
 *
 * Output matrix must be the transpose of the input matrix (so make sure
 * you initialize properly!).
 *
 *
 */
void pinv_matrix(double** matc, int nrows, int ncols, double** out) {

	double* bowels = new double[nrows*ncols];

	// fill the bowels
	int i,j;
	for(i = 1; i <= nrows; i++) {
		for(j = 1; j <= ncols; j++) {
			bowels[(j-1)*nrows + (i-1)] = matc[i][j]; // mat from SL
		}
	}

	// initialize ARMA mat
	mat A(bowels,nrows,ncols);
	mat B = pinv(A);

	for(i = 1; i <= ncols; i++) {
		for(j = 1; j <= nrows; j++) {
			out[i][j] = B(i-1,j-1); // mat from ARMADILLO
		}
	}
}

