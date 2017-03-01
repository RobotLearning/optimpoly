#include <armadillo>
#include "kalman.h"
#include "player.hpp"
#include "tabletennis.h"
#include <cmath>
#include <sys/time.h>
#include "carma.h"

using namespace arma;

/* The data structures */
typedef struct { /*!< joint space state for each DOF */
  double   th;   /*!< theta */
  double   thd;  /*!< theta-dot */
  double   thdd; /*!< theta-dot-dot */
  double   ufb;  /*!< feedback portion of command */
  double   u;    /*!< torque command */
  double   load; /*!< sensed torque */
} SL_Jstate;

typedef struct { /*!< desired values for controller */
  double   th;   /*!< theta */
  double   thd;  /*!< theta-dot */
  double   thdd; /*!< theta-dot-dot */
  double   uff;  /*!< feedforward torque command */
  double   uex;  /*!< externally imposed torque */
} SL_DJstate;

typedef struct { /*!< Cartesian state */
  double   x[NCART+1];    /*!< Position [x,y,z] */
  double   xd[NCART+1];   /*!< Velocity */
  double   xdd[NCART+1];  /*!< Acceleration */
} SL_Cstate;

typedef struct { /*!< Vision Blob */
  char       status;
  SL_Cstate  blob;
} SL_VisionBlob;

/*
 *
 * Interface to the PLAYER class that generates desired hitting trajectories.
 * First initializes the player and then starts calling play() interface function.
 *
 */
void play(const SL_Jstate joint_state[NDOF+1],
		  const SL_VisionBlob blobs[4],
		  SL_DJstate joint_des_state[NDOF+1]) {

	static bool firsttime = true;
	static vec7 q0;
	static vec3 ball_obs;
	static joint qact;
	static joint qdes;
	static Player cp; // centered player

	if (firsttime) {
		for (int i = 0; i < NDOF; i++) {
			q0(i) = joint_state[i+1].th;
		}
		EKF filter = init_filter();
		cp = Player(q0,filter);
		firsttime = false;
	}
	else {
		for (int i = 0; i < NDOF; i++) {
			qact.q(i) = joint_state[i+1].th;
			qact.qd(i) = joint_state[i+1].thd;
			qact.qdd(i) = joint_state[i+1].thdd;
		}
		for (int i = 0; i < NCART; i++)
			ball_obs(i) = blobs[1].blob.x[i+1];
		cp.play(qact,ball_obs,qdes);
	}

	// update desired joint state
	for (int i = 0; i < NDOF; i++) {
		joint_des_state[i].th = qdes.q(i);
		joint_des_state[i].thd = qdes.qd(i);
		joint_des_state[i].thdd = qdes.qdd(i);
	}

}

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
	bool new_ball = check_new_carma_obs(obs);

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
	std::cout << "Data matrix:" << std::endl << M << std::endl;
	mat Beta = solve(M,observations.t());
	std::cout << "Parameters:" << endl << Beta << std::endl;
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
bool check_new_carma_obs(const vec3 & obs) {

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

