/**
 * @file extkalman.cpp
 *
 * @brief (Discrete) Extended Kalman filtering class EKF in C++
 *
 * Kalman Filter class for basic filtering in C++ using ARMADILLO linear
 * algebra library.
 * Compatible with continuous models by using discretize() method.
 *
 *  Created on: Feb 2, 2017
 *      Author: okoc
 */

#include <armadillo>
#include "kalman.h"

using namespace arma;

/**
 * @brief Initialize the Discrete Extended Kalman filter with given noise covariance
 * matrices and const observation matrix C
 *
 * Function pointer is typically a function that integrates an (underlying hidden)
 * continuous model by dt seconds supplied in its second argument
 *
 * Sets D matrix to zero, sets A and B to effectively to uninitialized
 * (inf in armadillo)
 *
 * This is useful in scenarios where model is time-varying and/or continuous
 * so needs to be discretized to apply predict/update equations
 *
 * TODO: not added pointers to functions with also control input dependence!
 *
 * @param fp (Nonlinear) function pointer.
 * @param Cin Observation matrix C.
 * @param Qin Process noise covariance Q.
 * @param Rin Observation noise covariance R.
 * @param out_rej_mult Outlier rejection standard deviation multiplier
 */
EKF::EKF(vec (*fp)(const vec &, const double, const void *p), mat & Cin, mat & Qin, mat & Rin, double out_rej_mult) : KF(Cin,Qin,Rin) {

	this->f = fp;
	outlier_reject_mult = out_rej_mult;
	/*if (x(0) == datum::inf) {
		std::cout
		<< "EKF not initialized! Call set_prior before filtering!"
		<< std::endl;
	}*/
}

/*
 * Linearize the discrete function (that integrates a continuous functions dt seconds)
 * to get Ad matrix
 *
 * Using 'TWO-SIDED-SECANT' to do a stable linearization
 *
 */
mat EKF::linearize(const double dt, const double h) const {

	int dimx = x.n_elem;
	static mat delta = h * eye<mat>(dimx,dimx);
	static mat dfdx = zeros<mat>(dimx,dimx);
	static mat fPlus = zeros<mat>(dimx,dimx);
	static mat fMinus = zeros<mat>(dimx,dimx);
	for (int i = 0; i < dimx; i++) {
		fPlus.col(i) = this->f(x + delta.col(i),dt, fparams);
		fMinus.col(i) = this->f(x - delta.col(i),dt, fparams);
	}
	return (fPlus - fMinus) / (2*h);
}

/**
 * @brief Predict dt seconds for mean x and (if flag is true) variance P.
 *
 * @param dt Prediction horizon.
 * @param lin_flag If true, will linearize the nonlinear function
 * around current x and make the covariance update. Useful to turn off for debugging.
 *
 */
void EKF::predict(const double dt, const bool lin_flag) {

	x = this->f(x,dt,fparams);
	if (lin_flag) {
		//cout << "A = \n" << linearize(dt,0.01);
		mat A = linearize(dt,0.0001);
		P = A * P * A.t() + Q;
		//cout << "P = \n" << P << "A = \n" << A;
	}
}

/**
 * @brief Predict future path of the estimated object
 * up to a horizon size N
 *
 * Does not update covariances!
 * @param dt Prediction step
 * @param N Number of (small) prediction steps.
 * @return Matrix of ball mean and variances as columns.
 *
 */
mat EKF::predict_path(const double dt, const int N) {

	int dimx = x.n_elem;
	mat X(dimx,N);

	// save mean and covariance
	vec x0 = x;
	mat P0 = P;

	for (int i = 0; i < N; i++) {
		predict(dt,false);
		X.col(i) = x;
	}
	// load mean and variance
	this->set_prior(x0,P0);
	return X;
}

/**
 * @brief Checks to see if the ball observation could be an
 * outlier.
 *
 * Using the covariance matrix estimate to detect such an outlier
 * that possibly escaped the elimination from check_blob_validity function
 *
 * The new obs has to be located a certain standard deviations away from last obs
 *
 * @param y Observation to check.
 * @param verbose If flag is TRUE, will print info when outlier is detected
 * @param TRUE if outlier is detected
 *
 */
bool EKF::check_outlier(const vec & y, const bool verbose) const {

	static int dim_y = y.n_elem;
	bool outlier = true;
	vec threshold = outlier_reject_mult * arma::sqrt(P.diag());
	vec inno = clamp(y - C * x, -10, 10);
	//std::cout << "INNO:" << inno.t() << "THRESH:" << threshold.t();

	if (all(abs(inno) < threshold.head(dim_y))) {
		outlier = false;
	}
	else {
		if (verbose) {
			std::cout << "Outlier:" << y.t()
					  //<< "State: \t" << x.t()
			          << "Inno:  " << inno.t()
					  << "Thresh: " << threshold.head(dim_y).t();
		}
	}
	return outlier;
}
