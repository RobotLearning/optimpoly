/*
 * extkalman.cpp
 *
 *  Created on: Feb 2, 2017
 *      Author: okoc
 */

#include <armadillo>
#include "kalman.h"

using namespace std;
using namespace arma;

/*
 *
 * Initialize the Discrete Extended Kalman filter with given noise covariance
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
 */
EKF::EKF(vec (*fp)(const vec &, double), mat & Cin, mat & Qin, mat & Rin) : KF(Cin,Qin,Rin) {

	this->f = fp;

	if (x(0) == datum::inf) {
		cout << "EKF not initialized! Call set_prior before filtering!" << endl;
	}
}

/*
 * Linearize the discrete function (that integrates a continuous functions dt seconds)
 * to get Ad matrix
 *
 * Using 'TWO-SIDED-SECANT' to do a stable linearization
 *
 */
mat EKF::linearize(double dt) const {

	int dimx = x.n_elem;
	mat delta = dt * eye<mat>(dimx,dimx);
	mat dfdx = zeros<mat>(dimx,dimx);
	mat fPlus = zeros<mat>(dimx,dimx);
	mat fMinus = zeros<mat>(dimx,dimx);
	for (int i = 0; i < dimx; i++) {
		fPlus.col(i) = this->f(x + delta.col(i),dt);
		fMinus.col(i) = this->f(x - delta.col(i),dt);
	}
	return (fPlus - fMinus) / (2*dt);
}

/*
 * Predict dt seconds for mean x and variance P.
 * Linearizes the nonlinear function around current x
 * to make the covariance update
 *
 */
void EKF::predict(double dt) {

	mat A = linearize(dt);
	x = this->f(x,dt);
	P = A * P * A.t() + Q;
}

/*
 * Predict future path of the estimated object
 * up to a horizon size N
 *
 */
mat EKF::predict_path(double dt, int N) {

	int dimx = x.n_elem;
	mat X(dimx,N,fill::zeros);

	// save mean and covariance
	vec x0 = x;
	mat P0 = P;

	for (int i = 0; i < N; i++) {
		predict(dt);
		X.col(i) = x;
	}
	// load mean and variance
	this->set_prior(x0,P0);
	return X;
}


