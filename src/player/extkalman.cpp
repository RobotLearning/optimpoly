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

EKF::EKF(vec (*fp)(const vec &, const double, const void *p),
        mat & Cin,
        mat & Qin,
        mat & Rin,
        double out_rej_mult) : KF(Cin,Qin,Rin) {

	this->f = fp;
	outlier_reject_mult = out_rej_mult;
	/*if (x(0) == datum::inf) {
		std::cout
		<< "EKF not initialized! Call set_prior before filtering!"
		<< std::endl;
	}*/
}

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

void EKF::predict(const double dt, const bool lin_flag) {

	x = this->f(x,dt,fparams);
	if (lin_flag) {
		//cout << "A = \n" << linearize(dt,0.01);
		mat A = linearize(dt,0.0001);
		P = A * P * A.t() + Q;
		//cout << "P = \n" << P << "A = \n" << A;
	}
}

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
