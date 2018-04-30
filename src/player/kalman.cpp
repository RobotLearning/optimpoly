/**
 * @file kalman.cpp
 *
 * @brief (Discrete) Kalman filtering class KF in C++
 *
 * Kalman Filter class for basic filtering in C++ using ARMADILLO linear
 * algebra library.
 * Compatible with continuous models by using discretize() method.
 * Can be extended easily (e.g. see EKF).
 *
 * TODO: add square-root form for stability!
 *
 *  Created on: Jan 25, 2017
 *      Author: okoc
 */

#include <iostream>
#include <string>
#include <armadillo>
#include "kalman.h"

using namespace arma;

namespace player {

KF::KF(mat & Cin, mat & Qin, mat & Rin) {

	// checking for correct noise covariances
	check_spd(Qin);
	check_spd(Rin);
	int dimx = Cin.n_cols;

	this->x = datum::inf * ones<vec>(dimx);
	this->P = datum::inf * ones<mat>(dimx,dimx);

	// set values
	//this->A = zeros<mat>(dimx,dimx);
	this->A = datum::inf * ones<mat>(dimx,dimx);
	this->B = datum::inf * ones<mat>(dimx,dimx); // we dont know B's dimensions
	this->C = Cin;
	this->Q = Qin;
	this->R = Rin;
}

KF::KF(mat & Ain, mat & Bin, mat & Cin, mat & Qin, mat & Rin) {


	check_models(Ain,Bin,Cin);
	// checking for correct noise covariances
	check_spd(Qin);
	check_spd(Rin);

	int dimx = Ain.n_rows;
	//int dimu = Bin.n_cols;
	//int dimy = Cin.n_rows;

	// initialize state mean and covariance
	this->x = datum::inf * ones<vec>(dimx);
	this->P = datum::inf * ones<mat>(dimx,dimx);
	//P.eye();

	// set values
	//this->A = zeros<mat>(dimx,dimx);
	this->A = Ain;
	this->B = Bin;
	this->C = Cin;
	this->Q = Qin;
	this->R = Rin;
}

KF::KF(vec & x0, mat & P0, mat & Ain, mat & Bin, mat & Cin, mat & Qin, mat & Rin) {

	check_models(Ain,Bin,Cin);
	// checking for correct noise covariances
	check_spd(Qin);
	check_spd(Rin);
	set_prior(x0,P0);

	// set values
	this->A = Ain;
	this->B = Bin;
	this->C = Cin;
	this->Q = Qin;
	this->R = Rin;
}

void KF::set_prior(const vec & x0, const mat & P0) {

	x = x0;
	P = P0;
}

void KF::check_models(const mat & Ain,
                    const mat & Bin,
                    const mat & Cin) const {

	// checking for correct model matrix sizes
	if (Ain.n_rows != Ain.n_cols)
		cerr << "A must be square!" << endl;
	if (Bin.n_rows != Ain.n_rows)
		cerr << "A and B must have same row size!" << endl;
	if (Cin.n_cols != Ain.n_cols)
		cerr << "C and A must have same column size!" << endl;
	//if (Cin.n_rows != Din.n_rows)
	//	cerr << "C and D must have same row size!" << endl;
}

void KF::check_spd(const mat & M) const {

	if (M.n_cols != M.n_rows)
		cerr << "Covariance matrix must be square!" << endl;
	if (any(vectorise(M - M.t())))
		cerr << "Covariance matrix must be symmetric!" << endl;
	const std::string form = "sa";
	vec eigvals = eig_sym(M);
	if (eigvals(0) < 0.0) {
		throw "Covariance matrix must be positive semidefinite!";
	}
}

mat KF::chol_semi(const mat & M) const {

	int size;
	mat Y;
	try {
		size = M.n_cols;
		check_spd(M);
		Y = chol(M);
	}
	catch (char * err1) {
		std::cout << "Covariance is reset to zero!" << std::endl;
		Y = zeros<mat>(size,size);
	}
	catch (std::runtime_error & err2) {
		vec eigval;
		mat eigvec;
		eig_sym(eigval, eigvec, M);
		Y = eigvec * diagmat(sqrt(eigval));
	}
	return Y;

}

void KF::discretize(const mat & Ac, const mat & Bc, double dt) {

	/*! MATLAB code
	    trick to get discrete time versions
	    Mat = [obj.A, obj.B; zeros(dimu, dimx + dimu)];
	    MD = expm(h * Mat);
	    Ad = MD(1:dimx,1:dimx);
	    Bd = MD(1:dimx,dimx+1:end); */

	int dimu = Bc.n_cols;
	int dimx = x.n_elem;
	mat M = join_vert(join_horiz(Ac, Bc),zeros<mat>(dimu, dimx + dimu));
	mat MD = expmat(dt * M);

	A = MD(span(0, dimx-1), span(0, dimx-1));
	B = MD(span(0, dimx-1), span(dimx, dimx + dimu -1));

}

mat KF::smoothen(const mat & observations) {

	// TODO: Only one iteration?
	throw std::runtime_error("Unimplemented!");

}

vec KF::get_mean() const {

	// quick and dirty check for initialization
	if (!x.is_finite()) {
		throw std::runtime_error("KF not initialized! Please set prior!");
	}

	return x;
}

mat KF::get_covar() const {

	// quick and dirty check for initialization
	if (P(0,0) == datum::inf) {
		throw std::runtime_error("KF not initialized! Please set prior!");
	}

	return P;
}

mat KF::get_model(int idx) const {

	if (idx > 4 || idx < 1) {
		throw std::runtime_error("Index must be 1 to 4 only!");
	}
	mat out;
	switch (idx) {
		case 1: out = A;
				break;
		case 2: out = B;
				break;
		case 3: out = C;
				break;
		case 4:	throw std::runtime_error("D matrix unimplemented!");
				break;
		default: throw std::runtime_error("This shouldn't happen!");
	}
	if (out(0,0) == datum::inf) {
		throw std::runtime_error("Matrix is not initialized!");
	}
	return out;
}

void KF::predict() {

	x = A * x;
	P = A * P * A + Q;
}

void KF::predict(const vec & u) {

	if (u.n_elem != B.n_cols)
		cerr << "Control input u should have the right size!" << endl;

	x = A * x + B * u;
	P = A * P * A + Q;
}

void KF::update(const vec & y) {

	if (y.n_elem != C.n_rows)
		cerr << "Observation y should have the right size!" << endl;

	// innovation sequence
	vec z = y - C * x;
	// innovation covariance
	mat Theta = C * P * C.t() + R;
	// optimal Kalman gain
	mat K = P * C.t() * inv(Theta);

	// update state mean and covariance
	P = P - K * C * P;
	//cout << "obs:" << "\t" << y.t();
	//cout << "x_pre:" << "\t" << x.t();
	x = x + K * z;
	//cout << "x_post:" << "\t" << x.t();
}

mat KF::sample_observations(int N) const {

	// disable messages being printed to the err2 stream
	std::ostream nullstream(0);
	set_stream_err2(nullstream);

	int dimy = C.n_rows;
	int dimx = x.n_elem;
	vec x_next = x;
	mat Y(dimy,N,fill::zeros);
	vec w(dimy,fill::randn);

	Y.col(0) = C * x + chol(R) * w;
	for (int i = 1; i < N; i++) {
		w = randn<vec>(dimy);
		Y.col(i) = C * x_next + chol_semi(R) * w;
		w = randn<vec>(dimx);
		x_next = A * x_next + chol_semi(Q) * w;
	}

	return Y;
}

}
