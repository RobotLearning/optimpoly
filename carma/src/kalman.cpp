/*
 * kalman.cpp
 *
 * (Discrete) Kalman filtering class in C++
 *
 * Compatible with continuous models by using discretize() method.
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

using namespace std;
using namespace arma;

/*
 *
 * Initialize the Kalman filter with given noise covariance
 * matrices and const observation matrix C
 *
 * Sets A and B to effectively to uninitialized
 * (inf in armadillo)
 *
 * Leaves state vector x and covariance matrix P uninitialized
 * by setting them Inf
 *
 * This is useful in scenarios where model is time-varying and/or continuous
 * so needs to be discretized to apply predict/update equations
 *
 */
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

/*
 * Initialize the Kalman filter with given
 * discrete model and noise covariance matrices
 *
 * Model matrices must have the right size and
 * noise covariances must be symmetric positive definite!
 *
 * Leaves state vector x and covariance matrix P uninitialized
 * by setting them Inf
 *
 */
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

/*
 * Initialize the Kalman filter with given
 * discrete model and noise matrices
 *
 * Overloaded constructor with given
 * initial state mean x0 and variance P0
 *
 */
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

KF::~KF() {
	// dont do anything special
}

/*
 * Initialize the filter state and the covariance
 *
 */
void KF::set_prior(const vec & x0, const mat & P0) {

	x = x0;
	P = P0;
}

/*
 * Checks if the model matrices have the right size
 * for initialization
 */
void KF::check_models(const mat & Ain, const mat & Bin, const mat & Cin) const {

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

/*
 * Check if the matrix is symmetric positive semidefinite
 *
 * For positive definiteness we look at the value of smallest
 * eigenvalue
 *
 * TODO: surely there must be a faster way!
 *
 */
void KF::check_spd(const mat & M) const {

	if (M.n_cols != M.n_rows)
		cerr << "Covariance matrix must be square!" << endl;
	if (any(vectorise(M - M.t())))
		cerr << "Covariance matrix must be symmetric!" << endl;
	const string form = "sa";
	vec eigvals = eig_sym(M);
	if (eigvals(0) < 0.0) {
		throw "Covariance matrix must be positive semidefinite!";
	}
}

/*
 * Takes cholesky of a positive semi-definite matrix
 * if an eigenvalue is 0 then it leaves 0 be!
 *
 * TODO: maybe move into matrix utilities?
 */
mat KF::chol_semi(const mat & M) const {

	int size;
	mat Y;
	try {
		size = M.n_cols;
		check_spd(M);
		Y = chol(M);
	}
	catch (char * err1) {
		cout << "Covariance is reset to zero!" << endl;
		Y = zeros<mat>(size,size);
	}
	catch (runtime_error & err2) {
		vec eigval;
		mat eigvec;
		eig_sym(eigval, eigvec, M);
		Y = eigvec * diagmat(sqrt(eigval));
	}
	return Y;

}

/*
 * Discretizing continuous model matrices
 *
 */
void KF::discretize(const mat & Ac, const mat & Bc, double dt) {

	// MATLAB code
	//    % trick to get discrete time versions
	//    Mat = [obj.A, obj.B; zeros(dimu, dimx + dimu)];
	//    MD = expm(h * Mat);
	//    Ad = MD(1:dimx,1:dimx);
	//    Bd = MD(1:dimx,dimx+1:end);

	int dimu = Bc.n_cols;
	int dimx = x.n_elem;
	mat M = join_vert(join_horiz(Ac, Bc),zeros<mat>(dimu, dimx + dimu));
	mat MD = expmat(dt * M);

	A = MD(span(0, dimx-1), span(0, dimx-1));
	B = MD(span(0, dimx-1), span(dimx, dimx + dimu -1));

}

/*
 *
 * Kalman smoothing using a backwards Rauch-recursion
 * Assumes that the KF has been initialized and priors set
 *
 * Returns smoothened observations.
 *
 * TODO: Only one iteration?
 *
 *
 */
mat KF::smoothen(const mat & observations) {

	// TODO:
	throw "Unimplemented!";

}

/*
 * Get state mean
 *
 * Throws error message if prior is not set
 *
 */
vec KF::get_mean() const {

	// quick and dirty check for initialization
	if (x(0) == datum::inf) {
		throw "KF not initialized! Please set prior!";
	}

	return x;
}

/*
 * Get state covariance
 *
 * Throws error message if prior is not set
 *
 */
mat KF::get_covar() const {

	// quick and dirty check for initialization
	if (P(0,0) == datum::inf) {
		throw "KF not initialized! Please set prior!";
	}

	return P;
}

/*
 *
 * Get model matrices based on numeric input
 * 1 - A, 2 - B, 3 - C, 4 - D respectively.
 *
 * Throws error if uninitialized (inf in ARMADILLO).
 * D is not implemented!
 *
 */
mat KF::get_model(int idx) const {

	if (idx > 4 || idx < 1) {
		throw "Index must be 1 to 4 only!";
	}
	mat out;
	switch (idx) {
		case 1: out = A;
				break;
		case 2: out = B;
				break;
		case 3: out = C;
				break;
		case 4:	throw "D matrix unimplemented!";
				break;
		default: throw "This shouldn't happen!";
	}
	if (out(0,0) == datum::inf) {
		throw "Matrix is not initialized!";
	}
	return out;
}

/*
 * Predict next state mean and covariance
 *
 */
void KF::predict() {

	x = A * x;
	P = A * P * A + Q;
}

/*
 * Predict next state mean and covariance
 *
 */
void KF::predict(const vec & u) {

	if (u.n_elem != B.n_cols)
		cerr << "Control input u should have the right size!" << endl;

	x = A * x + B * u;
	P = A * P * A + Q;
}

/*
 * Update the mean and variance of the state
 * after making an observation
 *
 * Simple form without any control input (the usual case)
 *
 */
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
	x = x + K * z;
}

/*
 * Sample observations up to a horizon size N
 */
mat KF::sample_observations(int N) const {

	// disable messages being printed to the err2 stream
	ostream nullstream(0);
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
