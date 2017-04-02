/**
 * @file kalman.h
 *
 * @brief Kalman Filter(s) declarations.
 *
 *  Created on: Jan 25, 2017
 *      Author: okoc
 */

#ifndef KALMAN_H_
#define KALMAN_H_

#include <string>

using namespace arma;

/**
 *
 * @brief Kalman filter class with some additional checks.
 */
class KF {

private:

	mat A; // matrices for the linear system
	mat B;

	// checks and balances
	void check_models(const mat & A, const mat & B, const mat & C) const;

protected:

	vec x; // state
	mat P; // covariance of the state
	mat C; // observation matrix that EKF can also borrow

	mat Q; // covariance of the process noise
	mat R; // covariance of the observation noise

	void check_spd(const mat & M) const;
	mat chol_semi(const mat & M) const;

public:

	// constructor useful in continuous/time-varying filtering scenarios
	KF(mat & Cin, mat & Qin, mat & Rin);

	// constructor for matrices and zero state
	KF(mat & A, mat & B, mat & C, mat & Q, mat & R);

	// initialization for matrices and state
	KF(vec & x0, mat & P0, mat & A, mat & B, mat & C, mat & Q, mat & R);

	// set initial mean and covariance
	void set_prior(const vec & x0, const mat & P0);

	// return mean and covariance
	vec get_mean() const;
	mat get_covar() const;

	// return model (used for testing / comparison)
	mat get_model(int idx) const;

	// discretize (time varying) continuous matrices
	void discretize(const mat & Ac, const mat & Bc, double dt);

	// predict update equations
	void predict();
	void predict(const vec & u);
	void update(const vec & y);

	// sample observations
	mat sample_observations(int N) const;

	// Kalman smoother, not implemented yet!
	mat smoothen(const mat & observations);
};

/**
 * @brief Extended Kalman Filter.
 *
 * Inherits the usual update method, but uses its own prediction method
 * based on the (nonlinear) function pointer data member.
 * predict() calls linearize if linearization flag is TRUE.
 * predict_path() does NOT call linearize!
 *
 */
class EKF : public KF {

private:

	// function pointer
	vec (*f)(const vec &,double);

	// linearize function to get Ad matrix
	mat linearize(double dt, double h) const;

public:

	// initializing with a function pointer
	EKF(vec (*fp)(const vec &,double), mat & Cin, mat & Qin, mat & Rin);

	// overriding predict() of KF superclass
	void predict(double dt, bool lin_flag = true);

	// predict future path
	mat predict_path(double dt, int N);

	// check observation for outlier (used after prediction)
	bool check_outlier(const vec & obs, const bool verbose = false) const;
};


#endif /* KALMAN_H_ */
