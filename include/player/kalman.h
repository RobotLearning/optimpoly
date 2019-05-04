/**
 * @file kalman.h
 *
 * @brief Kalman Filter(s) declarations.
 *
 *  Created on: Jan 25, 2017
 *      Author: okoc
 */

#pragma once
#include <armadillo>
#include <string>

using arma::mat;
using arma::vec;

namespace player {

/**
 * @brief Discrete Kalman filter class with some additional checks.
 *
 * Discrete Kalman Filter class which can discretize continuous model
 * matrices as well as sampling observations.
 *
 */
class KF {

private:

  /**
   * @brief Checks if the model matrices have the right size
   * for initialization
   */
  void check_models(const mat &A, const mat &B, const mat &C) const;

protected:
  mat A_; //!< matrix for the linear(ized) system (discrete)
  mat B_;
  mat C_; //!< observation matrix that EKF can also borrow
  mat Q_; //!< covariance of the process noise (discrete)
  mat R_; //!< covariance of the observation noise (discrete)
  vec x_; //!< state
  mat P_; //!< covariance of the state

  /**
   * @brief Check if the matrix is symmetric positive semidefinite
   *
   * For positive definiteness we look at the value of smallest
   * eigenvalue
   *
   * TODO: surely there must be a faster way!
   *
   */
  void check_spd(const mat &M) const;

  /**
   * Takes cholesky of a positive semi-definite matrix
   * if an eigenvalue is 0 then it leaves 0 be!
   *
   * TODO: maybe move into matrix utilities?
   */
  mat chol_semi(const mat &M) const;

public:
  /**
   * @brief Initialize the Kalman filter with given noise covariance
   * matrices and const observation matrix C.
   *
   * Useful KF initialization in case model is time-varying and/or continuous
   * so needs to be discretized to apply predict/update equations.
   *
   * Sets A and B to effectively to uninitialized
   * (inf in armadillo). Leaves state vector x and covariance matrix P
   * uninitialized by setting them Inf.
   *
   * @param Cin Observation matrix C
   * @param Qin Process (model) covariance Q
   * @param Rin Noise (observation) covariance R
   */
  KF(mat &Cin, mat &Qin, mat &Rin);

  /**
   * @brief Initialize the Kalman filter with given
   * discrete model and noise covariance matrices.
   *
   * Model matrices must have the right size and
   * noise covariances must be symmetric positive definite!
   *
   * Leaves state vector x and covariance matrix P uninitialized
   * by setting them Inf.
   *
   * @param Ain Drift matrix A
   * @param Bin Control matrix B
   * @param Cin Observation matrix C
   * @param Qin Process (model) covariance Q
   * @param Rin Noise (observation) covariance R
   *
   */
  KF(mat &A, mat &B, mat &C, mat &Q, mat &R);

  /**
   * @brief Initialize the Kalman filter with given
   * discrete model and noise matrices
   *
   * Overloaded constructor for KF with additional given
   * initial state mean x0 and variance P0
   *
   * @param x0 Initial state
   * @param P0 Initial state covariance
   * @param Ain Drift matrix A
   * @param Bin Control matrix B
   * @param Cin Observation matrix C
   * @param Qin Process (model) covariance Q
   * @param Rin Noise (observation) covariance R
   *
   */
  KF(vec &x0, mat &P0, mat &A, mat &B, mat &C, mat &Q, mat &R);

  /**
   * @brief Initialize the filter state and the covariance.
   *
   * This is very important function to call, as trying to get
   * filter mean and/or variance before initializing throws an uninitialized
   * error (as const char *)
   *
   * @param x0 initial mean.
   * @param P0 initial covariance.
   *
   */
  void set_prior(const vec &x0, const mat &P0);

  /**
   * @brief Get state mean.
   *
   * @return State mean.
   * @throw Exception if prior was not set before!
   *
   */
  vec get_mean() const;

  /**
   * @brief Get state covariance.
   *
   * @return State covariance.
   * @throw Exception if prior was not set before!
   *
   */
  mat get_covar() const;

  /**
   * @brief Get model matrices based on numeric input
   *
   * @param idx 1 - A, 2 - B, 3 - C, 4 - D respectively.
   * @return model matrix A to C, D is not implemented!
   *
   * @throw Exception if uninitialized (inf in ARMADILLO).
   * D is not implemented!
   *
   */
  mat get_model(int idx) const;

  /**
   * @brief Discretizing continuous model matrices.
   *
   * In case we have a continuous linear(ized) model, we can call
   * this function to first discretize the model.
   *
   * @param Ac continuous (instantaneous) drift model matrix.
   * @param Bc continuous control matrix.
   * @param dt discretize over dt horizon.
   *
   */
  void discretize(const mat &Ac, const mat &Bc, double dt);

  /**
   * @brief Predict next state mean and covariance.
   *
   * Does not use control matrix or control inputs. Useful
   * for table tennis ball.
   *
   */
  void predict();

  /**
   * @brief Predict next state mean and covariance.
   *
   * @param u Control inputs. Must have the same size as columns of B.
   *
   */
  void predict(const vec &u);

  /**
   * @brief Update the mean and variance of the state
   * after making an observation.
   *
   * Simple form without any control input (the usual case).
   *
   * @param y observations. Must have the same size as rows of C.
   *
   */
  void update(const vec &y);

  /**
   * @brief Sample observations up to a horizon size N.
   * @param N Horizon size.
   * @return Sampled future observations.
   */
  mat sample_observations(int N) const;
};

using model = vec (*)(const vec &, const double, const void *p);

/**
 * @brief Extended Kalman Filter.
 *
 * Inherits the usual update method, but uses its own prediction method
 * based on the (nonlinear) function pointer data member.
 */
class EKF : public KF {

private:
  void *fparams_ = nullptr; //!< parameters to function
  double
      outlier_reject_mult_; //!< standard deviation multiplier to reject outliers
  model f_;                 //!< pointer to function for prediction

  /**
   * @brief Linearize the discrete function (that integrates a continuous
   * functions dt seconds) to get Ad matrix
   *
   * Using 'TWO-SIDED-SECANT' to do a stable linearization
   */
  mat linearize(const double dt, const double h) const;

public:
  /**
   * @brief Initialize the Discrete Extended Kalman filter with given noise
   * covariance matrices and const observation matrix C
   *
   * Function pointer is typically a function that integrates an (underlying
   * hidden) continuous model by dt seconds supplied in its second argument
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
  EKF(model fp, mat &Cin, mat &Qin, mat &Rin, double rej_mult = 2.0);

  EKF(const EKF &ekf);            //!< copy constructor
  EKF &operator=(const EKF &ekf); //!< assignment operator

  /** @brief Set function co-parameters for predicting */
  void set_fun_params(void *params) { fparams_ = params; };

  /**
   * @brief Predict dt seconds for mean x and (if flag is true) variance P.
   *
   * @param dt Prediction horizon.
   * @param lin_flag If true, will linearize the nonlinear function
   * around current x and make the covariance update. Useful to turn off for
   * debugging.
   *
   */
  void predict(const double dt, const bool lin_flag = true);

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
  mat predict_path(const double dt, const int N) const;

  /**
   * @brief Checks to see if the ball observation could be an
   * outlier.
   *
   * Using the covariance matrix estimate to detect such an outlier
   * that possibly escaped the elimination from check_blob_validity function
   *
   * The new obs has to be located a certain standard deviations away from last
   * obs
   *
   * @param y Observation to check.
   * @param verbose If flag is TRUE, will print info when outlier is detected
   * @param TRUE if outlier is detected
   *
   */
  bool check_outlier(const vec &obs, const bool verbose = false) const;
};

/**
 * @brief Initialize an EKF
 *
 * Called generally when ball state needs to be reset
 * Useful for passing to Player constructor.
 * @param var_model Process noise multiplier for identity matrix.
 * @param var_noise Observation noise mult. for identity matrix.
 * @param spin Use spin model if true
 * @param out_reject_mult Mult. for outlier rejection.
 * @param topspin Set topspin parameter (NOT state!) for kalman filter
 * prediction if spin is TRUE
 * @return EKF Extended Kalman Filter (state uninitialized!)
 */
EKF init_ball_filter(const double var_model = 0.001,
                     const double var_noise = 0.001, const bool spin = false,
                     const double out_reject_mult = 2.0,
                     const double *topspin = nullptr);

/**
 * @brief Checks to see if the observation is new (updated)
 *
 * The blobs need to be at least tol apart from each other in distance
 *
 */
bool check_new_obs(const arma::vec3 &obs, double tol);

/**
 * @brief Check to see if we want to reset the filter.
 *
 * Basically if a new ball appears 300 ms later than the last new ball
 * we reset the filter.
 *
 */
bool check_reset_filter(const bool newball, const int verbose,
                        const double threshold);

} // namespace player
