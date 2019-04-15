/*
 * utils.c
 *
 *  Created on: Feb 7, 2017
 *      Author: okoc
 */

#include "utils.h"
#include "constants.h"
#include "json.hpp"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "sys/time.h"
#include <armadillo>

using namespace const_tt;
using json = nlohmann::json;

arma::mat json2mat(const json &jobj) {

  int n = jobj.size();
  int m = jobj[0].size();
  arma::mat arma_mat(n, m, arma::fill::zeros);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      arma_mat(i, j) = jobj[i][j];
    }
  }
  return arma_mat;
}

arma::vec json2vec(const json &jobj) {

  int n = jobj.size();
  arma::vec arma_vec(n, arma::fill::zeros);

  for (int i = 0; i < n; i++) {
    arma_vec(i) = jobj[i];
  }
  return arma_vec;
}

/*
 * Prints the 2*DOF + 1 dimensional solution in user-friendly format
 */
void print_optim_vec(const double *x) {

  int i;
  printf("qf = [");
  for (i = 0; i < NDOF; i++) {
    printf("%.2f  ", x[i]);
  }
  printf("]\n");
  printf("qfdot = [");
  for (i = 0; i < NDOF; i++) {
    printf("%.2f  ", x[i + NDOF]);
  }
  printf("]\n");
  printf("T = %.2f\n", x[2 * NDOF]);
}

/*
 * Return the maximum absolute value of an array with length size
 */
double max_abs_array(const double *x, const int length) {

  double val = x[0];
  for (int i = 1; i < length; i++) {
    if (fabs(x[i]) > val) {
      val = fabs(x[i]);
    }
  }
  return val;
}

/*
 * Return the maximum value of an array with length size
 */
double max_array(const double *x, const int length) {

  double val = x[0];
  for (int i = 1; i < length; i++) {
    if (x[i] > val) {
      val = x[i];
    }
  }
  return val;
}

/*
 * Compare elements of both arrays of same length n
 * and return TRUE if they are all the same
 */
bool vec_is_equal(const int n, const double *x1, const double *x2) {

  for (int i = 0; i < n; i++) {
    if (x1[i] != x2[i]) {
      return false;
    }
  }
  return true;
}

/*
 * Returns constant vector of val value from 1 to n
 */
void const_vec(const int n, const double val, double *vec) {

  int i;
  for (i = 0; i < n; i++) {
    vec[i] = val;
  }
}

/*
 * Returns the inner product between two vectors of size n
 */
double inner_prod(const int n, const double *a1, const double *a2) {

  int i;
  double val = 0.0;
  for (i = 0; i < n; i++) {
    val += a1[i] * a2[i];
  }

  return val;
}

/*
 * Returns the weighted inner product between two vectors of size given in last
 * argument
 */
double inner_w_prod(const int size, const double *w, const double *a1,
                    const double *a2) {

  int i;
  double val = 0.0;
  for (i = 0; i < size; i++) {
    val += a1[i] * w[i] * a2[i];
  }
  return val;
}

/*
 * Returns the inverse weighted inner product between two vectors of size given
 * in last argument
 */
double inner_winv_prod(const int size, const double *w, const double *a1,
                       const double *a2) {

  int i;
  double val = 0.0;
  for (i = 0; i < size; i++) {
    val += a1[i] * a2[i] / w[i];
  }
  return val;
}

/*
 * Makes a2 array equal to a1
 *
 */
void make_equal(const int n, const double *a1, double *a2) {

  for (int i = 0; i < n; i++)
    a2[i] = a1[i];
}

/*
 * Returns a1 + a2 vector into a1, assuming both have the same length n
 */
void vec_plus(const int n, const double *a2, double *a1) {

  for (int i = 0; i < n; i++) {
    a1[i] += a2[i];
  }
}

/*
 * Returns a1 - a2 vector into a1, assuming both have the same length n
 */
void vec_minus(const int n, const double *a2, double *a1) {

  for (int i = 0; i < n; i++) {
    a1[i] -= a2[i];
  }
}

/*
 * Return time of day as micro seconds
 */
long get_time() {
  struct timeval tv;
  if (gettimeofday(&tv, (struct timezone *)0) == 0)
    return (tv.tv_sec * 1000 * 1000 + tv.tv_usec); // us

  return 0.;
}

/*!*****************************************************************************
 *******************************************************************************
 \note  sign
 \date  02/25/91

 \remarks

 calculates the SIGN function

 *******************************************************************************
 Function Parameters: [in]=input,[out]=output

 \param[in]     expr : argument for SIGN
 return	: SIGN(expr)

 ******************************************************************************/
double sign(double expr) {
  if (expr > 0)
    return (1.0);
  else if (expr < 0)
    return (-1.0);
  else
    return (0.0);
}
