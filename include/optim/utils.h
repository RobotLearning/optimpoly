/**
 * @file utils.h
 *
 * @brief All C-legacy utility methods that we are using for
 * running NLOPT optimizations are stored here.
 *
 *  Created on: Feb 7, 2017
 *      Author: okoc
 */

#ifndef UTILS_H_
#define UTILS_H_

#include <armadillo>
#include "json.hpp"

// utility methods, zero indexed
long get_time();
void vec_minus(const int n, const double *a2, double *a1);
void vec_plus(const int n, const double *a2, double *a1);
void make_equal(const int n, const double *a1, double *a2);
bool vec_is_equal(const int n, const double *x1, const double *x2);
double inner_prod(const int size, const double *a1, const double *a2);
double inner_w_prod(const int size, const double *w, const double *a1, const double *a2);
double inner_winv_prod(const int size, const double *w, const double *a1, const double *a2);
void const_vec(const int n, const double val, double * vec);
void print_optim_vec(const double *x);
double max_abs_array(const double *x, const int length);
double max_array(const double *x, const int length);
double sign(double expr);

// new functions
arma::mat json2mat(const nlohmann::json & jobj);
arma::vec json2vec(const nlohmann::json & jobj);

#endif /* UTILS_H_ */
