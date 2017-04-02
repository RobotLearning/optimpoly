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

#include "stdio.h"

#define N_MAT_INFO 3
#define NR       0
#define NC       1

#define TRUE 1
#define FALSE 0

#define sqr(x)  ((x)*(x))

typedef double* Vector;
typedef double** Matrix;

double sign(double expr);

// utility methods, zero indexed
long get_time();
void vec_minus(const int n, const double *a2, double *a1);
void vec_plus(const int n, const double *a2, double *a1);
void make_equal(const int n, const double *a1, double *a2);
int vec_is_equal(const int n, const double *x1, const double *x2);
double inner_prod(const int size, const double *a1, const double *a2);
double inner_w_prod(const int size, const double *w, const double *a1, const double *a2);
double inner_winv_prod(const int size, const double *w, const double *a1, const double *a2);
void const_vec(const int n, const double val, double * vec);
void print_optim_vec(const double *x);
void print_mat_size(const char *comm, Matrix a, int nr, int nc);
double max_abs_array(const double *x, const int length);
double max_array(const double *x, const int length);

// utility methods, one indexed
Vector my_vector(int nl, int nh);
Matrix my_matrix(int nrl, int nrh, int ncl, int nch);
void my_free_vector(Vector vec, int nl, int nh);
void my_free_matrix(Matrix mat, int nrl, int nrh, int ncl, int nch);

// io method from SL, one indexed
int find_keyword(FILE *fp, char *name);

#endif /* UTILS_H_ */
