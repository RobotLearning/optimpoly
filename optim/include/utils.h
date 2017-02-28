/*
 * utils.h
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

#ifndef NDOF
#define NDOF 7
#endif

#ifndef NCART
#define NCART 3
#endif

double sign(double expr);

// utility methods, zero indexed
long get_time();
void vec_minus(double *a1, const double *a2);
void vec_plus(double *a1, const double *a2);
double inner_prod(const double *a1, const double *a2);
void const_vec(const int n, const double val, double * vec);
void print_optim_vec(const double *x);
void print_mat_size(const char *comm, Matrix a, int nr, int nc);
double max_abs_array(const double *x, const int length);
double max_array(const double *x, const int length);

// utility methods, one indexed
Vector my_vector(int nl, int nh);
double vec_mult_inner(Vector a, Vector b);
void vec_mult_scalar(Vector a, double scalar, Vector c);
int vec_sub(Vector a, Vector b, Vector c);
int vec_add(Vector a, Vector b, Vector c);
Matrix my_matrix(int nrl, int nrh, int ncl, int nch);
int mat_vec_mult(Matrix a, Vector b, Vector c);
void my_free_vector(Vector vec, int nl, int nh);
int vec_equal(Vector a, Vector c);
int mat_mult(Matrix a, Matrix b, Matrix c);
int mat_equal(Matrix a, Matrix c);
void my_free_matrix(Matrix mat, int nrl, int nrh, int ncl, int nch);

// io method from SL, one indexed
int find_keyword(FILE *fp, char *name);

#endif /* UTILS_H_ */
