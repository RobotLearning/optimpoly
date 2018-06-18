/*
 * utils.h
 *
 *  Created on: Feb 7, 2017
 *      Author: okoc
 */

#ifndef UTILS_H_
#define UTILS_H_

#define N_MAT_INFO 3
#define NR       0
#define NC       1

#define TRUE 1
#define FALSE 0

#define sqr(x)  ((x)*(x))

typedef double* Vector;
typedef double** Matrix;

long get_time();
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
double sign(double expr);

#endif /* UTILS_H_ */
