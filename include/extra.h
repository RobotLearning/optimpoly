/*
 * extra.h
 *
 *  Created on: Jun 21, 2016
 *      Author: okoc
 */

#ifndef EXTRA_H_
#define EXTRA_H_

#define DOF 7 //TODO: put ifndef macro preprocessor command

#include <stdio.h>
#include <stdlib.h>
#include "SL.h"

// utility methods, zero indexed
long get_time();
void vec_minus(double *a1, const double *a2);
void vec_plus(double *a1, const double *a2);
double inner_prod(const double *a1, const double *a2);
void const_vec(const int n, const double val, double * vec);
void print_optim_vec(double *x);

// io method from SL, one indexed
void load_vec_into_mat(Matrix mat, int m, int n, char name[]);

#endif /* EXTRA_H_ */
