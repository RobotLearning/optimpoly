/*
 * utils.c
 *
 *  Created on: Feb 7, 2017
 *      Author: okoc
 */

#include "constants.h"
#include "utils.h"
#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "sys/time.h"
#include "math.h"

using namespace const_tt;

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
		printf("%.2f  ", x[i+NDOF]);
	}
	printf("]\n");
	printf("T = %.2f\n", x[2*NDOF]);
}


/*!*****************************************************************************
 *******************************************************************************
 \note  print_mat_size
 \date  August 17, 92

 \remarks

 just prints the given matrix with limits

 *******************************************************************************
 Parameters:  (i/o = input/output)

 \param[in]     com   : comment
 \param[in]     a     : the matrix to be printed
 \param[in]     nr    : number of rows
 \param[in]     nc    : number of columns

 ******************************************************************************/
void print_mat_size(const char *comm, Matrix a, int nr, int nc) {

	printf("Matrix %s :\n",comm);

	for (int i = 0; i < nr; i++) {
		printf("          ");
		for (int j = 0; j < nc; j++) {
			printf("% 8.4f ",a[i][j]);
		}
		printf("\n");
	}
	printf("\n");

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

/*!*****************************************************************************
 *******************************************************************************
 \note  find_keyword
 \date  October, 1995

 \remarks

 find a key word in the file given by fp and positions the read cursor
 right after the keyword

 *******************************************************************************
 Function Parameters: [in]=input,[out]=output

 \param[in]     fp       : file pointer
 \param[in]     word     : keyword to be found

 ******************************************************************************/
int find_keyword(FILE *fp, char *name) {

	int  i;
	int  rc = TRUE;
	char string[strlen(name)*2];
	int  l;
	char sep[]={' ','\n',':',',',';','=','\t','\0'};

	rewind(fp);
	l  = strlen(name);
	i  = 0;

	while (rc != EOF) {

		rc=fgetc(fp);
		if ( rc != EOF ) {

			string[i++] = rc;
			string[i]   = '\0';

			if ( strstr(string,name) != NULL) {
				// wait for one more character to judge whether this string
				// has the correct end delimiters
				if (strchr(sep,string[i-1]) != NULL) {
					// now check for preceeding delimiter

					if (i-l-2 < 0) // this means "name" was the first string in file
						return TRUE;
					else if (strchr(sep,string[i-l-2]) != NULL) //otherwise check delim
						return TRUE;
				}
			}

			if (i >= 2*l-1) {
				strcpy(string,&(string[i-l]));
				i = strlen(string);
			}

		}

	}

	return FALSE;

}

/*
 * Compare elements of both arrays of same length n
 * and return TRUE if they are all the same
 */
int vec_is_equal(const int n, const double *x1, const double *x2) {

	for (int i = 0; i < n; i++) {
		if (x1[i] != x2[i]) {
			return FALSE;
		}
	}
	return TRUE;
}

/*
 * Returns constant vector of val value from 1 to n
 */
void const_vec(const int n, const double val, double * vec) {

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
		val += a1[i]*a2[i];
	}

	return val;
}

/*
 * Returns the weighted inner product between two vectors of size given in last argument
 */
double inner_w_prod(const int size, const double *w, const double *a1, const double *a2) {

	int i;
	double val = 0.0;
	for (i = 0; i < size; i++) {
		val += a1[i]*w[i]*a2[i];
	}
	return val;
}

/*
 * Returns the inverse weighted inner product between two vectors of size given in last argument
 */
double inner_winv_prod(const int size, const double *w, const double *a1, const double *a2) {

	int i;
	double val = 0.0;
	for (i = 0; i < size; i++) {
		val += a1[i]*a2[i]/w[i];
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
		return (tv.tv_sec*1000*1000 + tv.tv_usec);  //us

	return 0.;
}


/*****************************************************************************
  utility program vector: allocates a double vector with range nl...nh
 ****************************************************************************/
Vector my_vector(int nl, int nh) {
	double *v;

	if (nl == 1) {

		v = (double *) calloc((size_t) (nh-nl+1+1),sizeof(double));
		if (v == NULL)
			printf("allocation failure in vector()");
		v[0] = nh-nl+1;
		return v;

	} else {

		v = (double *) calloc((size_t) (nh-nl+1),sizeof(double));
		if (v == NULL)
			printf("allocation failure in vector()");
		return v-nl;

	}

}

/******************************************************************************
  utility program my_matrix: note: this is a modified version of
  the normal matrix() numerical recipe version. It allocates
  the memory for the matrix in one chunk, such that it can
  be treated as consecutive memory in write() or DSP programs

 *****************************************************************************/
Matrix my_matrix(int nrl, int nrh, int ncl, int nch) {

	double **m;
	double  *chunk;
	int      info = FALSE;

	if (nrl==1 && ncl == 1) {
		info = TRUE;
	}

	m = (double **) calloc((size_t) (nrh-nrl+1+info),sizeof(double*));

	if (info) {

		m[0] = (double *) calloc((size_t) N_MAT_INFO,sizeof(double));
		m[0][NR]       = nrh-nrl+1;
		m[0][NC]       = nch-ncl+1;

	} else {

		m -= nrl;

	}

	chunk = (double *) calloc( (size_t) (nrh-nrl+1) * (nch-ncl+1),sizeof(double));

	for(int i = nrl; i <= nrh;i ++) {
		m[i] = (double *) &(chunk[(i-nrl)*(nch-ncl+1)]);
		m[i] -= ncl;
	}

	return m;

}

void my_free_vector(Vector vec, int nl, int nh) {
	if (nl == 1) {
		free((char*) (vec));
	} else {
		free((char*) (vec+nl));
	}
}


/*****************************************************************************
  utility program my_free_matrix; adjusted to my special matrix() program
 ***************************************************************************/
void my_free_matrix(Matrix mat, int nrl, int nrh, int ncl, int nch) {

	free((char*) &(mat[nrl][ncl]));
	if (nrl==1 && ncl==1) {
		free((char*) mat[0]);
		free((char*) mat);
	} else {
		free((char*) (mat+nrl));
	}

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
	else
		if (expr < 0)
			return (-1.0);
		else
			return (0.0);
}
