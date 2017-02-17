/*
 * utils.c
 *
 *  Created on: Feb 7, 2017
 *      Author: okoc
 */

#include "utils.h"
#include "stdlib.h"
#include "stdio.h"

/*****************************************************************************
  utility program vector: allocates a double vector with range nl...nh
 ****************************************************************************/
Vector
my_vector(int nl, int nh)

{
	double *v;

	if (nl == 1) {

		v = (double *) calloc((size_t) (nh-nl+1+1),sizeof(double));
		if (v==NULL)
			printf("allocation failure in vector()");
		v[0] = nh-nl+1;
		return v;

	} else {

		v = (double *) calloc((size_t) (nh-nl+1),sizeof(double));
		if (v==NULL)
			printf("allocation failure in vector()");
		return v-nl;

	}

}

/*!*****************************************************************************
 *******************************************************************************
 \note  vec_mult_inner
 \date  August 17, 92

 \remarks

 inner product a' * b, return result
vector indices start at "1".

 *******************************************************************************
Parameters:  (i/o = input/output)

 \param[in]     a		 : vector a
 \param[in]     b		 : vector b

 ******************************************************************************/
double
vec_mult_inner(Vector a, Vector b)

{

	double aux = 0;
	int i;

	if (a[NR] != b[NR]) {
		printf("Incompatible vectors in vec_mult_inner\n");
		return 0.0;
	}

	for (i=1; i<=a[NR]; ++i) aux += a[i] * b[i];

	return aux;

}

/*!*****************************************************************************
 *******************************************************************************
\note  vec_mult_scalar
\date  August 17, 92

\remarks

product of vector a  * scalar = c
vector indices start at "1".

 *******************************************************************************
Parameters:  (i/o = input/output)

 \param[in]     a		 : vector a
 \param[in]     scalar           : value by which a is to be multiplied with
 \param[out]    c		 : vector c (result)

 ******************************************************************************/
void
vec_mult_scalar(Vector a, double scalar, Vector c)

{

	double aux = 0;
	int i;

	for (i=1; i<=a[NR]; ++i) c[i] = a[i] * scalar;

}

/*!*****************************************************************************
 *******************************************************************************
 \note  vec_sub
 \date  August 17, 92

 \remarks

 subtracts two arbitrary (compatible) vectors a - b, assuming the
 vector indices start at "1".

 *******************************************************************************
 Parameters:  (i/o = input/output)

 \param[in]     a		 : vector a
 \param[in]     b		 : vector b
 \param[out]    c		 : result of addition

 ******************************************************************************/
int
vec_sub(Vector a, Vector b, Vector c)
{

	int     i;

	if (a[NR] != b[NR] || a[NR] != c[NR]) {
		printf("Incompatible vectors in vec_sub\n");
		return FALSE;
	}

	for (i=1; i <= a[NR]; ++i) {
		c[i] = a[i] - b[i];
	}

	return TRUE;

}

/*!*****************************************************************************
 *******************************************************************************
 \note  vec_add
 \date  August 17, 92

 \remarks

 adds two arbitrary (compatible) vectors a + b, assuming the
 vector indices start at "1".

 *******************************************************************************
 Parameters:  (i/o = input/output)

 \param[in]     a		 : vector a
 \param[in]     b		 : vector b
 \param[out]    c		 : result of addition

 ******************************************************************************/
int
vec_add(Vector a, Vector b, Vector c)
{

	int     i;

	if (a[NR] != b[NR] || a[NR] != c[NR]) {
		printf("Incompatible vectors in vec_add\n");
		return FALSE;
	}

	for (i=1; i <= a[NR]; ++i) {
		c[i] = a[i] + b[i];
	}

	return TRUE;

}

/******************************************************************************
  utility program my_matrix: note: this is a modified version of
  the normal matrix() numerical recipe version. It allocates
  the memory for the matrix in one chunk, such that it can
  be treated as consecutive memory in write() or DSP programs

 *****************************************************************************/
Matrix
my_matrix(int nrl, int nrh, int ncl, int nch)
{

	int      i;
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

	for(i=nrl;i<=nrh;i++) {
		m[i]=(double *) &(chunk[(i-nrl)*(nch-ncl+1)]);
		m[i] -= ncl;
	}

	return m;

}

/*!*****************************************************************************
 *******************************************************************************
 \note  mat_vec_mult
 \date  August 17, 92

 \remarks

 multiplies a matrix with a vector: a * b, assuming indices
 start at "1".
 Note: The program can also cope with passing the same vector as
 factor and result.

 *******************************************************************************
 Parameters:  (i/o = input/output)

 \param[in]     a		 : matrix a
 \param[in]     b		 : vector b
 \param[out]    c		 : result of multipliciation

 ******************************************************************************/
int
mat_vec_mult(Matrix a, Vector b, Vector c) {

	int     i,j,m;
	Vector  temp;
	int     ac,ar,br,cr;

	ac = a[0][NC];
	ar = a[0][NR];
	br = b[NR];
	cr = c[NR];

	if (ac != br) {
		printf("Matrix and vector are incompatible.\n");
		return FALSE;
	}

	/* bugfix (adsouza July 10, 2002)
  if (cr != br) {
    printf("Input and output vector are incompatible.\n");
    return FALSE;
  }
  end of old version */

	if (cr != ar) {
		printf("Input and output vector are incompatible.\n");
		return FALSE;
	}


	if (b == c) {
		temp = my_vector(1,ar);
	} else {
		temp = c;
	}

	for (i=1; i <= ar; ++i) {
		temp[i]=0;
		for (j=1; j <= br; ++j){
			temp[i] += a[i][j] * b[j];
		}
	}



	if (b == c) {
		vec_equal(temp,c);
		my_free_vector(temp,1,ar);
	}

	return TRUE;

}

void
my_free_vector(Vector vec, int nl, int nh)
{
	if (nl == 1) {
		free((char*) (vec));
	} else {
		free((char*) (vec+nl));
	}
}

/*!*****************************************************************************
 *******************************************************************************
 \note  vec_equal
 \date  August 17, 92

 \remarks

 set vector c = a
 vector indices start at "1".

 *******************************************************************************
 Parameters:  (i/o = input/output)

 \param[in]     a		 : vector a
 \param[out]    c		 : result of assignment

 ******************************************************************************/
int
vec_equal(Vector a, Vector c)

{
	int i;

	if (a[NR] != c[NR]) {
		printf("Incompatible vectors in vec_equal\n");
		return FALSE;
	}

	for (i=1; i<=a[NR]; ++i) c[i] = a[i];

	return TRUE;
}

/*!*****************************************************************************
 *******************************************************************************
 \note  mat_mult
 \date  August 17, 92

 \remarks

 multiplies two arbitrary (compatible) matrices a * b

 *******************************************************************************
 Parameters:  (i/o = input/output)

 \param[in]     a		 : matrix a
 \param[in]     b		 : matrix b
 \param[out]    c		 : result of multipliciation

 ******************************************************************************/
int
mat_mult(Matrix a, Matrix b, Matrix c)
{
	int      i,j,m;
	Matrix   temp;
	int      ar,ac,br,bc;
	int      type;

	ar = a[0][NR];
	ac = a[0][NC];

	br = b[0][NR];
	bc = b[0][NC];


	/* check whether input and output matrices are different */

	if (a == c || b == c) {
		temp = my_matrix(1,ar,1,bc);
	}
	else {
		temp = c;
	}



	for (i=1; i <= ar; ++i) {
		for (j=1; j <= bc; ++j){
			temp[i][j]=0;
			for (m=1; m <= br; ++m){
				temp[i][j] += a[i][m] * b[m][j];
			}
		}
	}

	if (a == c || b == c) {
		mat_equal(temp,c);
		my_free_matrix(temp,1,ar,1,bc);
	}

	return TRUE;

}

/*****************************************************************************
  utility program my_free_matrix; adjusted to my special matrix() program
 ***************************************************************************/
void
my_free_matrix(Matrix mat, int nrl, int nrh, int ncl, int nch)

{
	int i;

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
 \note  mat_equal
 \date  August 17, 92

 \remarks

 set matrix c = a
 matrix indices start at "1".
 Note: if a vector is passed, it must be the pointer to the vector;
 everything is handled as a matrix!

 *******************************************************************************
 Parameters:  (i/o = input/output)

 \param[in]     a		 : matrix a
 \param[out]    c		 : result of addition

 ******************************************************************************/
int
mat_equal(Matrix a, Matrix c) {
	int i,j;

	for (i=1; i <= a[0][NR]; ++i) {
		for (j=1; j <= a[0][NC]; ++j){
			c[i][j] = a[i][j];
		}
	}

	return TRUE;
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
double
sign(double expr)
{
	if (expr > 0)
		return (1.0);
	else
		if (expr < 0)
			return (-1.0);
		else
			return (0.0);
}
