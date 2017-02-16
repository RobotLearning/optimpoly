/*
 * utils.c
 *
 *  Created on: Feb 7, 2017
 *      Author: okoc
 */

#include "utils.h"
#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "sys/time.h"

/*
 * Prints the 2*DOF + 1 dimensional solution in user-friendly format
 */
void print_optim_vec(double *x) {

	int i;
	printf("qf = [");
	for (i = 0; i < DOF; i++) {
		printf("%.2f  ", x[i]);
	}
	printf("]\n");
	printf("qfdot = [");
	for (i = 0; i < DOF; i++) {
		printf("%.2f  ", x[i+DOF]);
	}
	printf("]\n");
	printf("T = %.2f\n", x[2*DOF]);
}

/*
 * Read a vector into a matrix (e.g. saved from MATLAB)
 * Necessary for loading large matrices into memory
 *
 */
void load_vec_into_mat(Matrix mat, int m, int n, char name[]) {

	FILE * fid;
	int i, j;
	static char fileName[100];
	static char dirName[] = "../robolab/barrett/saveData";

	sprintf(fileName,"%s//%s.txt",dirName,name);
    fid = fopen(fileName,"r");

    for (i = 1; i <= n; i++) {
    	for (j = 1; j <= m; j++) {
    		fscanf(fid, "%lf", &(mat[j][i]));
    		//printf("M[%d][%d] = %f \n", j, i, mat[j][i]);
    	}
    }

	fclose(fid);
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
 * Returns constant vector of val value from 1 to n
 */
void const_vec(const int n, const double val, double * vec) {

	int i;
	for (i = 0; i < n; i++) {
		vec[i] = val;
	}
}

/*
 * Returns the inner product between two vectors of size DOF
 */
double inner_prod(const double *a1, const double *a2) {

	int i;
	double val = 0.0;
	for (i = 0; i < DOF; i++) {
		val += a1[i]*a2[i];
	}

	return val;
}

/*
 * Returns a1 + a2 vector into a1, assuming both have dof = 7 length
 */
void vec_plus(double *a1, const double *a2) {

	int i;
	for (i = 0; i < DOF; i++) {
		a1[i] = a1[i] + a2[i];
	}
}

/*
 * Returns a1 - a2 vector into a1, assuming both have dof = 7 length
 */
void vec_minus(double *a1, const double *a2) {

	int i;
	for (i = 0; i < DOF; i++) {
		a1[i] = a1[i] - a2[i];
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
vec_mult_inner(Vector a, Vector b) {

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
void vec_mult_scalar(Vector a, double scalar, Vector c) {

	for (int i = 1; i <= a[NR]; ++i)
		c[i] = a[i] * scalar;

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
int vec_sub(Vector a, Vector b, Vector c) {

	if (a[NR] != b[NR] || a[NR] != c[NR]) {
		printf("Incompatible vectors in vec_sub\n");
		return FALSE;
	}

	for (int i = 1; i <= a[NR]; ++i) {
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
int vec_add(Vector a, Vector b, Vector c) {

	if (a[NR] != b[NR] || a[NR] != c[NR]) {
		printf("Incompatible vectors in vec_add\n");
		return FALSE;
	}

	for (int i = 1; i <= a[NR]; ++i) {
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

	if (cr != ar) {
		printf("Input and output vector are incompatible.\n");
		return FALSE;
	}

	if (b == c) {
		temp = my_vector(1,ar);
	} else {
		temp = c;
	}

	for (int i = 1; i <= ar; ++i) {
		temp[i] = 0;
		for (int j = 1; j <= br; ++j){
			temp[i] += a[i][j] * b[j];
		}
	}

	if (b == c) {
		vec_equal(temp,c);
		my_free_vector(temp,1,ar);
	}

	return TRUE;

}

void my_free_vector(Vector vec, int nl, int nh) {
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
int vec_equal(Vector a, Vector c) {
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
int mat_mult(Matrix a, Matrix b, Matrix c) {

	Matrix   temp;
	int      ar, br, bc;

	ar = a[0][NR];
	br = b[0][NR];
	bc = b[0][NC];


	/* check whether input and output matrices are different */

	if (a == c || b == c) {
		temp = my_matrix(1,ar,1,bc);
	}
	else {
		temp = c;
	}

	for (int i = 1; i <= ar; ++i) {
		for (int j = 1; j <= bc; ++j) {
			temp[i][j] = 0;
			for (int m = 1; m <= br; ++m) {
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
int mat_equal(Matrix a, Matrix c) {
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
double sign(double expr) {
	if (expr > 0)
		return (1.0);
	else
		if (expr < 0)
			return (-1.0);
		else
			return (0.0);
}
