/*
 * extra.c
 *
 * Matrix manipulation, and i/o functions are located here for convenience
 *
 *  Created on: Jun 21, 2016
 *      Author: okoc
 */

#include "extra.h"

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
