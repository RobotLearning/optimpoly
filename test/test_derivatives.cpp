/*
 * test_derivatives.cpp
 *
 * Here we're testing the derivatives of our
 * optimization cost and constraint functions.
 *
 *  Created on: Jun 5, 2017
 *      Author: okan
 */

#ifndef BOOST_TEST_MODULE
#define BOOST_TEST_MODULE test_table_tennis
#endif

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/data/monomorphic.hpp>

#include <armadillo>
#include "optim.h"
#include "kinematics.h"
#include "player.hpp"

#define OPTIM_DIM 2*NDOF + 1

/*
 * This is the constraint that makes sure we hit the ball
 */
static void kinematics_eq_constr(unsigned m, double *result, unsigned n,
		                  const double *x, double *grad, void *my_function_data) {

	static double pos[NCART];
	static double qfdot[NDOF];
	static double vel[NCART];
	static double normal[NCART];
	static double qf[NDOF];

	HittingPlane * vhp = (HittingPlane*)my_function_data;

	// extract state information from optimization variables
	for (int i = 0; i < NDOF; i++) {
		qf[i] = x[i];
		qfdot[i] = x[i+NDOF];
	}

	// compute the actual racket pos,vel and normal
	calc_racket_state(qf,qfdot,pos,vel,normal);

	// deviations from the desired racket frame
	for (int i = 0; i < NCART; i++) {
		result[i] = pos[i] - vhp->param_des->racket_pos(i);
		result[i + NCART] = vel[i] - vhp->param_des->racket_vel(i);
		result[i + 2*NCART] = normal[i] - vhp->param_des->racket_normal(i);
	}
}

static double calc_max_diff(const double mat1[EQ_CONSTR_DIM][OPTIM_DIM], const double mat2[EQ_CONSTR_DIM][OPTIM_DIM],
		                    int m1, int m2, int n1, int n2) {

	double maxdiff = 0.0;
	double absdiff;
	for (int i = m1; i < m2; i++) {
		for (int j = n1; j < n2; j++) {
			//printf("jac[%d,%d] = %f, numjac[%d,%d] = %f\n", i,j, mat1[i][j], i,j, mat2[i][j]);
			absdiff = fabs(mat1[i][j] - mat2[i][j]);
			if (absdiff > maxdiff) {
				maxdiff = absdiff;
			}
		}
	}
	return maxdiff;
}

/*
 * Compare derivatives with numerical differentiation
 */
BOOST_AUTO_TEST_CASE( test_deriv_opt ) {

	double lb[OPTIM_DIM], ub[OPTIM_DIM];
	set_bounds(lb,ub,0.0,1.0);
	double q0[NDOF] = {1.0, -0.2, -0.1, 1.8, -1.57, 0.1, 0.3};
	HittingPlane opt = HittingPlane(q0,lb,ub);
	vec6 ball_pred = {0.0, -0.3, -1.0, -1.0, 3.0, -2.0};
	vec2 ball_land_des = {0.0, -3.0};
	optim_des pred_params;
	calc_racket_strategy(ball_pred,ball_land_des,0.8,pred_params);
	opt.set_des_params(&pred_params);

	double constr1[EQ_CONSTR_DIM], constr2[EQ_CONSTR_DIM];
	double x[OPTIM_DIM];
	double xdiff[OPTIM_DIM];
	double deriv[EQ_CONSTR_DIM][OPTIM_DIM];
	double num_deriv[EQ_CONSTR_DIM][OPTIM_DIM];
	double jac[NCART][NDOF];
	double h = 0.001;

	for (int i = 0; i < NDOF; i++) {
		x[i] = xdiff[i] = q0[i];
		x[i+NDOF] = xdiff[i+NDOF] = 0.0;
	}

	/*
	 * Calculate exact derivatives
	 */
	get_jacobian(x,jac);

	/*
	 * Calculate numerical derivatives
	 */
	for (int i = 0; i < OPTIM_DIM; i++) {
		xdiff[i] = x[i] + h;
		kinematics_eq_constr(EQ_CONSTR_DIM,constr1,OPTIM_DIM,xdiff,nullptr,(void*)&opt);
		xdiff[i] = x[i] - h;
		kinematics_eq_constr(EQ_CONSTR_DIM,constr2,OPTIM_DIM,xdiff,nullptr,(void*)&opt);
		for (int j = 0; j < EQ_CONSTR_DIM; j++) {
			num_deriv[j][i] = (constr1[j] - constr2[j]) / (2*h);
		}
		xdiff[i] = x[i];
	}

	/*
	 * Fill the big exact derivative matrix
	 */
	for (int i = 0; i < NCART; i++) {
		for (int j = 0; j < NDOF; j++) {
			deriv[i][j] = jac[i][j];
		}
	}
	for (int i = NCART; i < 2*NCART; i++) {
		for (int j = NDOF; j < 2*NDOF; j++) {
			deriv[i][j] = jac[i-NCART][j-NDOF];
		}
	}

	/*
	 * Calculate the maximum difference between numerical and exact derivative matrices (big jacobian)
	 */
	BOOST_TEST(calc_max_diff(deriv,num_deriv,0,NCART,0,NDOF) < 1e-3);
	BOOST_TEST(calc_max_diff(deriv,num_deriv,NCART,2*NCART,NDOF,2*NDOF) < 1e-3);
}
