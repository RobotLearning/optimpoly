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
#include "player.hpp"

static double penalize_dist_to_limits(unsigned n, const double *x, double *grad, void *my_func_params) {

	double cost = 0.0;
	HittingPlane * vhp = (HittingPlane*)my_func_params;

	if (grad) {
		for (int i = 0; i < NDOF; i++) {
			grad[i] = 2 * (x[i] - vhp->limit_avg[i]);
			grad[i+NDOF] = 2 * x[i+NDOF];
		}
	}

	for (int i = 0; i < NDOF; i++) {
		cost += pow(x[i] - vhp->limit_avg[i],2);
		cost += pow(x[i + NDOF], 2);
	}
	return cost;
}

/*
 * Compare derivatives with numerical differentiation
 */
BOOST_AUTO_TEST_CASE( test_deriv_opt ) {

	const int OPTIM_DIM = 2*NDOF;
	double lb[2*NDOF+1], ub[2*NDOF+1];
	set_bounds(lb,ub,0.0,1.0);
	double q0[NDOF] = {1.0, -0.2, -0.1, 1.8, -1.57, 0.1, 0.3};
	HittingPlane opt = HittingPlane(q0,lb,ub);

	double x[OPTIM_DIM];
	double xdiff[OPTIM_DIM];
	double grad[OPTIM_DIM];
	double numgrad[OPTIM_DIM];
	double h = 0.001;

	for (int i = 0; i < NDOF; i++) {
		x[i] = xdiff[i] = q0[i];
		x[i+NDOF] = xdiff[i+NDOF] = 0.0;
	}
	penalize_dist_to_limits(OPTIM_DIM,x,grad,(void*)&opt);

	double val1,val2;
	for (int i = 0; i < OPTIM_DIM; i++) {
		xdiff[i] = x[i] + h;
		val1 = penalize_dist_to_limits(OPTIM_DIM,xdiff,nullptr,(void*)&opt);
		xdiff[i] = x[i] - h;
		val2 = penalize_dist_to_limits(OPTIM_DIM,xdiff,nullptr,(void*)&opt);
		numgrad[i] = (val1 - val2)/(2*h);
		xdiff[i] = x[i];
	}

	double maxdiff = 0.0;
	double absdiff;
	for (int i = 0; i < OPTIM_DIM; i++) {
		//printf("grad[%d] = %f, numgrad[%d] = %f\n", i, grad[i], i, numgrad[i]);
		absdiff = fabs(grad[i] - numgrad[i]);
		if (absdiff > maxdiff) {
			maxdiff = absdiff;
		}
	}
	BOOST_TEST(maxdiff < 1e-3);
}





