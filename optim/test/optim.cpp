/*
 * optim.cpp
 *
 * Unit Tests for polynomial optimization
 *
 *  Created on: Feb 17, 2017
 *      Author: okoc
 */

#define BOOST_TEST_MODULE test_optim
#include <boost/test/unit_test.hpp>
#include "SL.h"
#include "constants.h"
#include "optimpoly.h"
#include "utils.h"
#include "kinematics.h"


BOOST_AUTO_TEST_CASE(test_nlopt_optim) {

	Vector ballLand = my_vector(1,CART);
	double landTime;
	double b0[CART];
	double v0[CART];
	double x[OPTIM_DIM]; /* initial guess for optim */

	Matrix lookupTable = my_matrix(1, LOOKUP_TABLE_SIZE, 1, LOOKUP_COLUMN_SIZE);
	load_lookup_table(lookupTable);
	set_des_land_param(ballLand,&landTime);

	/* initialize ball and racket */
	// predict for T_pred seconds
	init_ball_state(b0,v0);
	predict_ball_state(b0,v0);
	calc_racket_strategy(ballLand,landTime);

	/* run NLOPT opt algorithm here */
	lookup(lookupTable,b0,v0,x); // initialize solution
	nlopt_optim_poly_run(x);

	// test the lookup value to see if constraint is not violated
	printf("================== TEST ==================\n");
	printf("Lookup values:\n");
	lookup(lookupTable,b0,v0,x);
	test_optim(x);


}




