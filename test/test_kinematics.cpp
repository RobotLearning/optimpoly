/*
 * kinematics.cpp
 *
 * Testing the kinematics inside the optimization
 *
 *  Created on: Feb 20, 2017
 *      Author: okoc
 */

#include <boost/test/unit_test.hpp>
#include <iostream>
#include <armadillo>
#include "kinematics.h"

using namespace arma;

/*
 * Comparing with MATLAB the racket states calculated given q, qd
 *
 */
BOOST_AUTO_TEST_CASE( test_kinematics_calculations ) {

	cout << "Comparing racket state calculations with MATLAB..." << endl;

	static double q[NDOF], qdot[NDOF], pos[NCART], vel[NCART], normal[NCART];
	for (int i = 0; i < NDOF; i++) {
		q[i] = 1.0;
	}

	calc_racket_state(q,qdot,pos,vel,normal);

	vec3 x(pos);
	vec3 xdot(vel);
	vec3 n(normal);

	vec3 x_MATLAB = {-0.8309, 0.0861, -0.6672};
	vec3 xdot_MATLAB = zeros<vec>(3);
	vec3 n_MATLAB = {0.3857, 0.0174, -0.9225};

	/*cout << "Kin values:" << endl;
	cout << x << xdot << n << endl;
	cout << "MATLAB values:" << endl;
	cout << x_MATLAB << xdot_MATLAB << n_MATLAB << endl;*/

	BOOST_TEST(approx_equal(x,x_MATLAB,"absdiff", 0.002));
	BOOST_TEST(approx_equal(xdot,xdot_MATLAB,"absdiff", 0.002));
	BOOST_TEST(approx_equal(n,n_MATLAB,"absdiff", 0.002));

}



