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
#include "optim.h"
#include "kinematics.h"
#include "kinematics.hpp"
#include "player.hpp"

using namespace arma;
using namespace optim;

#define OPTIM_DIM 2*NDOF + 1
static void print_mat(const double mat[NCART][NDOF]);
static void cross_prods(const double mat[NCART][NDOF],
                          const double v[NCART],
                          double out[NCART][NDOF]);
static double calc_max_diff(const double mat1[EQ_CONSTR_DIM][OPTIM_DIM],
                            const double mat2[EQ_CONSTR_DIM][OPTIM_DIM],
                            int m1, int m2, int n1, int n2);
static void kinematics_eq_constr(double *result,
                                 const double *x,
                                 void *my_function_data);

/*
 * Compare derivatives with numerical differentiation
 */
void test_kin_deriv() {

    BOOST_TEST_MESSAGE("\nComparing kinematics derivatives with numerical diff...");
	static double lb[OPTIM_DIM];
	static double ub[OPTIM_DIM];
	set_bounds(lb,ub,0.0,1.0);
	double q0[NDOF] = {1.0, -0.2, -0.1, 1.8, -1.57, 0.1, 0.3};
	HittingPlane opt = HittingPlane(q0,lb,ub);
	vec6 ball_pred = {0.0, -0.3, -1.0, -1.0, 3.0, -2.0};
	vec2 ball_land_des = {0.0, -3.0};
	optim_des pred_params;
	calc_racket_strategy(ball_pred,ball_land_des,0.8,pred_params);
	opt.set_des_params(&pred_params);

	static double constr1[EQ_CONSTR_DIM];
	static double constr2[EQ_CONSTR_DIM];
	static double racket_pos[NCART];
	static double racket_normal[NCART];
	static double x[OPTIM_DIM];
	static double xdiff[OPTIM_DIM];
	static double deriv[EQ_CONSTR_DIM][OPTIM_DIM];
	static double num_deriv[EQ_CONSTR_DIM][OPTIM_DIM];
	static double jac[2*NCART][NDOF];
	static double jac_w[NCART][NDOF];
	double h = 0.001;

	for (int i = 0; i < NDOF; i++) {
		x[i] = xdiff[i] = q0[i];
		x[i+NDOF] = xdiff[i+NDOF] = randn();
	}

	/*
	 * Calculate exact derivatives
	 */
	calc_racket_state(x,racket_pos,racket_normal,jac);

	printf("racket_normal = [%f,%f,%f]\n",racket_normal[X],racket_normal[Y],racket_normal[Z]);
	// copy geometric jacobians angular velocity subblock to jac_w
	for (int i = 0; i < NDOF; i++) {
		for (int j = 0; j < NCART; j++) {
			jac_w[j][i] = jac[j+NCART][i];
		}
	}
	double dndq[NCART][NDOF];
	printf("jac_w = \n");
	print_mat(jac_w);
	cross_prods(jac_w,racket_normal,dndq);
	printf("dndq = \n");
	print_mat(dndq);

	// copy geometric jacobians angular velocity subblock to jac_w
	for (int i = 0; i < NDOF; i++) {
		for (int j = 0; j < NCART; j++) {
			jac_w[j][i] = jac[j+NCART][i];
		}
	}

	/*
	 * Calculate numerical derivatives
	 */
	for (int i = 0; i < OPTIM_DIM; i++) {
		xdiff[i] = x[i] + h;
		kinematics_eq_constr(constr1,xdiff,(void*)&opt);
		xdiff[i] = x[i] - h;
		kinematics_eq_constr(constr2,xdiff,(void*)&opt);
		for (int j = 0; j < EQ_CONSTR_DIM; j++) {
			num_deriv[j][i] = (constr1[j] - constr2[j]) / (2*h);
		}
		xdiff[i] = x[i];
	}
	printf("dndq_num = \n");
	for (int i = NCART; i < 2*NCART; i++) {
		for (int j = 0; j < NDOF; j++) {
			printf("%.2f\t", num_deriv[i][j]);
		}
		printf("\n");
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
	for (int i = NCART; i < 2*NCART; i++) {
		for (int j = 0; j < NDOF; j++) {
			deriv[i][j] = dndq[i-NCART][j];
		}
	}

	/*
	 * Calculate the maximum difference between numerical and exact derivative matrices (big jacobian)
	 */
	BOOST_TEST(calc_max_diff(deriv,num_deriv,0,2*NCART,0,2*NDOF) < 1e-3);
}


/**
 * @brief Comparing with MATLAB the racket states calculated given q, qd
 *
 * Comparing both C and C++ versions of calc_racket_state function
 * against MATLAB.
 *
 */
void test_kinematics_calculations() {

	BOOST_TEST_MESSAGE("\nComparing racket state calculations with MATLAB...");

	static double q[NDOF], qdot[NDOF], pos[NCART], vel[NCART], normal[NCART];
	for (int i = 0; i < NDOF; i++) {
		q[i] = 1.0;
		qdot[i] = 0.0;
	}
	vec7 q0_cpp = ones<vec>(NDOF);
	joint q_cpp;
	q_cpp.q = q0_cpp;
	player::racket robot_racket;

	// C version
	calc_racket_state(q,qdot,pos,vel,normal);
	// C++ version
	calc_racket_state(q_cpp,robot_racket);

	vec3 x_c(pos);
	vec3 xdot_c(vel);
	vec3 n_c(normal);

	vec3 x_MATLAB = {-0.8309, 0.0861, -0.6672};
	vec3 xdot_MATLAB = zeros<vec>(3);
	vec3 n_MATLAB = {0.3857, 0.0174, -0.9225};

	/*cout << "Kin values:" << endl;
	cout << x << xdot << n << endl;
	cout << "MATLAB values:" << endl;
	cout << x_MATLAB << xdot_MATLAB << n_MATLAB << endl;*/

	BOOST_TEST(approx_equal(x_c,x_MATLAB,"absdiff", 0.002));
	BOOST_TEST(approx_equal(xdot_c,xdot_MATLAB,"absdiff", 0.002));
	BOOST_TEST(approx_equal(n_c,n_MATLAB,"absdiff", 0.002));
	BOOST_TEST(approx_equal(robot_racket.pos,x_c,"absdiff", 0.002));
	BOOST_TEST(approx_equal(robot_racket.vel,xdot_c,"absdiff", 0.002));
	BOOST_TEST(approx_equal(robot_racket.normal,n_c,"absdiff", 0.002));


}

/*
 * This is the constraint that makes sure we hit the ball
 */
static void kinematics_eq_constr(double *result,
                                 const double *x,
                                 void *my_function_data) {

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

static double calc_max_diff(const double mat1[EQ_CONSTR_DIM][OPTIM_DIM],
                            const double mat2[EQ_CONSTR_DIM][OPTIM_DIM],
                            int m1, int m2, int n1, int n2) {

    double maxdiff = 0.0;
    double absdiff = 0.0;
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
 * Find the cross products between columns of mat matrix and the given vector
 */
static void cross_prods(const double mat[NCART][NDOF],
                          const double v[NCART],
                          double out[NCART][NDOF]) {

    for (int i = 0; i < NDOF; i++) {
        out[X][i] = mat[Y][i] * v[Z] - mat[Z][i] * v[Y];
        out[Y][i] = mat[Z][i] * v[X] - mat[X][i] * v[Z];
        out[Z][i] = mat[X][i] * v[Y] - mat[Y][i] * v[X];
    }

}


static void print_mat(const double mat[NCART][NDOF]) {

    for (int i = 0; i < NCART; i++) {
        for (int j = 0; j < NDOF; j++) {
            printf("%.2f\t", mat[i][j]);
        }
        printf("\n");
    }
}

