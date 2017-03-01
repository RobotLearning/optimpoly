/*
 * kinematics.cpp
 *
 * Kinematics functions are stored here
 *
 *  Created on: Feb 12, 2017
 *      Author: okoc
 */

#include <iostream>
#include <armadillo>
#include "constants.h"
#include "kinematics.hpp"
#include "player.hpp"
#include "tabletennis.h"

using namespace arma;

/*
 * Calculates cartesian racket pos, vel and normal
 * given joint positions and velocities
 *
 */
void calc_racket_state(const joint & robot_joint,
		               racket & robot_racket) {

	static mat::fixed<3,7> origin = zeros<mat>(3,7);
	static mat::fixed<3,7> axis = zeros<mat>(3,7);
	static mat::fixed<3,6> link = zeros<mat>(3,6);
	static mat::fixed<6,7> jac = zeros<mat>(6,7);
	static cube::fixed<4,4,7> amats = zeros<cube>(4,4,7);

	kinematics(robot_joint.q,link,origin,axis,amats);
	//rotate_to_quat(amats.slice(PALM)(span(X,Z),span(X,Z)),quat);
	//calc_racket_orient(quat);
	jacobian(link,origin,axis,jac);
	robot_racket.pos = link.col(PALM);
	robot_racket.vel = jac.rows(X,Z) * robot_joint.qd;
	robot_racket.normal = amats.slice(PALM).col(1).head(3);
}

/*
 * Rotate racket by 90 degrees to get
 * racket orientation from endeffector orientation
 *
 */
void calc_racket_orient(vec4 & quat) {

	double pi = datum::pi;
	vec4 rot_quat = {cos(pi/4), -sin(pi/4), 0, 0};
	vec4 quat_new;
	mult_two_quats(rot_quat,quat,quat_new);
	quat = quat_new;
}

/* Function to multiply two quaternions
 *
 */
void mult_two_quats(const vec4 & q1, const vec4 & q2, vec4 & q3) {

	q3(0) = q1(0)*q2(0) - q1(1)*q1(1) - q1(2)*q1(2) - q1(3)*q1(3);
	q3(1) = q1(1)*q2(0) + q1(0)*q2(1) + q1(2)*q2(3) - q1(3)*q2(2);
	q3(2) = q1(0)*q2(2) - q1(1)*q2(3) + q1(2)*q2(0) + q1(3)*q2(1);
	q3(3) = q1(0)*q2(3) + q1(1)*q2(2) - q1(2)*q2(1) + q1(3)*q2(0);
}


/*
 * Form the quaternion out of the rotation matrix
 *
 */
void rotate_to_quat(const mat33 & R, vec4 & quat) {

	double T,S;
	T = 1.0 + R(0,0) + R(1,1) + R(2,2);

	if (T > 0.00000001) {
		S  = 0.5 / sqrt(T);
		quat(0) = 0.25 / S;
		quat(1) = (R(2,1) - R(1,2)) * S;
		quat(2) = (R(0,2) - R(2,0)) * S;
		quat(3) = (R(1,0) - R(0,1)) * S;
	}
	else {
		if ((R(0,0) > R(1,1)) && (R(0,0) > R(2,2))) {
			S = sqrt(1.0 + R(0,0) - R(1,1) - R(2,2) ) * 2;
			quat(0) = 0.25 * S;
			quat(1) = (R(0,1) + R(1,0) ) / S;
			quat(2) = (R(0,2) + R(2,0) ) / S;
			quat(3) = (R(1,2) - R(2,1) ) / S;
		}
		else if (R(1,1) > R(2,2)) {
			S = sqrt(1.0 + R(1,1) - R(0,0) - R(2,2) ) * 2;
			quat(0) = (R(0,1) + R(1,0) ) / S;
			quat(1) = 0.25 * S;
			quat(2) = (R(1,2) + R(2,1) ) / S;
			quat(3) = (R(0,2) - R(2,0) ) / S;
		}
		else {
			S = sqrt(1.0 + R(2,2) - R(0,0) - R(1,1)) * 2;
			quat(0) = (R(0,2) + R(2,0)) / S;
			quat(1) = (R(1,2) + R(2,1)) / S;
			quat(2) = 0.25 * S;
			quat(3) = (R(0,1) - R(1,0)) / S;
		}
	}
}

/*
 * Computes the jacobian. Taken from SL and simplified.
 *
 * Function Parameters: [in]=input,[out]=output
 *
 * \param[in]     lp      : the link positions
 * \param[in]     jop     : joint origin positions
 * \param[in]     jap     : joint axix unit vectors
 * \param[out]    Jac     : the jacobian
 *
 */
void jacobian(const mat & lp, const mat & jop, const mat & jap, mat & jac) {

	vec6 col;
	for (int i = 0; i < NDOF; ++i) {
		revolute_jac_col(lp.col(PALM), jop.col(i), jap.col(i), col);
		jac.col(i) = col;
	}

}

/*
 *
 * Copied from SL_common.
 *
 * computes one column for the geometric jacobian of a revolute joint
 * from the given input vectors
 *
 *
 * p    : position of endeffector
 * pi   : position of joint origin
 * zi   : unit vector of joint axis
 * col  : column vector of Jacobian [out]
 *
 */
void revolute_jac_col(const vec3 & p, const vec3 & pi, const vec3 & zi, vec6 & col) {

	col(span(X,Z)) = cross(zi, p-pi);
	col(span(DX,DZ)) = zi;

}

/* Barrett WAM forward kinematics
 * used to show trajectories in cartesian space
 *
 * Function taken from SL:
 * shared/barrett/math/LInfo_declare.h
 * shared/barrett/math/LInfo_math.h
 *
 * \param[out]    Xaxis   : array of rotation axes (z)
 * \param[out]    Xorigin : array of coord.sys. origin vectors
 * \param[out]    Xlink   : array of link position
 * \param[out]    Amats   : homogeneous transformation matrices of each link
 *
 */
void kinematics(const vec7 & q, mat & Xlink, mat & Xorigin, mat & Xaxis, cube & Amats) {

	 static double  ss0th;
	 static double  cs0th;
	 static double  ss1th;
	 static double  cs1th;
	 static double  ss2th;
	 static double  cs2th;
	 static double  ss3th;
	 static double  cs3th;
	 static double  ss4th;
	 static double  cs4th;
	 static double  ss5th;
	 static double  cs5th;
	 static double  ss6th;
	 static double  cs6th;

	 static double  rseff0a0;
	 static double  rceff0a0;
	 static double  rseff0a1;
	 static double  rceff0a1;
	 static double  rseff0a2;
	 static double  rceff0a2;

	 static mat::fixed<3,4>  Hi00 = zeros<mat>(3,4);
	 static mat::fixed<3,4>  Hi01 = zeros<mat>(3,4);
	 static mat::fixed<3,2>  Hi12 = zeros<mat>(3,2);
	 static mat::fixed<3,4>  Hi23 = zeros<mat>(3,4);
	 static mat::fixed<3,4>  Hi34 = zeros<mat>(3,4);
	 static mat::fixed<3,4>  Hi45 = zeros<mat>(3,4);
	 static mat::fixed<3,4>  Hi56 = zeros<mat>(3,4);
	 static mat::fixed<3,2>  Hi67 = zeros<mat>(3,2);
	 static mat::fixed<3,4>  Hi78 = zeros<mat>(3,4);

	 static mat::fixed<3,4>  Ai01 = zeros<mat>(3,4);
	 static mat::fixed<3,4>  Ai02 = zeros<mat>(3,4);
	 static mat::fixed<3,4>  Ai03 = zeros<mat>(3,4);
	 static mat::fixed<3,4>  Ai04 = zeros<mat>(3,4);
	 static mat::fixed<3,4>  Ai05 = zeros<mat>(3,4);
	 static mat::fixed<3,4>  Ai06 = zeros<mat>(3,4);
	 static mat::fixed<3,4>  Ai07 = zeros<mat>(3,4);
	 static mat::fixed<3,4>  Ai08 = zeros<mat>(3,4);

	 vec3 pos = {0.0, 0.0, 0.30};
	 vec3 orient = zeros<vec>(3);
	 eff racket = {pos, orient};
	 vec3 basec = zeros<vec>(3);
	 vec4 baseo = {0.0, 1.0, 0.0, 0.0};
	 pose base = {basec, baseo};


	// sine and cosine precomputation
	ss0th = sin(q(0));
	cs0th = cos(q(0));
	ss1th = sin(q(1));
	cs1th = cos(q(1));
	ss2th = sin(q(2));
	cs2th = cos(q(2));
	ss3th = sin(q(3));
	cs3th = cos(q(3));
	ss4th = sin(q(4));
	cs4th = cos(q(4));
	ss5th = sin(q(5));
	cs5th = cos(q(5));
	ss6th = sin(q(6));
	cs6th = cos(q(6));

	// endeffector orientations


	rseff0a0 = sin(racket.o(0));
	rceff0a0 = cos(racket.o(0));
	rseff0a1 = sin(racket.o(1));
	rceff0a1 = cos(racket.o(1));
	rseff0a2 = sin(racket.o(2));
	rceff0a2 = cos(racket.o(2));

	// Calculations are done here

	// inverse homogeneous rotation matrices
	Hi00(0,0) = -1 + 2*pow(base.q(0),2) + 2*pow(base.q(1),2);
	Hi00(0,1) = 2*(base.q(1)*base.q(2) - base.q(0)*base.q(3));
	Hi00(0,2) = 2*(base.q(0)*base.q(2) + base.q(1)*base.q(3));
	Hi00(0,3) = base.x(0);
	Hi00(1,0) = 2*(base.q(1)*base.q(2) + base.q(0)*base.q(3));
	Hi00(1,1) = -1 + 2*pow(base.q(0),2) + 2*pow(base.q(2),2);
	Hi00(1,2) = 2*(-(base.q(0)*base.q(1)) + base.q(2)*base.q(3));
	Hi00(1,3) = base.x(1);
	Hi00(2,0) = 2*(-(base.q(0)*base.q(2)) + base.q(1)*base.q(3));
	Hi00(2,1) = 2*(base.q(0)*base.q(1) + base.q(2)*base.q(3));
	Hi00(2,2) = -1 + 2*pow(base.q(0),2) + 2*pow(base.q(3),2);
	Hi00(2,3) = base.x(2);

	Hi01(0,0) = cs0th;
	Hi01(0,1) = -ss0th;
	Hi01(1,0) = ss0th;
	Hi01(1,1) = cs0th;
	Hi01(2,3) = ZSFE;
	Hi12(1,0) = ss1th;
	Hi12(1,1) = cs1th;
	Hi12(2,0) = cs1th;
	Hi12(2,1) = -ss1th;
	Hi23(0,3) = ZHR;
	Hi23(1,0) = ss2th;
	Hi23(1,1) = cs2th;
	Hi23(2,0) = -cs2th;
	Hi23(2,1) = ss2th;
	Hi34(1,0) = ss3th;
	Hi34(1,1) = cs3th;
	Hi34(1,3) = YEB;
	Hi34(2,0) = cs3th;
	Hi34(2,1) = -ss3th;
	Hi34(2,3) = ZEB;
	Hi45(0,3) = ZWR;
	Hi45(1,0) = ss4th;
	Hi45(1,1) = cs4th;
	Hi45(1,3) = YWR;
	Hi45(2,0) = -cs4th;
	Hi45(2,1) = ss4th;
	Hi56(1,0) = ss5th;
	Hi56(1,1) = cs5th;
	Hi56(2,0) = cs5th;
	Hi56(2,1) = -ss5th;
	Hi56(2,3) = ZWFE;
	Hi67(1,0) = ss6th;
	Hi67(1,1) = cs6th;
	Hi67(2,0) = -cs6th;
	Hi67(2,1) = ss6th;
	Hi78(0,0) = rceff0a1*rceff0a2;
	Hi78(0,1) = -(rceff0a1*rseff0a2);
	Hi78(0,2) = rseff0a1;
	Hi78(0,3) = racket.x(0);
	Hi78(1,0) = rceff0a2*rseff0a0*rseff0a1 + rceff0a0*rseff0a2;
	Hi78(1,1) = rceff0a0*rceff0a2 - rseff0a0*rseff0a1*rseff0a2;
	Hi78(1,2) = -(rceff0a1*rseff0a0);
	Hi78(1,3) = racket.x(1);
	Hi78(2,0) = -(rceff0a0*rceff0a2*rseff0a1) + rseff0a0*rseff0a2;
	Hi78(2,1) = rceff0a2*rseff0a0 + rceff0a0*rseff0a1*rseff0a2;
	Hi78(2,2) = rceff0a0*rceff0a1;
	Hi78(2,3) = racket.x(2);

	/*cout << "Hi00" << endl << Hi00 << endl;
	cout << "Hi01" << endl << Hi01 << endl;
	cout << "Hi12" << endl << Hi12 << endl;
	cout << "Hi23" << endl << Hi23 << endl;
	cout << "Hi34" << endl << Hi34 << endl;
	cout << "Hi45" << endl << Hi45 << endl;
	cout << "Hi56" << endl << Hi56 << endl;
	cout << "Hi67" << endl << Hi67 << endl;
	cout << "Hi78" << endl << Hi78 << endl;*/

	// per link inverse homogeneous rotation matrices
	Ai01(0,0) = Hi00(0,0)*Hi01(0,0) + Hi00(0,1)*Hi01(1,0);
	Ai01(0,1) = Hi00(0,0)*Hi01(0,1) + Hi00(0,1)*Hi01(1,1);
	Ai01(0,2) = Hi00(0,2);
	Ai01(0,3) = Hi00(0,3) + Hi00(0,2)*Hi01(2,3);
	Ai01(1,0) = Hi00(1,0)*Hi01(0,0) + Hi00(1,1)*Hi01(1,0);
	Ai01(1,1) = Hi00(1,0)*Hi01(0,1) + Hi00(1,1)*Hi01(1,1);
	Ai01(1,2) = Hi00(1,2);
	Ai01(1,3) = Hi00(1,3) + Hi00(1,2)*Hi01(2,3);
	Ai01(2,0) = Hi00(2,0)*Hi01(0,0) + Hi00(2,1)*Hi01(1,0);
	Ai01(2,1) = Hi00(2,0)*Hi01(0,1) + Hi00(2,1)*Hi01(1,1);
	Ai01(2,2) = Hi00(2,2);
	Ai01(2,3) = Hi00(2,3) + Hi00(2,2)*Hi01(2,3);
	Ai02(0,0) = Ai01(0,1)*Hi12(1,0) + Ai01(0,2)*Hi12(2,0);
	Ai02(0,1) = Ai01(0,1)*Hi12(1,1) + Ai01(0,2)*Hi12(2,1);
	Ai02(0,2) = -Ai01(0,0);
	Ai02(0,3) = Ai01(0,3);
	Ai02(1,0) = Ai01(1,1)*Hi12(1,0) + Ai01(1,2)*Hi12(2,0);
	Ai02(1,1) = Ai01(1,1)*Hi12(1,1) + Ai01(1,2)*Hi12(2,1);
	Ai02(1,2) = -Ai01(1,0);
	Ai02(1,3) = Ai01(1,3);
	Ai02(2,0) = Ai01(2,1)*Hi12(1,0) + Ai01(2,2)*Hi12(2,0);
	Ai02(2,1) = Ai01(2,1)*Hi12(1,1) + Ai01(2,2)*Hi12(2,1);
	Ai02(2,2) = -Ai01(2,0);
	Ai02(2,3) = Ai01(2,3);
	Ai03(0,0) = Ai02(0,1)*Hi23(1,0) + Ai02(0,2)*Hi23(2,0);
	Ai03(0,1) = Ai02(0,1)*Hi23(1,1) + Ai02(0,2)*Hi23(2,1);
	Ai03(0,2) = Ai02(0,0);
	Ai03(0,3) = Ai02(0,3) + Ai02(0,0)*Hi23(0,3);
	Ai03(1,0) = Ai02(1,1)*Hi23(1,0) + Ai02(1,2)*Hi23(2,0);
	Ai03(1,1) = Ai02(1,1)*Hi23(1,1) + Ai02(1,2)*Hi23(2,1);
	Ai03(1,2) = Ai02(1,0);
	Ai03(1,3) = Ai02(1,3) + Ai02(1,0)*Hi23(0,3);
	Ai03(2,0) = Ai02(2,1)*Hi23(1,0) + Ai02(2,2)*Hi23(2,0);
	Ai03(2,1) = Ai02(2,1)*Hi23(1,1) + Ai02(2,2)*Hi23(2,1);
	Ai03(2,2) = Ai02(2,0);
	Ai03(2,3) = Ai02(2,3) + Ai02(2,0)*Hi23(0,3);
	Ai04(0,0) = Ai03(0,1)*Hi34(1,0) + Ai03(0,2)*Hi34(2,0);
	Ai04(0,1) = Ai03(0,1)*Hi34(1,1) + Ai03(0,2)*Hi34(2,1);
	Ai04(0,2) = -Ai03(0,0);
	Ai04(0,3) = Ai03(0,3) + Ai03(0,1)*Hi34(1,3) + Ai03(0,2)*Hi34(2,3);
	Ai04(1,0) = Ai03(1,1)*Hi34(1,0) + Ai03(1,2)*Hi34(2,0);
	Ai04(1,1) = Ai03(1,1)*Hi34(1,1) + Ai03(1,2)*Hi34(2,1);
	Ai04(1,2) = -Ai03(1,0);
	Ai04(1,3) = Ai03(1,3) + Ai03(1,1)*Hi34(1,3) + Ai03(1,2)*Hi34(2,3);
	Ai04(2,0) = Ai03(2,1)*Hi34(1,0) + Ai03(2,2)*Hi34(2,0);
	Ai04(2,1) = Ai03(2,1)*Hi34(1,1) + Ai03(2,2)*Hi34(2,1);
	Ai04(2,2) = -Ai03(2,0);
	Ai04(2,3) = Ai03(2,3) + Ai03(2,1)*Hi34(1,3) + Ai03(2,2)*Hi34(2,3);
	Ai05(0,0) = Ai04(0,1)*Hi45(1,0) + Ai04(0,2)*Hi45(2,0);
	Ai05(0,1) = Ai04(0,1)*Hi45(1,1) + Ai04(0,2)*Hi45(2,1);
	Ai05(0,2) = Ai04(0,0);
	Ai05(0,3) = Ai04(0,3) + Ai04(0,0)*Hi45(0,3) + Ai04(0,1)*Hi45(1,3);
	Ai05(1,0) = Ai04(1,1)*Hi45(1,0) + Ai04(1,2)*Hi45(2,0);
	Ai05(1,1) = Ai04(1,1)*Hi45(1,1) + Ai04(1,2)*Hi45(2,1);
	Ai05(1,2) = Ai04(1,0);
	Ai05(1,3) = Ai04(1,3) + Ai04(1,0)*Hi45(0,3) + Ai04(1,1)*Hi45(1,3);
	Ai05(2,0) = Ai04(2,1)*Hi45(1,0) + Ai04(2,2)*Hi45(2,0);
	Ai05(2,1) = Ai04(2,1)*Hi45(1,1) + Ai04(2,2)*Hi45(2,1);
	Ai05(2,2) = Ai04(2,0);
	Ai05(2,3) = Ai04(2,3) + Ai04(2,0)*Hi45(0,3) + Ai04(2,1)*Hi45(1,3);
	Ai06(0,0) = Ai05(0,1)*Hi56(1,0) + Ai05(0,2)*Hi56(2,0);
	Ai06(0,1) = Ai05(0,1)*Hi56(1,1) + Ai05(0,2)*Hi56(2,1);
	Ai06(0,2) = -Ai05(0,0);
	Ai06(0,3) = Ai05(0,3) + Ai05(0,2)*Hi56(2,3);
	Ai06(1,0) = Ai05(1,1)*Hi56(1,0) + Ai05(1,2)*Hi56(2,0);
	Ai06(1,1) = Ai05(1,1)*Hi56(1,1) + Ai05(1,2)*Hi56(2,1);
	Ai06(1,2) = -Ai05(1,0);
	Ai06(1,3) = Ai05(1,3) + Ai05(1,2)*Hi56(2,3);
	Ai06(2,0) = Ai05(2,1)*Hi56(1,0) + Ai05(2,2)*Hi56(2,0);
	Ai06(2,1) = Ai05(2,1)*Hi56(1,1) + Ai05(2,2)*Hi56(2,1);
	Ai06(2,2) = -Ai05(2,0);
	Ai06(2,3) = Ai05(2,3) + Ai05(2,2)*Hi56(2,3);
	Ai07(0,0) = Ai06(0,1)*Hi67(1,0) + Ai06(0,2)*Hi67(2,0);
	Ai07(0,1) = Ai06(0,1)*Hi67(1,1) + Ai06(0,2)*Hi67(2,1);
	Ai07(0,2) = Ai06(0,0);
	Ai07(0,3) = Ai06(0,3);
	Ai07(1,0) = Ai06(1,1)*Hi67(1,0) + Ai06(1,2)*Hi67(2,0);
	Ai07(1,1) = Ai06(1,1)*Hi67(1,1) + Ai06(1,2)*Hi67(2,1);
	Ai07(1,2) = Ai06(1,0);
	Ai07(1,3) = Ai06(1,3);
	Ai07(2,0) = Ai06(2,1)*Hi67(1,0) + Ai06(2,2)*Hi67(2,0);
	Ai07(2,1) = Ai06(2,1)*Hi67(1,1) + Ai06(2,2)*Hi67(2,1);
	Ai07(2,2) = Ai06(2,0);
	Ai07(2,3) = Ai06(2,3);
	Ai08(0,0) = Ai07(0,0)*Hi78(0,0) + Ai07(0,1)*Hi78(1,0) + Ai07(0,2)*Hi78(2,0);
	Ai08(0,1) = Ai07(0,0)*Hi78(0,1) + Ai07(0,1)*Hi78(1,1) + Ai07(0,2)*Hi78(2,1);
	Ai08(0,2) = Ai07(0,0)*Hi78(0,2) + Ai07(0,1)*Hi78(1,2) + Ai07(0,2)*Hi78(2,2);
	Ai08(0,3) = Ai07(0,3) + Ai07(0,0)*Hi78(0,3) + Ai07(0,1)*Hi78(1,3) + Ai07(0,2)*Hi78(2,3);
	Ai08(1,0) = Ai07(1,0)*Hi78(0,0) + Ai07(1,1)*Hi78(1,0) + Ai07(1,2)*Hi78(2,0);
	Ai08(1,1) = Ai07(1,0)*Hi78(0,1) + Ai07(1,1)*Hi78(1,1) + Ai07(1,2)*Hi78(2,1);
	Ai08(1,2) = Ai07(1,0)*Hi78(0,2) + Ai07(1,1)*Hi78(1,2) + Ai07(1,2)*Hi78(2,2);
	Ai08(1,3) = Ai07(1,3) + Ai07(1,0)*Hi78(0,3) + Ai07(1,1)*Hi78(1,3) + Ai07(1,2)*Hi78(2,3);
	Ai08(2,0) = Ai07(2,0)*Hi78(0,0) + Ai07(2,1)*Hi78(1,0) + Ai07(2,2)*Hi78(2,0);
	Ai08(2,1) = Ai07(2,0)*Hi78(0,1) + Ai07(2,1)*Hi78(1,1) + Ai07(2,2)*Hi78(2,1);
	Ai08(2,2) = Ai07(2,0)*Hi78(0,2) + Ai07(2,1)*Hi78(1,2) + Ai07(2,2)*Hi78(2,2);
	Ai08(2,3) = Ai07(2,3) + Ai07(2,0)*Hi78(0,3) + Ai07(2,1)*Hi78(1,3) + Ai07(2,2)*Hi78(2,3);

	/*cout << "Ai01" << endl << Ai01 << endl;
	cout << "Ai02" << endl << Ai02 << endl;
	cout << "Ai03" << endl << Ai03 << endl;
	cout << "Ai04" << endl << Ai04 << endl;
	cout << "Ai05" << endl << Ai05 << endl;
	cout << "Ai06" << endl << Ai06 << endl;
	cout << "Ai07" << endl << Ai07 << endl;
	cout << "Ai08" << endl << Ai08 << endl;*/


	// joint ID: 1
	Xorigin.col(0) = Ai01.col(3);
	Xaxis.col(0)   = Ai01.col(2);
	// joint ID: 2
	Xorigin.col(1) = Ai02.col(3);
	Xaxis.col(1) = Ai02.col(2);
	// joint ID: 3
	Xorigin.col(2) = Ai03.col(3);
	Xaxis.col(2) = Ai03.col(2);
	// joint ID: 4
	Xorigin.col(3) = Ai04.col(3);
	Xaxis.col(3) = Ai04.col(2);
	// joint ID: 5
	Xorigin.col(4) = Ai05.col(3);
	Xaxis.col(4) = Ai05.col(2);
	// joint ID: 6
	Xorigin.col(5) = Ai06.col(3);
	Xaxis.col(5) = Ai06.col(2);
	// joint ID: 7
	Xorigin.col(6) = Ai07.col(3);
	Xaxis.col(6) = Ai07.col(2);

	// NOTES:
	Xlink.col(0) = Ai01.col(3);
	Xlink.col(1) = Ai03.col(3);
	Xlink.col(2) = Ai04.col(3);
	Xlink.col(3) = Ai05.col(3);
	Xlink.col(4) = Ai07.col(3);
	Xlink.col(5) = Ai08.col(3);

	Amats.slice(0) = join_vert(Hi00,zeros<rowvec>(4));
	Amats(3,3,0) = 1.0;
	Amats.slice(1) = join_vert(Ai02,zeros<rowvec>(4));
	Amats(3,3,1) = 1.0;
	Amats.slice(2) = join_vert(Ai03,zeros<rowvec>(4));
	Amats(3,3,2) = 1.0;
	Amats.slice(3) = join_vert(Ai04,zeros<rowvec>(4));
	Amats(3,3,3) = 1.0;
	Amats.slice(4) = join_vert(Ai05,zeros<rowvec>(4));
	Amats(3,3,4) = 1.0;
	Amats.slice(5) = join_vert(Ai07,zeros<rowvec>(4));
	Amats(3,3,5) = 1.0;
	Amats.slice(6) = join_vert(Ai08,zeros<rowvec>(4));
	Amats(3,3,6) = 1.0;
}



