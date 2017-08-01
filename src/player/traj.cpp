/*
 * traj.cpp
 *
 * This file includes trajectory generation functions
 * used by player class (player.cpp)
 *
 *  Created on: Jul 29, 2017
 *      Author: okoc
 */

#include <armadillo>
#include "constants.h"
#include "kinematics.h"
#include "player.hpp"
using namespace arma;


/**
 * @brief Generate BATCH 3rd order strike + return polynomials.
 *
 * Based on hitting and returning joint state parameters qf,qfdot
 * and hitting time T, calculates the relevant polynomial parameters
 * and generates BATCH polynomial values till time T.
 * @param qf Hitting joint pos.
 * @param qfdot Hitting joint vels.
 * @param T Hitting time
 * @param qact From actual joint state generate the joint des values
 * @param q_rest_des Desired resting posture
 * @param time2return Time to return after hit to rest posture
 * @param Q Generated joint pos.
 * @param Qd Generated joint vel.
 * @param Qdd Generated joint acc.
 */
void generate_strike(const vec7 & qf, const vec7 & qfdot, const double T,
		             const joint & qact,
		             const vec7 & q_rest_des,
					 const double time2return,
		             mat & Q, mat & Qd, mat & Qdd) {

	// first create hitting polynomials
	vec7 qnow = qact.q;
	vec7 qdnow = qact.qd;
	vec7 a3 = 2.0 * (qnow - qf) / pow(T,3) + (qfdot + qdnow) / pow(T,2);
	vec7 a2 = 3.0 * (qf - qnow) / pow(T,2) - (qfdot + 2.0*qdnow) / T;
	vec7 b3 = 2.0 * (qf - q_rest_des) / pow(time2return,3) + (qfdot) / pow(time2return,2);
	vec7 b2 = 3.0 * (q_rest_des - qf) / pow(time2return,2) - (2.0*qfdot) / time2return;

	int N_hit = T/DT;
	rowvec times_hit = linspace<rowvec>(DT,T,N_hit);
	int N_return = time2return/DT;
	rowvec times_ret = linspace<rowvec>(DT,time2return,N_return);

	mat Q_hit, Qd_hit, Qdd_hit, Q_ret, Qd_ret, Qdd_ret;
	Q_hit = Qd_hit = Qdd_hit = zeros<mat>(NDOF,N_hit);
	Q_ret = Qd_ret = Qdd_ret = zeros<mat>(NDOF,N_return);

	gen_3rd_poly(times_hit,a3,a2,qdnow,qnow,Q_hit,Qd_hit,Qdd_hit);
	gen_3rd_poly(times_ret,b3,b2,qfdot,qf,Q_ret,Qd_ret,Qdd_ret);
	Q = join_horiz(Q_hit,Q_ret);
	Qd = join_horiz(Qd_hit,Qd_ret);
	Qdd = join_horiz(Qdd_hit,Qdd_ret);
}


/**
 * @brief Generate strike and return traj. incrementally
 *
 * Given polynomial parameters saved in poly,
 * move on to the NEXT desired state only (joint pos,vel,acc).
 * @param poly Polynomial parameters updated in OPTIM classes
 * @param q_rest_des FIXED desired resting posture
 * @param time2return FIXED time to return to rest posture after hit
 * @param t The time passed already following trajectory
 * @param qdes Update pos,vel,acc values of this desired joints structure
 * @return
 */
bool update_next_state(const spline_params & poly,
		           const vec7 & q_rest_des,
				   const double time2return,
				   double & t,
				   joint & qdes) {
	mat a,b;
	double tbar;
	bool flag = true;

	if (t <= poly.time2hit) {
		a = poly.a;
		qdes.q = a.col(0)*t*t*t + a.col(1)*t*t + a.col(2)*t + a.col(3);
		qdes.qd = 3*a.col(0)*t*t + 2*a.col(1)*t + a.col(2);
		qdes.qdd = 6*a.col(0)*t + 2*a.col(1);
		t += DT;
		//cout << qdes.q << qdes.qd << qdes.qdd << endl;
	}
	else if (t <= poly.time2hit + time2return) {
		b = poly.b;
		tbar = t - poly.time2hit;
		qdes.q = b.col(0)*tbar*tbar*tbar + b.col(1)*tbar*tbar + b.col(2)*tbar + b.col(3);
		qdes.qd = 3*b.col(0)*tbar*tbar + 2*b.col(1)*tbar + b.col(2);
		qdes.qdd = 6*b.col(0)*tbar + 2*b.col(1);
		t += DT;
	}
	else {
		//printf("Hitting finished!\n");
		t = 0.0;
		flag = false;
		qdes.q = q_rest_des;
		qdes.qd = zeros<vec>(NDOF);
		qdes.qdd = zeros<vec>(NDOF);
	}
	return flag;
}

/*
 * Generate matrix of joint angles, velocities and accelerations
 */
void gen_3rd_poly(const rowvec & times, const vec7 & a3, const vec7 & a2, const vec7 & a1, const vec7 & a0,
		     mat & Q, mat & Qd, mat & Qdd) {

	// IN MATLAB:
	//	qStrike(m,:) = a(1)*t.^3 + a(2)*t.^2 + a(3)*t + a(4);
	//	qdStrike(m,:) = 3*a(1)*t.^2 + 2*a(2)*t + a(3);
	//	qddStrike(m,:) = 6*a(1)*t + 2*a(2);

	for(int i = 0; i < NDOF; i++) {
		Q.row(i) = a3(i) * pow(times,3) + a2(i) * pow(times,2) + a1(i) * times + a0(i);
		Qd.row(i) = 3*a3(i) * pow(times,2) + 2*a2(i) * times + a1(i);
		Qdd.row(i) = 6*a3(i) * times + 2*a2(i);
	}
}

/*
 * Set upper and lower bounds on the optimization.
 * First loads the joint limits and then puts some slack
 */
void set_bounds(double *lb, double *ub, double SLACK, double Tmax) {

	read_joint_limits(lb,ub);
	// lower bounds and upper bounds for qf are the joint limits
	for (int i = 0; i < NDOF; i++) {
		ub[i] -= SLACK;
		lb[i] += SLACK;
		ub[i+NDOF] = MAX_VEL;
		lb[i+NDOF] = -MAX_VEL;
	}
	// constraints on final time
	ub[2*NDOF] = Tmax;
	lb[2*NDOF] = 0.01;
}



