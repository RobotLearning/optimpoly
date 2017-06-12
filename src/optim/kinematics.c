/**
 * @file kinematics.c
 *
 * @brief Here we include the kinematics related functions taken from SL
 *
 * Optimizators call these functions to calculate the Cartesian
 * racket positions, velocities and normals
 *
 *  Created on: Jun 22, 2016
 *      Author: okoc
 */

#include "stdio.h"
#include "stdlib.h"
#include "constants.h"
#include "utils.h"
#include "string.h"
#include "kinematics.h"

// internal functions used in calculating racket quantities
static void get_cart_velocity(double jac[2*NCART][NDOF],
		               const double qdot[NDOF],
					   double vel[NCART]);
static void kinematics(const double state[NDOF],
		        double Xlink[NLINK+1][4],
				double Xorigin[NDOF+1][4],
				double Xaxis[NDOF+1][4],
		        double Ahmat[NDOF+1][5][5]);
static void jacobian(const double link[NLINK+1][4],
		const double origin[NDOF+1][4],
		const double axis[NDOF+1][4],
		double jac[2*NCART][NDOF]);
static void set_endeffector(double endeff_pos[NCART]);

/**
 * @brief Calculates cartesian racket pos, vel and normal
 * given joint positions and velocities.
 *
 * Calls the kinematics and the (geometric) jacobian functions
 * to compute the racket quantities.
 *
 * @param q Joint positions.
 * @param qdot Joint velocities.
 * @param pos Racket positions.
 * @param vel Racket velocities.
 * @param normal Racket normals.
 */
void calc_racket_state(const double q[NDOF],
		               const double qdot[NDOF],
					   double pos[NCART],
					   double vel[NCART],
					   double normal[NCART]) {

	static const int PALM = 6;

	static double link[NLINK+1][3+1];
	static double origin[NDOF+1][3+1];
	static double axis[NDOF+1][3+1];
	static double amats[NDOF+1][4+1][4+1];
	static double jacobi[2*NCART][NDOF];

	kinematics(q,link,origin,axis,amats);
	//rotate_to_quat(amats.slice(PALM)(span(X,Z),span(X,Z)),quat);
	//calc_racket_orient(quat);
	jacobian(link,origin,axis,jacobi);

	for (int i = 0; i < NCART; i++) {
		pos[i] = link[PALM][i+1];
		normal[i] = amats[PALM][i+1][2];
	}
	get_cart_velocity(jacobi,qdot,vel);
}

/**
 * @brief Return racket positions, normal and the jacobian
 *
 * Makes the input pos vector equal to racket positions for given joints q
 * Makes the input matrix equal to the jacobian at input q
 *
 * Useful to test derivatives of kinematics
 */
void get_racket_state(const double q[NDOF], double pos[NCART], double normal[NCART], double jacobi[2*NCART][NDOF]) {

	static const int PALM = 6;
	static double link[NLINK+1][3+1];
	static double origin[NDOF+1][3+1];
	static double axis[NDOF+1][3+1];
	static double amats[NDOF+1][4+1][4+1];
	kinematics(q,link,origin,axis,amats);
	jacobian(link,origin,axis,jacobi);
	for (int i = 0; i < NCART; i++) {
		pos[i] = link[PALM][i+1];
		normal[i] = amats[PALM][i+1][2];
	}
}

/**
 * @brief Returns the cartesian racket positions
 */
void get_position(double q[NDOF]) {
	static double link[NLINK+1][3+1];
	static double origin[NDOF+1][3+1];
	static double axis[NDOF+1][3+1];
	static double amats[NDOF+1][4+1][4+1];
	kinematics(q,link,origin,axis,amats);
}

/*
 * Find cartesian racket velocity given
 * joint velocities and Jacobian
 */
static void get_cart_velocity(double jac[2*NCART][NDOF],
		               const double qdot[NDOF],
					   double vel[NCART]) {

	for (int i = 0; i < NCART; i++) {
		vel[i] = 0.0;
		for (int j = 0; j < NDOF; j++) {
			vel[i] += jac[i][j] * qdot[j];
		}
	}
}

/*
 *
 * Kinematics from SL
 * TODO: Get rid of 1-based indexing for the 2-5th arguments!
 *
 */
static void kinematics(const double state[NDOF],
		        double Xlink[NLINK+1][4],
				double Xorigin[NDOF+1][4],
				double Xaxis[NDOF+1][4],
		        double Ahmat[NDOF+1][5][5]) {

	static int firsttime = TRUE;
	static double basec[3+1] = {0.0};
	static double baseo[4+1] = {0.0};
	static double eff_a[NCART+1];
	static double eff_x[NCART+1];

	if (firsttime) {
		firsttime = FALSE;
		eff_x[3] = 0.3; // attach the racket
		baseo[2] = 1.0;
	}

	static double  sstate1th;
	static double  cstate1th;
	static double  sstate2th;
	static double  cstate2th;
	static double  sstate3th;
	static double  cstate3th;
	static double  sstate4th;
	static double  cstate4th;
	static double  sstate5th;
	static double  cstate5th;
	static double  sstate6th;
	static double  cstate6th;
	static double  sstate7th;
	static double  cstate7th;

	static double  rseff1a1;
	static double  rceff1a1;
	static double  rseff1a2;
	static double  rceff1a2;
	static double  rseff1a3;
	static double  rceff1a3;

	static double  Hi00[4+1][4+1];
	static double  Hi01[4+1][4+1];
	static double  Hi12[4+1][4+1];
	static double  Hi23[4+1][4+1];
	static double  Hi34[4+1][4+1];
	static double  Hi45[4+1][4+1];
	static double  Hi56[4+1][4+1];
	static double  Hi67[4+1][4+1];
	static double  Hi78[4+1][4+1];

	static double  Ai01[4+1][4+1];
	static double  Ai02[4+1][4+1];
	static double  Ai03[4+1][4+1];
	static double  Ai04[4+1][4+1];
	static double  Ai05[4+1][4+1];
	static double  Ai06[4+1][4+1];
	static double  Ai07[4+1][4+1];
	static double  Ai08[4+1][4+1];

	/* Need [n_joints+1]x[3+1] matrices: Xorigin,Xmcog,Xaxis, and Xlink[nLinks+1][3+1] */

	/* sine and cosine precomputation */
	sstate1th=Sin(state[0]);
	cstate1th=Cos(state[0]);

	sstate2th=Sin(state[1]);
	cstate2th=Cos(state[1]);

	sstate3th=Sin(state[2]);
	cstate3th=Cos(state[2]);

	sstate4th=Sin(state[3]);
	cstate4th=Cos(state[3]);

	sstate5th=Sin(state[4]);
	cstate5th=Cos(state[4]);

	sstate6th=Sin(state[5]);
	cstate6th=Cos(state[5]);

	sstate7th=Sin(state[6]);
	cstate7th=Cos(state[6]);

	/* rotation matrix sine and cosine precomputation */


	rseff1a1 = Sin(eff_a[1]);
	rceff1a1 = Cos(eff_a[1]);

	rseff1a2 = Sin(eff_a[2]);
	rceff1a2 = Cos(eff_a[2]);

	rseff1a3 = Sin(eff_a[3]);
	rceff1a3 = Cos(eff_a[3]);



	/* inverse homogeneous rotation matrices */
	Hi00[1][1]=-1 + 2*Power(baseo[1],2) + 2*Power(baseo[2],2);
	Hi00[1][2]=2*(baseo[2]*baseo[3] - baseo[1]*baseo[4]);
	Hi00[1][3]=2*(baseo[1]*baseo[3] + baseo[2]*baseo[4]);
	Hi00[1][4]=basec[1];

	Hi00[2][1]=2*(baseo[2]*baseo[3] + baseo[1]*baseo[4]);
	Hi00[2][2]=-1 + 2*Power(baseo[1],2) + 2*Power(baseo[3],2);
	Hi00[2][3]=2*(-(baseo[1]*baseo[2]) + baseo[3]*baseo[4]);
	Hi00[2][4]=basec[2];

	Hi00[3][1]=2*(-(baseo[1]*baseo[3]) + baseo[2]*baseo[4]);
	Hi00[3][2]=2*(baseo[1]*baseo[2] + baseo[3]*baseo[4]);
	Hi00[3][3]=-1 + 2*Power(baseo[1],2) + 2*Power(baseo[4],2);
	Hi00[3][4]=basec[3];


	Hi01[1][1]=cstate1th;
	Hi01[1][2]=-sstate1th;

	Hi01[2][1]=sstate1th;
	Hi01[2][2]=cstate1th;

	Hi01[3][4]=ZSFE;


	Hi12[2][1]=sstate2th;
	Hi12[2][2]=cstate2th;

	Hi12[3][1]=cstate2th;
	Hi12[3][2]=-sstate2th;


	Hi23[1][4]=ZHR;

	Hi23[2][1]=sstate3th;
	Hi23[2][2]=cstate3th;

	Hi23[3][1]=-cstate3th;
	Hi23[3][2]=sstate3th;


	Hi34[2][1]=sstate4th;
	Hi34[2][2]=cstate4th;
	Hi34[2][4]=YEB;

	Hi34[3][1]=cstate4th;
	Hi34[3][2]=-sstate4th;
	Hi34[3][4]=ZEB;


	Hi45[1][4]=ZWR;

	Hi45[2][1]=sstate5th;
	Hi45[2][2]=cstate5th;
	Hi45[2][4]=YWR;

	Hi45[3][1]=-cstate5th;
	Hi45[3][2]=sstate5th;


	Hi56[2][1]=sstate6th;
	Hi56[2][2]=cstate6th;

	Hi56[3][1]=cstate6th;
	Hi56[3][2]=-sstate6th;
	Hi56[3][4]=ZWFE;


	Hi67[2][1]=sstate7th;
	Hi67[2][2]=cstate7th;

	Hi67[3][1]=-cstate7th;
	Hi67[3][2]=sstate7th;


	Hi78[1][1]=rceff1a2*rceff1a3;
	Hi78[1][2]=-(rceff1a2*rseff1a3);
	Hi78[1][3]=rseff1a2;
	Hi78[1][4]=eff_x[1];

	Hi78[2][1]=rceff1a3*rseff1a1*rseff1a2 + rceff1a1*rseff1a3;
	Hi78[2][2]=rceff1a1*rceff1a3 - rseff1a1*rseff1a2*rseff1a3;
	Hi78[2][3]=-(rceff1a2*rseff1a1);
	Hi78[2][4]=eff_x[2];

	Hi78[3][1]=-(rceff1a1*rceff1a3*rseff1a2) + rseff1a1*rseff1a3;
	Hi78[3][2]=rceff1a3*rseff1a1 + rceff1a1*rseff1a2*rseff1a3;
	Hi78[3][3]=rceff1a1*rceff1a2;
	Hi78[3][4]=eff_x[3];


	/*print_mat("Hi00", Hi00);
	print_mat("Hi01", Hi01);
	print_mat("Hi12", Hi12);
	print_mat("Hi23", Hi23);
	print_mat("Hi34", Hi34);
	print_mat("Hi45", Hi45);
	print_mat("Hi56", Hi56);
	print_mat("Hi67", Hi67);
	print_mat("Hi78", Hi78);*/


	/* per link inverse homogeneous rotation matrices */
	Ai01[1][1]=Hi00[1][1]*Hi01[1][1] + Hi00[1][2]*Hi01[2][1];
	Ai01[1][2]=Hi00[1][1]*Hi01[1][2] + Hi00[1][2]*Hi01[2][2];
	Ai01[1][3]=Hi00[1][3];
	Ai01[1][4]=Hi00[1][4] + Hi00[1][3]*Hi01[3][4];

	Ai01[2][1]=Hi00[2][1]*Hi01[1][1] + Hi00[2][2]*Hi01[2][1];
	Ai01[2][2]=Hi00[2][1]*Hi01[1][2] + Hi00[2][2]*Hi01[2][2];
	Ai01[2][3]=Hi00[2][3];
	Ai01[2][4]=Hi00[2][4] + Hi00[2][3]*Hi01[3][4];

	Ai01[3][1]=Hi00[3][1]*Hi01[1][1] + Hi00[3][2]*Hi01[2][1];
	Ai01[3][2]=Hi00[3][1]*Hi01[1][2] + Hi00[3][2]*Hi01[2][2];
	Ai01[3][3]=Hi00[3][3];
	Ai01[3][4]=Hi00[3][4] + Hi00[3][3]*Hi01[3][4];


	Ai02[1][1]=Ai01[1][2]*Hi12[2][1] + Ai01[1][3]*Hi12[3][1];
	Ai02[1][2]=Ai01[1][2]*Hi12[2][2] + Ai01[1][3]*Hi12[3][2];
	Ai02[1][3]=-Ai01[1][1];
	Ai02[1][4]=Ai01[1][4];

	Ai02[2][1]=Ai01[2][2]*Hi12[2][1] + Ai01[2][3]*Hi12[3][1];
	Ai02[2][2]=Ai01[2][2]*Hi12[2][2] + Ai01[2][3]*Hi12[3][2];
	Ai02[2][3]=-Ai01[2][1];
	Ai02[2][4]=Ai01[2][4];

	Ai02[3][1]=Ai01[3][2]*Hi12[2][1] + Ai01[3][3]*Hi12[3][1];
	Ai02[3][2]=Ai01[3][2]*Hi12[2][2] + Ai01[3][3]*Hi12[3][2];
	Ai02[3][3]=-Ai01[3][1];
	Ai02[3][4]=Ai01[3][4];


	Ai03[1][1]=Ai02[1][2]*Hi23[2][1] + Ai02[1][3]*Hi23[3][1];
	Ai03[1][2]=Ai02[1][2]*Hi23[2][2] + Ai02[1][3]*Hi23[3][2];
	Ai03[1][3]=Ai02[1][1];
	Ai03[1][4]=Ai02[1][4] + Ai02[1][1]*Hi23[1][4];

	Ai03[2][1]=Ai02[2][2]*Hi23[2][1] + Ai02[2][3]*Hi23[3][1];
	Ai03[2][2]=Ai02[2][2]*Hi23[2][2] + Ai02[2][3]*Hi23[3][2];
	Ai03[2][3]=Ai02[2][1];
	Ai03[2][4]=Ai02[2][4] + Ai02[2][1]*Hi23[1][4];

	Ai03[3][1]=Ai02[3][2]*Hi23[2][1] + Ai02[3][3]*Hi23[3][1];
	Ai03[3][2]=Ai02[3][2]*Hi23[2][2] + Ai02[3][3]*Hi23[3][2];
	Ai03[3][3]=Ai02[3][1];
	Ai03[3][4]=Ai02[3][4] + Ai02[3][1]*Hi23[1][4];


	Ai04[1][1]=Ai03[1][2]*Hi34[2][1] + Ai03[1][3]*Hi34[3][1];
	Ai04[1][2]=Ai03[1][2]*Hi34[2][2] + Ai03[1][3]*Hi34[3][2];
	Ai04[1][3]=-Ai03[1][1];
	Ai04[1][4]=Ai03[1][4] + Ai03[1][2]*Hi34[2][4] + Ai03[1][3]*Hi34[3][4];

	Ai04[2][1]=Ai03[2][2]*Hi34[2][1] + Ai03[2][3]*Hi34[3][1];
	Ai04[2][2]=Ai03[2][2]*Hi34[2][2] + Ai03[2][3]*Hi34[3][2];
	Ai04[2][3]=-Ai03[2][1];
	Ai04[2][4]=Ai03[2][4] + Ai03[2][2]*Hi34[2][4] + Ai03[2][3]*Hi34[3][4];

	Ai04[3][1]=Ai03[3][2]*Hi34[2][1] + Ai03[3][3]*Hi34[3][1];
	Ai04[3][2]=Ai03[3][2]*Hi34[2][2] + Ai03[3][3]*Hi34[3][2];
	Ai04[3][3]=-Ai03[3][1];
	Ai04[3][4]=Ai03[3][4] + Ai03[3][2]*Hi34[2][4] + Ai03[3][3]*Hi34[3][4];


	Ai05[1][1]=Ai04[1][2]*Hi45[2][1] + Ai04[1][3]*Hi45[3][1];
	Ai05[1][2]=Ai04[1][2]*Hi45[2][2] + Ai04[1][3]*Hi45[3][2];
	Ai05[1][3]=Ai04[1][1];
	Ai05[1][4]=Ai04[1][4] + Ai04[1][1]*Hi45[1][4] + Ai04[1][2]*Hi45[2][4];

	Ai05[2][1]=Ai04[2][2]*Hi45[2][1] + Ai04[2][3]*Hi45[3][1];
	Ai05[2][2]=Ai04[2][2]*Hi45[2][2] + Ai04[2][3]*Hi45[3][2];
	Ai05[2][3]=Ai04[2][1];
	Ai05[2][4]=Ai04[2][4] + Ai04[2][1]*Hi45[1][4] + Ai04[2][2]*Hi45[2][4];

	Ai05[3][1]=Ai04[3][2]*Hi45[2][1] + Ai04[3][3]*Hi45[3][1];
	Ai05[3][2]=Ai04[3][2]*Hi45[2][2] + Ai04[3][3]*Hi45[3][2];
	Ai05[3][3]=Ai04[3][1];
	Ai05[3][4]=Ai04[3][4] + Ai04[3][1]*Hi45[1][4] + Ai04[3][2]*Hi45[2][4];


	Ai06[1][1]=Ai05[1][2]*Hi56[2][1] + Ai05[1][3]*Hi56[3][1];
	Ai06[1][2]=Ai05[1][2]*Hi56[2][2] + Ai05[1][3]*Hi56[3][2];
	Ai06[1][3]=-Ai05[1][1];
	Ai06[1][4]=Ai05[1][4] + Ai05[1][3]*Hi56[3][4];

	Ai06[2][1]=Ai05[2][2]*Hi56[2][1] + Ai05[2][3]*Hi56[3][1];
	Ai06[2][2]=Ai05[2][2]*Hi56[2][2] + Ai05[2][3]*Hi56[3][2];
	Ai06[2][3]=-Ai05[2][1];
	Ai06[2][4]=Ai05[2][4] + Ai05[2][3]*Hi56[3][4];

	Ai06[3][1]=Ai05[3][2]*Hi56[2][1] + Ai05[3][3]*Hi56[3][1];
	Ai06[3][2]=Ai05[3][2]*Hi56[2][2] + Ai05[3][3]*Hi56[3][2];
	Ai06[3][3]=-Ai05[3][1];
	Ai06[3][4]=Ai05[3][4] + Ai05[3][3]*Hi56[3][4];


	Ai07[1][1]=Ai06[1][2]*Hi67[2][1] + Ai06[1][3]*Hi67[3][1];
	Ai07[1][2]=Ai06[1][2]*Hi67[2][2] + Ai06[1][3]*Hi67[3][2];
	Ai07[1][3]=Ai06[1][1];
	Ai07[1][4]=Ai06[1][4];

	Ai07[2][1]=Ai06[2][2]*Hi67[2][1] + Ai06[2][3]*Hi67[3][1];
	Ai07[2][2]=Ai06[2][2]*Hi67[2][2] + Ai06[2][3]*Hi67[3][2];
	Ai07[2][3]=Ai06[2][1];
	Ai07[2][4]=Ai06[2][4];

	Ai07[3][1]=Ai06[3][2]*Hi67[2][1] + Ai06[3][3]*Hi67[3][1];
	Ai07[3][2]=Ai06[3][2]*Hi67[2][2] + Ai06[3][3]*Hi67[3][2];
	Ai07[3][3]=Ai06[3][1];
	Ai07[3][4]=Ai06[3][4];


	Ai08[1][1]=Ai07[1][1]*Hi78[1][1] + Ai07[1][2]*Hi78[2][1] + Ai07[1][3]*Hi78[3][1];
	Ai08[1][2]=Ai07[1][1]*Hi78[1][2] + Ai07[1][2]*Hi78[2][2] + Ai07[1][3]*Hi78[3][2];
	Ai08[1][3]=Ai07[1][1]*Hi78[1][3] + Ai07[1][2]*Hi78[2][3] + Ai07[1][3]*Hi78[3][3];
	Ai08[1][4]=Ai07[1][4] + Ai07[1][1]*Hi78[1][4] + Ai07[1][2]*Hi78[2][4] + Ai07[1][3]*Hi78[3][4];

	Ai08[2][1]=Ai07[2][1]*Hi78[1][1] + Ai07[2][2]*Hi78[2][1] + Ai07[2][3]*Hi78[3][1];
	Ai08[2][2]=Ai07[2][1]*Hi78[1][2] + Ai07[2][2]*Hi78[2][2] + Ai07[2][3]*Hi78[3][2];
	Ai08[2][3]=Ai07[2][1]*Hi78[1][3] + Ai07[2][2]*Hi78[2][3] + Ai07[2][3]*Hi78[3][3];
	Ai08[2][4]=Ai07[2][4] + Ai07[2][1]*Hi78[1][4] + Ai07[2][2]*Hi78[2][4] + Ai07[2][3]*Hi78[3][4];

	Ai08[3][1]=Ai07[3][1]*Hi78[1][1] + Ai07[3][2]*Hi78[2][1] + Ai07[3][3]*Hi78[3][1];
	Ai08[3][2]=Ai07[3][1]*Hi78[1][2] + Ai07[3][2]*Hi78[2][2] + Ai07[3][3]*Hi78[3][2];
	Ai08[3][3]=Ai07[3][1]*Hi78[1][3] + Ai07[3][2]*Hi78[2][3] + Ai07[3][3]*Hi78[3][3];
	Ai08[3][4]=Ai07[3][4] + Ai07[3][1]*Hi78[1][4] + Ai07[3][2]*Hi78[2][4] + Ai07[3][3]*Hi78[3][4];



	/* joint ID: 0 */
	Xorigin[0][1]=Hi00[1][4];
	Xorigin[0][2]=Hi00[2][4];
	Xorigin[0][3]=Hi00[3][4];

	/* link: {basec$0$$x[[1]], basec$0$$x[[2]], basec$0$$x[[3]]} */
	Xlink[0][1]=Hi00[1][4];
	Xlink[0][2]=Hi00[2][4];
	Xlink[0][3]=Hi00[3][4];

	Ahmat[0][1][1]=Hi00[1][1];
	Ahmat[0][1][2]=Hi00[1][2];
	Ahmat[0][1][3]=Hi00[1][3];
	Ahmat[0][1][4]=Hi00[1][4];

	Ahmat[0][2][1]=Hi00[2][1];
	Ahmat[0][2][2]=Hi00[2][2];
	Ahmat[0][2][3]=Hi00[2][3];
	Ahmat[0][2][4]=Hi00[2][4];

	Ahmat[0][3][1]=Hi00[3][1];
	Ahmat[0][3][2]=Hi00[3][2];
	Ahmat[0][3][3]=Hi00[3][3];
	Ahmat[0][3][4]=Hi00[3][4];

	Ahmat[0][4][4]=1;


	/* joint ID: 1 */
	Xorigin[1][1]=Ai01[1][4];
	Xorigin[1][2]=Ai01[2][4];
	Xorigin[1][3]=Ai01[3][4];

	Xaxis[1][1]=Ai01[1][3];
	Xaxis[1][2]=Ai01[2][3];
	Xaxis[1][3]=Ai01[3][3];

	/* link: {0, 0, ZSFE} */
	Xlink[1][1]=Ai01[1][4];
	Xlink[1][2]=Ai01[2][4];
	Xlink[1][3]=Ai01[3][4];

	Ahmat[1][1][1]=Ai02[1][1];
	Ahmat[1][1][2]=Ai02[1][2];
	Ahmat[1][1][3]=Ai02[1][3];
	Ahmat[1][1][4]=Ai02[1][4];

	Ahmat[1][2][1]=Ai02[2][1];
	Ahmat[1][2][2]=Ai02[2][2];
	Ahmat[1][2][3]=Ai02[2][3];
	Ahmat[1][2][4]=Ai02[2][4];

	Ahmat[1][3][1]=Ai02[3][1];
	Ahmat[1][3][2]=Ai02[3][2];
	Ahmat[1][3][3]=Ai02[3][3];
	Ahmat[1][3][4]=Ai02[3][4];

	Ahmat[1][4][4]=1;


	/* joint ID: 2 */
	Xorigin[2][1]=Ai02[1][4];
	Xorigin[2][2]=Ai02[2][4];
	Xorigin[2][3]=Ai02[3][4];

	Xaxis[2][1]=Ai02[1][3];
	Xaxis[2][2]=Ai02[2][3];
	Xaxis[2][3]=Ai02[3][3];

	/* joint ID: 3 */
	Xorigin[3][1]=Ai03[1][4];
	Xorigin[3][2]=Ai03[2][4];
	Xorigin[3][3]=Ai03[3][4];

	Xaxis[3][1]=Ai03[1][3];
	Xaxis[3][2]=Ai03[2][3];
	Xaxis[3][3]=Ai03[3][3];

	/* link: {ZHR, 0, 0} */
	Xlink[2][1]=Ai03[1][4];
	Xlink[2][2]=Ai03[2][4];
	Xlink[2][3]=Ai03[3][4];

	Ahmat[2][1][1]=Ai03[1][1];
	Ahmat[2][1][2]=Ai03[1][2];
	Ahmat[2][1][3]=Ai03[1][3];
	Ahmat[2][1][4]=Ai03[1][4];

	Ahmat[2][2][1]=Ai03[2][1];
	Ahmat[2][2][2]=Ai03[2][2];
	Ahmat[2][2][3]=Ai03[2][3];
	Ahmat[2][2][4]=Ai03[2][4];

	Ahmat[2][3][1]=Ai03[3][1];
	Ahmat[2][3][2]=Ai03[3][2];
	Ahmat[2][3][3]=Ai03[3][3];
	Ahmat[2][3][4]=Ai03[3][4];

	Ahmat[2][4][4]=1;


	/* joint ID: 4 */
	Xorigin[4][1]=Ai04[1][4];
	Xorigin[4][2]=Ai04[2][4];
	Xorigin[4][3]=Ai04[3][4];

	Xaxis[4][1]=Ai04[1][3];
	Xaxis[4][2]=Ai04[2][3];
	Xaxis[4][3]=Ai04[3][3];

	/* link: {0, YEB, ZEB} */
	Xlink[3][1]=Ai04[1][4];
	Xlink[3][2]=Ai04[2][4];
	Xlink[3][3]=Ai04[3][4];

	Ahmat[3][1][1]=Ai04[1][1];
	Ahmat[3][1][2]=Ai04[1][2];
	Ahmat[3][1][3]=Ai04[1][3];
	Ahmat[3][1][4]=Ai04[1][4];

	Ahmat[3][2][1]=Ai04[2][1];
	Ahmat[3][2][2]=Ai04[2][2];
	Ahmat[3][2][3]=Ai04[2][3];
	Ahmat[3][2][4]=Ai04[2][4];

	Ahmat[3][3][1]=Ai04[3][1];
	Ahmat[3][3][2]=Ai04[3][2];
	Ahmat[3][3][3]=Ai04[3][3];
	Ahmat[3][3][4]=Ai04[3][4];

	Ahmat[3][4][4]=1;


	/* joint ID: 5 */
	Xorigin[5][1]=Ai05[1][4];
	Xorigin[5][2]=Ai05[2][4];
	Xorigin[5][3]=Ai05[3][4];

	Xaxis[5][1]=Ai05[1][3];
	Xaxis[5][2]=Ai05[2][3];
	Xaxis[5][3]=Ai05[3][3];

	/* link: {ZWR, YWR, 0} */
	Xlink[4][1]=Ai05[1][4];
	Xlink[4][2]=Ai05[2][4];
	Xlink[4][3]=Ai05[3][4];

	Ahmat[4][1][1]=Ai05[1][1];
	Ahmat[4][1][2]=Ai05[1][2];
	Ahmat[4][1][3]=Ai05[1][3];
	Ahmat[4][1][4]=Ai05[1][4];

	Ahmat[4][2][1]=Ai05[2][1];
	Ahmat[4][2][2]=Ai05[2][2];
	Ahmat[4][2][3]=Ai05[2][3];
	Ahmat[4][2][4]=Ai05[2][4];

	Ahmat[4][3][1]=Ai05[3][1];
	Ahmat[4][3][2]=Ai05[3][2];
	Ahmat[4][3][3]=Ai05[3][3];
	Ahmat[4][3][4]=Ai05[3][4];

	Ahmat[4][4][4]=1;


	/* joint ID: 6 */
	Xorigin[6][1]=Ai06[1][4];
	Xorigin[6][2]=Ai06[2][4];
	Xorigin[6][3]=Ai06[3][4];

	Xaxis[6][1]=Ai06[1][3];
	Xaxis[6][2]=Ai06[2][3];
	Xaxis[6][3]=Ai06[3][3];

	/* link: {0, 0, ZWFE} */
	Xlink[5][1]=Ai06[1][4];
	Xlink[5][2]=Ai06[2][4];
	Xlink[5][3]=Ai06[3][4];

	Ahmat[5][1][1]=Ai07[1][1];
	Ahmat[5][1][2]=Ai07[1][2];
	Ahmat[5][1][3]=Ai07[1][3];
	Ahmat[5][1][4]=Ai07[1][4];

	Ahmat[5][2][1]=Ai07[2][1];
	Ahmat[5][2][2]=Ai07[2][2];
	Ahmat[5][2][3]=Ai07[2][3];
	Ahmat[5][2][4]=Ai07[2][4];

	Ahmat[5][3][1]=Ai07[3][1];
	Ahmat[5][3][2]=Ai07[3][2];
	Ahmat[5][3][3]=Ai07[3][3];
	Ahmat[5][3][4]=Ai07[3][4];

	Ahmat[5][4][4]=1;


	/* joint ID: 7 */
	Xorigin[7][1]=Ai07[1][4];
	Xorigin[7][2]=Ai07[2][4];
	Xorigin[7][3]=Ai07[3][4];

	Xaxis[7][1]=Ai07[1][3];
	Xaxis[7][2]=Ai07[2][3];
	Xaxis[7][3]=Ai07[3][3];

	/* link: {eff$1$$x[[1]], eff$1$$x[[2]], eff$1$$x[[3]]} */
	Xlink[6][1]=Ai08[1][4];
	Xlink[6][2]=Ai08[2][4];
	Xlink[6][3]=Ai08[3][4];

	Ahmat[6][1][1]=Ai08[1][1];
	Ahmat[6][1][2]=Ai08[1][2];
	Ahmat[6][1][3]=Ai08[1][3];
	Ahmat[6][1][4]=Ai08[1][4];

	Ahmat[6][2][1]=Ai08[2][1];
	Ahmat[6][2][2]=Ai08[2][2];
	Ahmat[6][2][3]=Ai08[2][3];
	Ahmat[6][2][4]=Ai08[2][4];

	Ahmat[6][3][1]=Ai08[3][1];
	Ahmat[6][3][2]=Ai08[3][2];
	Ahmat[6][3][3]=Ai08[3][3];
	Ahmat[6][3][4]=Ai08[3][4];

	Ahmat[6][4][4]=1;

}

static void jacobian(const double link[NLINK+1][4],
		const double origin[NDOF+1][4],
		const double axis[NDOF+1][4],
		double jac[2*NCART][NDOF]) {

	static const int PALM = 6;
	static double c[2*NCART];
	for (int j = 1; j <= NDOF; ++j) {
		c[0] = axis[j][2] * (link[PALM][3] - origin[j][3]) - axis[j][3] * (link[PALM][2]-origin[j][2]);
		c[1] = axis[j][3] * (link[PALM][1] - origin[j][1]) - axis[j][1] * (link[PALM][3]-origin[j][3]);
		c[2] = axis[j][1] * (link[PALM][2] - origin[j][2]) - axis[j][2] * (link[PALM][1]-origin[j][1]);
		c[3] = axis[j][1];
		c[4] = axis[j][2];
		c[5] = axis[j][3];
		//rev_geo_jac_col(lp[PALM], jop[j], jap[j], c);
		for (int i = 0; i < 2*NCART; i++)
			jac[i][j-1] = c[i];
	}
}

/*
 * Copied from SL_user_common.c for convenience
 *
 */
static void set_endeffector(double endeff_pos[NCART]) {

	endeff_pos[X]  = 0.0;
	endeff_pos[Y]  = 0.0;
	endeff_pos[Z]  = 0.30;
	// attach the racket

}

/**
 *
 * @brief Reads joint limits from file.
 *
 * @param lb Array of lower bound values to be loaded.
 * @param ub Array of upper bound values to be loaded.
 * @return If can load the joint limits successfully then returns 1.
 *
 */
int read_joint_limits(double *lb, double *ub) {

	char joint_names[][20] = {
			{"R_SFE"},
			{"R_SAA"},
			{"R_HR"},
			{"R_EB"},
			{"R_WR"},
			{"R_WFE"},
			{"R_WAA"}
	};
	char fname[] = "SensorOffset.cf";

	/* find all joint variables and read them into the appropriate array */

	char string[100];
	FILE *in;

	/* get the max, min of the position sensors */

	sprintf(string,"%s/robolab/barrett/%s%s",getenv("HOME"),CONFIG,fname);
	in = fopen(string,"r");
	if (in == NULL) {
		printf("ERROR: Cannot open file >%s<!\n",string);
		return FALSE;
	}

	/* find all joint variables and read them into the appropriate array */

	for (int i = 0; i < NDOF; i++) {
		if (!find_keyword(in, &(joint_names[i][0]))) {
			printf("ERROR: Cannot find offset for %s!\n",joint_names[i]);
			fclose(in);
			return TRUE;
		}
		fscanf(in,"%lf %lf", &lb[i], &ub[i]);
	}
	fclose(in);

	return TRUE;

}
