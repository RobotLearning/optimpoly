/*
 * kinematics.c
 *
 * Here we include the kinematics related functions taken from SL
 *
 *  Created on: Jun 22, 2016
 *      Author: okoc
 */

#include "stdio.h"
#include "stdlib.h"
#include "SL.h"
#include "constants.h"
#include "kinematics.h"
#include "utils.h"
#include "string.h"

/*
 * Calculates cartesian racket pos, vel and normal
 * given joint positions and velocities
 *
 */
void calc_racket_state(const double q[NDOF],
		               const double qdot[NDOF],
					   double pos[NCART],
					   double vel[NCART],
					   double normal[NCART]) {

	static int firsttime = TRUE;
	static const int PALM = 5;
	static double base_orient[4]; // quat
	static double base_state[3]; // pos
	static double eff_angles[3]; // initial euler angles for racket
	static double eff_pos[3]; // initial racket positions

	static double link[N_LINKS][3];
	static double origin[NDOF+1][3];
	static double axis[NDOF][3];
	static double amats[N_LINKS][4][4];
	static double jacobi[NCART][NDOF];

	/* initialization of static variables */
	if (firsttime) {
		firsttime = FALSE;
		base_orient[1] = 1.0;
		set_endeffector(eff_pos); // set racket
	}

	kinematics(q,base_state,base_orient,eff_angles,eff_pos,
			     link,origin,axis,amats);
	//rotate_to_quat(amats.slice(PALM)(span(X,Z),span(X,Z)),quat);
	//calc_racket_orient(quat);
	jacobian(link[PALM],origin,axis,jacobi);
	pos = link[PALM];
	get_cart_velocity(jacobi,qdot,vel);
	normal = amats[PALM][1];
}

/*
 * Find cartesian racket velocity given
 * joint velocities and Jacobian
 */
void get_cart_velocity(const double jac[6][7],
		                  const double qdot[NDOF],
						  double vel[NCART]) {

	for (int i = 0; i < NCART; i++) {
		vel[i] = 0.0;
		for (int j = 0; j < NDOF; j++) {
			vel[i] += jac[i][j] * qdot[j];
		}
	}
}


/*!*****************************************************************************
 *******************************************************************************
\note  linkInformationDes
\date  March 2005

\remarks

        computes the m*cog, rotation axis, and local coord.sys. orgin for
        every link. This information can be used to compute the COG and
        COG jacobians, assuming the mass and center of mass parameters are
        properly identified.

 *******************************************************************************
 Function Parameters: [in]=input,[out]=output

 \param[in]     state   : the state containing th, thd, thdd, and receiving the
                  appropriate u
 \param[in]     basec   : the position state of the base
 \param[in]     baseo   : the orientational state of the base
 \param[in]     endeff  : the endeffector parameters
 \param[out]    Xmcog   : array of mass*cog vectors
 \param[out]    Xaxis   : array of rotation axes
 \param[out]    Xorigin : array of coord.sys. origin vectors
 \param[out]    Xlink   : array of link position
 \param[out]    Ahmat   : homogeneous transformation matrices of each link

 ******************************************************************************/
void kinematics(const double state[NDOF],
		        const double basec[NCART], const double baseo[NQUAT],
		        const double eff_a[NCART], const double eff_x[NCART],
		        double Xaxis[NDOF][3],
				double Xorigin[NDOF+1][3],
				double Xlink[N_LINKS][3],
				double Ahmat[N_LINKS][4][4]) {

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


	rseff1a1=Sin(eff_a[0]);
	rceff1a1=Cos(eff_a[0]);

	rseff1a2=Sin(eff_a[1]);
	rceff1a2=Cos(eff_a[1]);

	rseff1a3=Sin(eff_a[2]);
	rceff1a3=Cos(eff_a[2]);



	/* inverse homogeneous rotation matrices */
	Hi00[1][1]=-1 + 2*Power(baseo[0],2) + 2*Power(baseo[1],2);
	Hi00[1][2]=2*(baseo[1]*baseo[2] - baseo[0]*baseo[3]);
	Hi00[1][3]=2*(baseo[0]*baseo[2] + baseo[1]*baseo[3]);
	Hi00[1][4]=basec[0];

	Hi00[2][1]=2*(baseo[1]*baseo[2] + baseo[0]*baseo[3]);
	Hi00[2][2]=-1 + 2*Power(baseo[0],2) + 2*Power(baseo[2],2);
	Hi00[2][3]=2*(-(baseo[0]*baseo[1]) + baseo[2]*baseo[3]);
	Hi00[2][4]=basec[1];

	Hi00[3][1]=2*(-(baseo[0]*baseo[2]) + baseo[1]*baseo[3]);
	Hi00[3][2]=2*(baseo[0]*baseo[1] + baseo[2]*baseo[3]);
	Hi00[3][3]=-1 + 2*Power(baseo[0],2) + 2*Power(baseo[3],2);
	Hi00[3][4]=basec[2];


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
	Hi78[1][4]=eff_x[0];

	Hi78[2][1]=rceff1a3*rseff1a1*rseff1a2 + rceff1a1*rseff1a3;
	Hi78[2][2]=rceff1a1*rceff1a3 - rseff1a1*rseff1a2*rseff1a3;
	Hi78[2][3]=-(rceff1a2*rseff1a1);
	Hi78[2][4]=eff_x[1];

	Hi78[3][1]=-(rceff1a1*rceff1a3*rseff1a2) + rseff1a1*rseff1a3;
	Hi78[3][2]=rceff1a3*rseff1a1 + rceff1a1*rseff1a2*rseff1a3;
	Hi78[3][3]=rceff1a1*rceff1a2;
	Hi78[3][4]=eff_x[2];


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
	Xorigin[0][X]=Hi00[1][4];
	Xorigin[0][Y]=Hi00[2][4];
	Xorigin[0][Z]=Hi00[3][4];

	/* link: {basec$0$$x[[1]], basec$0$$x[[2]], basec$0$$x[[3]]} */
	Xlink[0][X]=Hi00[1][4];
	Xlink[0][Y]=Hi00[2][4];
	Xlink[0][Z]=Hi00[3][4];

	Ahmat[0][X][0]=Hi00[1][1];
	Ahmat[0][X][1]=Hi00[1][2];
	Ahmat[0][X][2]=Hi00[1][3];
	Ahmat[0][X][3]=Hi00[1][4];

	Ahmat[0][Y][0]=Hi00[2][1];
	Ahmat[0][Y][1]=Hi00[2][2];
	Ahmat[0][Y][2]=Hi00[2][3];
	Ahmat[0][Y][3]=Hi00[2][4];

	Ahmat[0][Z][0]=Hi00[3][1];
	Ahmat[0][Z][1]=Hi00[3][2];
	Ahmat[0][Z][2]=Hi00[3][3];
	Ahmat[0][Z][3]=Hi00[3][4];

	Ahmat[0][3][3]=1;


	/* joint ID: 1 */
	Xorigin[1][X]=Ai01[1][4];
	Xorigin[1][Y]=Ai01[2][4];
	Xorigin[1][Z]=Ai01[3][4];

	Xaxis[0][X]=Ai01[1][3];
	Xaxis[0][Y]=Ai01[2][3];
	Xaxis[0][Z]=Ai01[3][3];

	/* link: {0, 0, ZSFE} */
	Xlink[1][X]=Ai01[1][4];
	Xlink[1][Y]=Ai01[2][4];
	Xlink[1][Z]=Ai01[3][4];

	Ahmat[1][X][0]=Ai02[1][1];
	Ahmat[1][X][1]=Ai02[1][2];
	Ahmat[1][X][2]=Ai02[1][3];
	Ahmat[1][X][3]=Ai02[1][4];

	Ahmat[1][Y][0]=Ai02[2][1];
	Ahmat[1][Y][1]=Ai02[2][2];
	Ahmat[1][Y][2]=Ai02[2][3];
	Ahmat[1][Y][3]=Ai02[2][4];

	Ahmat[1][Z][0]=Ai02[3][1];
	Ahmat[1][Z][1]=Ai02[3][2];
	Ahmat[1][Z][2]=Ai02[3][3];
	Ahmat[1][Z][3]=Ai02[3][4];

	Ahmat[1][3][3]=1;


	/* joint ID: 2 */
	Xorigin[2][X]=Ai02[1][4];
	Xorigin[2][Y]=Ai02[2][4];
	Xorigin[2][Z]=Ai02[3][4];

	Xaxis[1][X]=Ai02[1][3];
	Xaxis[1][Y]=Ai02[2][3];
	Xaxis[1][Z]=Ai02[3][3];

	/* joint ID: 3 */
	Xorigin[3][X]=Ai03[1][4];
	Xorigin[3][Y]=Ai03[2][4];
	Xorigin[3][Z]=Ai03[3][4];

	Xaxis[2][X]=Ai03[1][3];
	Xaxis[2][Y]=Ai03[2][3];
	Xaxis[2][Z]=Ai03[3][3];

	/* link: {ZHR, 0, 0} */
	Xlink[2][X]=Ai03[1][4];
	Xlink[2][Y]=Ai03[2][4];
	Xlink[2][Z]=Ai03[3][4];

	Ahmat[2][X][0]=Ai03[1][1];
	Ahmat[2][X][1]=Ai03[1][2];
	Ahmat[2][X][2]=Ai03[1][3];
	Ahmat[2][X][3]=Ai03[1][4];

	Ahmat[2][Y][0]=Ai03[2][1];
	Ahmat[2][Y][1]=Ai03[2][2];
	Ahmat[2][Y][2]=Ai03[2][3];
	Ahmat[2][Y][3]=Ai03[2][4];

	Ahmat[2][Z][0]=Ai03[3][1];
	Ahmat[2][Z][1]=Ai03[3][2];
	Ahmat[2][Z][2]=Ai03[3][3];
	Ahmat[2][Z][3]=Ai03[3][4];

	Ahmat[2][3][3]=1;


	/* joint ID: 4 */
	Xorigin[4][X]=Ai04[1][4];
	Xorigin[4][Y]=Ai04[2][4];
	Xorigin[4][Z]=Ai04[3][4];

	Xaxis[3][X]=Ai04[1][3];
	Xaxis[3][Y]=Ai04[2][3];
	Xaxis[3][Z]=Ai04[3][3];

	/* link: {0, YEB, ZEB} */
	Xlink[3][X]=Ai04[1][4];
	Xlink[3][Y]=Ai04[2][4];
	Xlink[3][Z]=Ai04[3][4];

	Ahmat[3][X][0]=Ai04[1][1];
	Ahmat[3][X][1]=Ai04[1][2];
	Ahmat[3][X][2]=Ai04[1][3];
	Ahmat[3][X][3]=Ai04[1][4];

	Ahmat[3][Y][0]=Ai04[2][1];
	Ahmat[3][Y][1]=Ai04[2][2];
	Ahmat[3][Y][2]=Ai04[2][3];
	Ahmat[3][Y][3]=Ai04[2][4];

	Ahmat[3][Z][0]=Ai04[3][1];
	Ahmat[3][Z][1]=Ai04[3][2];
	Ahmat[3][Z][2]=Ai04[3][3];
	Ahmat[3][Z][3]=Ai04[3][4];

	Ahmat[3][3][3]=1;


	/* joint ID: 5 */
	Xorigin[5][X]=Ai05[1][4];
	Xorigin[5][Y]=Ai05[2][4];
	Xorigin[5][Z]=Ai05[3][4];

	Xaxis[4][X]=Ai05[1][3];
	Xaxis[4][Y]=Ai05[2][3];
	Xaxis[4][Z]=Ai05[3][3];

	/* link: {ZWR, YWR, 0} */
	Xlink[4][X]=Ai05[1][4];
	Xlink[4][Y]=Ai05[2][4];
	Xlink[4][Z]=Ai05[3][4];

	Ahmat[4][X][0]=Ai05[1][1];
	Ahmat[4][X][1]=Ai05[1][2];
	Ahmat[4][X][2]=Ai05[1][3];
	Ahmat[4][X][3]=Ai05[1][4];

	Ahmat[4][Y][0]=Ai05[2][1];
	Ahmat[4][Y][1]=Ai05[2][2];
	Ahmat[4][Y][2]=Ai05[2][3];
	Ahmat[4][Y][3]=Ai05[2][4];

	Ahmat[4][Z][0]=Ai05[3][1];
	Ahmat[4][Z][1]=Ai05[3][2];
	Ahmat[4][Z][2]=Ai05[3][3];
	Ahmat[4][Z][3]=Ai05[3][4];

	Ahmat[4][3][3]=1;


	/* joint ID: 6 */
	Xorigin[6][X]=Ai06[1][4];
	Xorigin[6][Y]=Ai06[2][4];
	Xorigin[6][Z]=Ai06[3][4];

	Xaxis[5][X]=Ai06[1][3];
	Xaxis[5][Y]=Ai06[2][3];
	Xaxis[5][Z]=Ai06[3][3];

	/* link: {0, 0, ZWFE} */
	Xlink[5][X]=Ai06[1][4];
	Xlink[5][Y]=Ai06[2][4];
	Xlink[5][Z]=Ai06[3][4];

	Ahmat[5][X][0]=Ai07[1][1];
	Ahmat[5][X][1]=Ai07[1][2];
	Ahmat[5][X][2]=Ai07[1][3];
	Ahmat[5][X][3]=Ai07[1][4];

	Ahmat[5][Y][0]=Ai07[2][1];
	Ahmat[5][Y][1]=Ai07[2][2];
	Ahmat[5][Y][2]=Ai07[2][3];
	Ahmat[5][Y][3]=Ai07[2][4];

	Ahmat[5][Z][0]=Ai07[3][1];
	Ahmat[5][Z][1]=Ai07[3][2];
	Ahmat[5][Z][2]=Ai07[3][3];
	Ahmat[5][Z][3]=Ai07[3][4];

	Ahmat[5][3][3]=1;


	/* joint ID: 7 */
	Xorigin[7][X]=Ai07[1][4];
	Xorigin[7][Y]=Ai07[2][4];
	Xorigin[7][Z]=Ai07[3][4];

	Xaxis[6][X]=Ai07[1][3];
	Xaxis[6][Y]=Ai07[2][3];
	Xaxis[6][Z]=Ai07[3][3];

	/* link: {eff$1$$x[[1]], eff$1$$x[[2]], eff$1$$x[[3]]} */
	Xlink[6][X]=Ai08[1][4];
	Xlink[6][Y]=Ai08[2][4];
	Xlink[6][Z]=Ai08[3][4];

	Ahmat[6][X][0]=Ai08[1][1];
	Ahmat[6][X][1]=Ai08[1][2];
	Ahmat[6][X][2]=Ai08[1][3];
	Ahmat[6][X][3]=Ai08[1][4];

	Ahmat[6][Y][0]=Ai08[2][1];
	Ahmat[6][Y][1]=Ai08[2][2];
	Ahmat[6][Y][2]=Ai08[2][3];
	Ahmat[6][Y][3]=Ai08[2][4];

	Ahmat[6][Z][0]=Ai08[3][1];
	Ahmat[6][Z][1]=Ai08[3][2];
	Ahmat[6][Z][2]=Ai08[3][3];
	Ahmat[6][Z][3]=Ai08[3][4];

	Ahmat[6][3][3]=1;
}

/*
 * Computes the jacobian. Take from SL and simplified. *
 *
 * Function Parameters: [in]=input,[out]=output
 *
 * link      : the link positions
 * origin    : joint origin positions
 * axis      : joint axis unit vectors
 * jac     : the jacobian [out]
 *
 */
void jacobian(const double link[NCART],
		      const double origin[NDOF+1][NCART],
		      const double axis[NDOF][NCART],
			  double Jac[NCART][NDOF]) {

	int i,j;
	double c[2*NCART];
	for (i = 0; i < NDOF; ++i) {
		revolute_geo_jac_col(link, origin[i], axis[i], c);
		for (j = 0; j < NCART; ++j)
			Jac[j][i] = c[j];
	}

}

/*
 *
 * Copied from SL_common. Dont want to call the function from SL because
 * we have to include a lot of extra SL files
 *

 computes one column for the geometric jacobian of a revolute joint
 from the given input vectors

 *******************************************************************************
 Function Parameters: [in]=input,[out]=output

 \param[in]     p    : position of endeffector
 \param[in]     pi   : position of joint origin
 \param[in]     zi   : unit vector of joint axis
 \param[out]    c    : column vector of Jacobian

 ******************************************************************************/
void revolute_geo_jac_col(const double p[NCART],
		                const double pi[NCART],
						const double zi[NCART],
		                double c[NDOF]) {

  c[X] = zi[Y] * (p[Z]-pi[Z]) - zi[Z] * (p[Y]-pi[Y]);
  c[Y] = zi[Z] * (p[X]-pi[X]) - zi[X] * (p[Z]-pi[Z]);
  c[Z] = zi[X] * (p[Y]-pi[Y]) - zi[Y] * (p[X]-pi[X]);
  c[DX] = zi[X];
  c[DY] = zi[Y];
  c[DZ] = zi[Z];

}

/*
 * Copied from SL_user_common.c for convenience
 *
 */
void set_endeffector(double endeff_pos[NCART]) {

	endeff_pos[X]  = 0.0;
	endeff_pos[Y]  = 0.0;
	endeff_pos[Z]  = 0.30;
	// attach the racket

}

/*
 * Copied from SL_common.
 *
 * Reads joint limits from file.
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
