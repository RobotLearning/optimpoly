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
#include "kinematics.h"
#include "utils.h"
#include "string.h"

SL_Cstate     cart_state[N_ENDEFFS+1];        /* endeffector state */
SL_quat       cart_orient[N_ENDEFFS+1];       /* endeffector orientation */
SL_endeff     endeff[N_ENDEFFS+1];            /* endeffector structure */
SL_Cstate     base_state;                     /* cartesian state of base coordinate system */
SL_quat       base_orient;                    /* cartesian orientation of base coordinate system */
Matrix        Alink_des[N_LINKS+1];           /* homogeneous transformation matrices for all links */
Matrix        link_pos_des;                   /* desired cart. pos of links */
Matrix        joint_cog_mpos_des;             /* vector of mass*COG of each joint based on desireds*/
Matrix       joint_origin_pos_des;           /* vector of pos. of local joint coord.sys based on des.*/
Matrix        joint_axis_pos_des;             /* unit vector of joint rotation axis based on des.*/
SL_DJstate    joint_default_state[N_DOFS+1];
SL_OJstate    joint_opt_state[N_DOFS+1];
SL_Jstate     joint_state[N_DOFS+1];
SL_DJstate    joint_des_state[N_DOFS+1];
double        joint_range[N_DOFS+1][3+1];
SL_link       links[N_DOFS+1];                /* specs of links: mass, inertia, cm */

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
void
kinematics(SL_DJstate *state,SL_Cstate *basec, SL_quat *baseo, SL_endeff *eff,
		   double **Xmcog, double **Xaxis, double **Xorigin, double **Xlink, double ***Ahmat) {

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
	sstate1th=Sin(state[1].th);
	cstate1th=Cos(state[1].th);

	sstate2th=Sin(state[2].th);
	cstate2th=Cos(state[2].th);

	sstate3th=Sin(state[3].th);
	cstate3th=Cos(state[3].th);

	sstate4th=Sin(state[4].th);
	cstate4th=Cos(state[4].th);

	sstate5th=Sin(state[5].th);
	cstate5th=Cos(state[5].th);

	sstate6th=Sin(state[6].th);
	cstate6th=Cos(state[6].th);

	sstate7th=Sin(state[7].th);
	cstate7th=Cos(state[7].th);


	/* rotation matrix sine and cosine precomputation */







	rseff1a1=Sin(eff[1].a[1]);
	rceff1a1=Cos(eff[1].a[1]);

	rseff1a2=Sin(eff[1].a[2]);
	rceff1a2=Cos(eff[1].a[2]);

	rseff1a3=Sin(eff[1].a[3]);
	rceff1a3=Cos(eff[1].a[3]);



	/* inverse homogeneous rotation matrices */
	Hi00[1][1]=-1 + 2*Power(baseo[0].q[1],2) + 2*Power(baseo[0].q[2],2);
	Hi00[1][2]=2*(baseo[0].q[2]*baseo[0].q[3] - baseo[0].q[1]*baseo[0].q[4]);
	Hi00[1][3]=2*(baseo[0].q[1]*baseo[0].q[3] + baseo[0].q[2]*baseo[0].q[4]);
	Hi00[1][4]=basec[0].x[1];

	Hi00[2][1]=2*(baseo[0].q[2]*baseo[0].q[3] + baseo[0].q[1]*baseo[0].q[4]);
	Hi00[2][2]=-1 + 2*Power(baseo[0].q[1],2) + 2*Power(baseo[0].q[3],2);
	Hi00[2][3]=2*(-(baseo[0].q[1]*baseo[0].q[2]) + baseo[0].q[3]*baseo[0].q[4]);
	Hi00[2][4]=basec[0].x[2];

	Hi00[3][1]=2*(-(baseo[0].q[1]*baseo[0].q[3]) + baseo[0].q[2]*baseo[0].q[4]);
	Hi00[3][2]=2*(baseo[0].q[1]*baseo[0].q[2] + baseo[0].q[3]*baseo[0].q[4]);
	Hi00[3][3]=-1 + 2*Power(baseo[0].q[1],2) + 2*Power(baseo[0].q[4],2);
	Hi00[3][4]=basec[0].x[3];


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
	Hi78[1][4]=eff[1].x[1];

	Hi78[2][1]=rceff1a3*rseff1a1*rseff1a2 + rceff1a1*rseff1a3;
	Hi78[2][2]=rceff1a1*rceff1a3 - rseff1a1*rseff1a2*rseff1a3;
	Hi78[2][3]=-(rceff1a2*rseff1a1);
	Hi78[2][4]=eff[1].x[2];

	Hi78[3][1]=-(rceff1a1*rceff1a3*rseff1a2) + rseff1a1*rseff1a3;
	Hi78[3][2]=rceff1a3*rseff1a1 + rceff1a1*rseff1a2*rseff1a3;
	Hi78[3][3]=rceff1a1*rceff1a2;
	Hi78[3][4]=eff[1].x[3];



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

	Xmcog[0][1]=links[0].mcm[1]*Hi00[1][1] + links[0].mcm[2]*Hi00[1][2] + links[0].mcm[3]*Hi00[1][3] + links[0].m*Hi00[1][4];
	Xmcog[0][2]=links[0].mcm[1]*Hi00[2][1] + links[0].mcm[2]*Hi00[2][2] + links[0].mcm[3]*Hi00[2][3] + links[0].m*Hi00[2][4];
	Xmcog[0][3]=links[0].mcm[1]*Hi00[3][1] + links[0].mcm[2]*Hi00[3][2] + links[0].mcm[3]*Hi00[3][3] + links[0].m*Hi00[3][4];

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

	Xmcog[1][1]=links[1].mcm[1]*Ai01[1][1] + links[1].mcm[2]*Ai01[1][2] + links[1].mcm[3]*Ai01[1][3] + links[1].m*Ai01[1][4];
	Xmcog[1][2]=links[1].mcm[1]*Ai01[2][1] + links[1].mcm[2]*Ai01[2][2] + links[1].mcm[3]*Ai01[2][3] + links[1].m*Ai01[2][4];
	Xmcog[1][3]=links[1].mcm[1]*Ai01[3][1] + links[1].mcm[2]*Ai01[3][2] + links[1].mcm[3]*Ai01[3][3] + links[1].m*Ai01[3][4];

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

	Xmcog[2][1]=links[2].mcm[1]*Ai02[1][1] + links[2].mcm[2]*Ai02[1][2] + links[2].mcm[3]*Ai02[1][3] + links[2].m*Ai02[1][4];
	Xmcog[2][2]=links[2].mcm[1]*Ai02[2][1] + links[2].mcm[2]*Ai02[2][2] + links[2].mcm[3]*Ai02[2][3] + links[2].m*Ai02[2][4];
	Xmcog[2][3]=links[2].mcm[1]*Ai02[3][1] + links[2].mcm[2]*Ai02[3][2] + links[2].mcm[3]*Ai02[3][3] + links[2].m*Ai02[3][4];

	Xaxis[2][1]=Ai02[1][3];
	Xaxis[2][2]=Ai02[2][3];
	Xaxis[2][3]=Ai02[3][3];

	/* joint ID: 3 */
	Xorigin[3][1]=Ai03[1][4];
	Xorigin[3][2]=Ai03[2][4];
	Xorigin[3][3]=Ai03[3][4];

	Xmcog[3][1]=links[3].mcm[1]*Ai03[1][1] + links[3].mcm[2]*Ai03[1][2] + links[3].mcm[3]*Ai03[1][3] + links[3].m*Ai03[1][4];
	Xmcog[3][2]=links[3].mcm[1]*Ai03[2][1] + links[3].mcm[2]*Ai03[2][2] + links[3].mcm[3]*Ai03[2][3] + links[3].m*Ai03[2][4];
	Xmcog[3][3]=links[3].mcm[1]*Ai03[3][1] + links[3].mcm[2]*Ai03[3][2] + links[3].mcm[3]*Ai03[3][3] + links[3].m*Ai03[3][4];

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

	Xmcog[4][1]=links[4].mcm[1]*Ai04[1][1] + links[4].mcm[2]*Ai04[1][2] + links[4].mcm[3]*Ai04[1][3] + links[4].m*Ai04[1][4];
	Xmcog[4][2]=links[4].mcm[1]*Ai04[2][1] + links[4].mcm[2]*Ai04[2][2] + links[4].mcm[3]*Ai04[2][3] + links[4].m*Ai04[2][4];
	Xmcog[4][3]=links[4].mcm[1]*Ai04[3][1] + links[4].mcm[2]*Ai04[3][2] + links[4].mcm[3]*Ai04[3][3] + links[4].m*Ai04[3][4];

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

	Xmcog[5][1]=links[5].mcm[1]*Ai05[1][1] + links[5].mcm[2]*Ai05[1][2] + links[5].mcm[3]*Ai05[1][3] + links[5].m*Ai05[1][4];
	Xmcog[5][2]=links[5].mcm[1]*Ai05[2][1] + links[5].mcm[2]*Ai05[2][2] + links[5].mcm[3]*Ai05[2][3] + links[5].m*Ai05[2][4];
	Xmcog[5][3]=links[5].mcm[1]*Ai05[3][1] + links[5].mcm[2]*Ai05[3][2] + links[5].mcm[3]*Ai05[3][3] + links[5].m*Ai05[3][4];

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

	Xmcog[6][1]=links[6].mcm[1]*Ai06[1][1] + links[6].mcm[2]*Ai06[1][2] + links[6].mcm[3]*Ai06[1][3] + links[6].m*Ai06[1][4];
	Xmcog[6][2]=links[6].mcm[1]*Ai06[2][1] + links[6].mcm[2]*Ai06[2][2] + links[6].mcm[3]*Ai06[2][3] + links[6].m*Ai06[2][4];
	Xmcog[6][3]=links[6].mcm[1]*Ai06[3][1] + links[6].mcm[2]*Ai06[3][2] + links[6].mcm[3]*Ai06[3][3] + links[6].m*Ai06[3][4];

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

	Xmcog[7][1]=links[7].mcm[1]*Ai07[1][1] + links[7].mcm[2]*Ai07[1][2] + links[7].mcm[3]*Ai07[1][3] + links[7].m*Ai07[1][4];
	Xmcog[7][2]=links[7].mcm[1]*Ai07[2][1] + links[7].mcm[2]*Ai07[2][2] + links[7].mcm[3]*Ai07[2][3] + links[7].m*Ai07[2][4];
	Xmcog[7][3]=links[7].mcm[1]*Ai07[3][1] + links[7].mcm[2]*Ai07[3][2] + links[7].mcm[3]*Ai07[3][3] + links[7].m*Ai07[3][4];

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

/*
 * Computes the jacobian. Take from SL and simplified. *
 *
 * Function Parameters: [in]=input,[out]=output
 *
 * \param[in]     lp      : the link positions
 * \param[in]     jop     : joint origin positions
 * \param[in]     jap     : joint axix unit vectors
 * \param[out]    Jac     : the jacobian
 *
 */
void jacobian(Matrix lp, Matrix jop, Matrix jap, Matrix Jac) {

	int i,j;
	static const int PALM = 6;
	double c[2*CART+1];
	for (i = 1; i <= DOF; ++i) {
		revoluteGJacColumn(lp[PALM], jop[i], jap[i], c);
		for (j = 1; j <= 2*CART; ++j)
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
void revoluteGJacColumn(Vector p, Vector pi, Vector zi, Vector c) {

  c[1] = zi[2] * (p[3]-pi[3]) - zi[3] * (p[2]-pi[2]);
  c[2] = zi[3] * (p[1]-pi[1]) - zi[1] * (p[3]-pi[3]);
  c[3] = zi[1] * (p[2]-pi[2]) - zi[2] * (p[1]-pi[1]);
  c[4] = zi[1];
  c[5] = zi[2];
  c[6] = zi[3];

}

/*
 * Copied from SL_user_common.c for convenience
 *
 */
void setDefaultEndeffector(void) {


	endeff[1].m       = 0.0;
	endeff[1].mcm[_X_]= 0.0;
	endeff[1].mcm[_Y_]= 0.0;
	endeff[1].mcm[_Z_]= 0.0;
	endeff[1].x[_X_]  = 0.0;
	endeff[1].x[_Y_]  = 0.0;
	endeff[1].x[_Z_]  = 0.08531+0.06;
	endeff[1].a[_X_]  = 0.0;
	endeff[1].a[_Y_]  = 0.0;
	endeff[1].a[_Z_]  = 0.0;

	// attach the racket
	endeff[1].x[_Z_] = .3;

}

/*
 * Load the joint limits from config/ file into joint_range array
 *
 */
void load_joint_limits() {

	char *fname = "SensorOffset.cf";
	read_sensor_offsets(fname);

}

/*
 * Copied from SL_common. Dont want to call the function from SL because
 * we have to include a lot of extra SL files
 *
 */
int read_sensor_offsets(char *fname) {

  char string[100];
  FILE  *in;

  char joint_names[][20]= {
    {"dummy"},
    {"R_SFE"},
    {"R_SAA"},
    {"R_HR"},
    {"R_EB"},
    {"R_WR"},
    {"R_WFE"},
    {"R_WAA"}
  };

  /* get the max, min, and offsets of the position sensors */

  sprintf(string,"%s/robolab/barrett/%s%s",getenv("HOME"),CONFIG,fname);
  in = fopen(string,"r");
  if (in == NULL) {
    printf("ERROR: Cannot open file >%s<!\n",string);
    return FALSE;
  }

  /* find all joint variables and read them into the appropriate array */

  for (int i=1; i<= N_DOFS; ++i) {
    if (!find_keyword(in, &(joint_names[i][0]))) {
      printf("ERROR: Cannot find offset for >%s<!\n",joint_names[i]);
      fclose(in);
      return FALSE;
    }
    fscanf(in,"%lf %lf %lf %lf %lf %lf",
	&joint_range[i][MIN_THETA], &joint_range[i][MAX_THETA],
	   &(joint_default_state[i].th),
	   &(joint_opt_state[i].th),
	   &(joint_opt_state[i].w),
	   &joint_range[i][THETA_OFFSET]);
    joint_default_state[i].thd = 0;
    joint_default_state[i].uff = 0;
  }

  fclose(in);

  return TRUE;

}
