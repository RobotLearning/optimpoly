/*
 * kinematics.h
 *
 * Here we include the kinematics related functions taken from SL
 *
 *  Created on: Jun 22, 2016
 *      Author: okoc
 */

#ifndef KINEMATICS_H_
#define KINEMATICS_H_

#include "math.h"
#include "SL.h"

#define DOF 7
#define CART 3

extern SL_Cstate     cart_state[N_ENDEFFS+1];        /* endeffector state */
extern SL_quat       cart_orient[N_ENDEFFS+1];       /* endeffector orientation */
extern SL_endeff     endeff[N_ENDEFFS+1];            /* endeffector structure */
extern SL_Cstate     base_state;                     /* cartesian state of base coordinate system */
extern SL_quat       base_orient;                    /* cartesian orientation of base coordinate system */
extern Matrix        Alink_des[N_LINKS+1];           /* homogeneous transformation matrices for all links */
extern Matrix        link_pos_des;                   /* desired cart. pos of links */
extern Matrix        joint_cog_mpos_des;             /* vector of mass*COG of each joint based on desireds*/
extern Matrix    joint_origin_pos_des;           /* vector of pos. of local joint coord.sys based on des.*/
extern Matrix        joint_axis_pos_des;             /* unit vector of joint rotation axis based on des.*/
extern SL_DJstate    joint_default_state[N_DOFS+1];
extern SL_OJstate    joint_opt_state[N_DOFS+1];
extern SL_Jstate     joint_state[N_DOFS+1];
extern SL_DJstate    joint_des_state[N_DOFS+1];
extern double        joint_range[N_DOFS+1][3+1];
extern SL_link       links[N_DOFS+1];                /* specs of links: mass, inertia, cm */


// FROM MDEFS.H FILE
#define Power(x, y)	(pow((double)(x), (double)(y)))
#define Sqrt(x)		(sqrt((double)(x)))

#define Abs(x)		(fabs((double)(x)))

#define Exp(x)		(exp((double)(x)))
#define Log(x)		(log((double)(x)))

#define Sin(x)		(sin((double)(x)))
#define Cos(x)		(cos((double)(x)))
#define Tan(x)		(tan((double)(x)))

#define ArcSin(x)       (asin((double)(x)))
#define ArcCos(x)       (acos((double)(x)))
#define ArcTan(x)       (atan((double)(x)))

#define Sinh(x)          (sinh((double)(x)))
#define Cosh(x)          (cosh((double)(x)))
#define Tanh(x)          (tanh((double)(x)))


#ifndef E
#define E		2.71828182845904523536029
#endif
#ifndef Pi
#define Pi		3.14159265358979323846264
#endif
#define Degree		0.01745329251994329576924

// kinematics functions from SL
void kinematics(SL_DJstate *state,SL_Cstate *basec, SL_quat *baseo, SL_endeff *eff,
		   double **Xmcog, double **Xaxis, double **Xorigin, double **Xlink, double ***Ahmat);
void jacobian(Matrix lp, Matrix jop, Matrix jap, Matrix Jac);

// loading joint limits from SL config files
void load_joint_limits();
int read_sensor_offsets(char *fname);

// SL kinematics functions copied for convenience
void revoluteGJacColumn(Vector p, Vector pi, Vector zi, Vector c);
void setDefaultEndeffector(void);

#endif /* KINEMATICS_H_ */
