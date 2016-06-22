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

#include <math.h>
#include "SL.h"
#include "SL_user.h"

#define DOF 7
#define CART 3

// kinematics functions from SL
void kinematics(SL_DJstate *state,SL_Cstate *basec, SL_quat *baseo, SL_endeff *eff,
		   double **Xmcog, double **Xaxis, double **Xorigin, double **Xlink, double ***Ahmat);
void jacobian(Matrix lp, Matrix jop, Matrix jap, Matrix Jac);

// loading joint limits from SL config files
int read_sensor_offsets(char *fname);

// SL kinematics functions copied for convenience
void revoluteGJacColumn(Vector p, Vector pi, Vector zi, Vector c);
void setDefaultEndeffector(void);


#endif /* KINEMATICS_H_ */
