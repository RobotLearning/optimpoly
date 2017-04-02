/*
 * kinematics.c
 *
 * Here we include the kinematics related functions taken from SL
 *
 *  Created on: Jun 22, 2016
 *      Author: okoc
 */

#include "../../../player/old/ccode/kinematics.h"

#include "stdio.h"
#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "../../../player/old/ccode/SL0.h"

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

 TODO: replace the header files! Does it become faster?

 ******************************************************************************/
void
kinematics(SL_DJstate *state,SL_Cstate *basec, SL_quat *baseo, SL_endeff *eff,
		   double **Xmcog, double **Xaxis, double **Xorigin, double **Xlink, double ***Ahmat) {

	#include "LInfo_declare.h"
	#include "LInfo_math.h"

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
  int i,j;

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
 * Copied from SL_common. Dont want to call the function from SL because
 * we have to include a lot of extra SL files
 *
 */
int read_sensor_offsets(char *fname) {

  int j,i,rc;
  char   string[100];
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

  for (i=1; i<= N_DOFS; ++i) {
    if (!find_keyword(in, &(joint_names[i][0]))) {
      printf("ERROR: Cannot find offset for >%s<!\n",joint_names[i]);
      fclose(in);
      return FALSE;
    }
    rc=fscanf(in,"%lf %lf %lf %lf %lf %lf",
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

/*!*****************************************************************************
 *******************************************************************************
 \note  find_keyword
 \date  October, 1995

 \remarks

 find a key word in the file given by fp and positions the read cursor
 right after the keyword

 *******************************************************************************
 Function Parameters: [in]=input,[out]=output

 \param[in]     fp       : file pointer
 \param[in]     word     : keyword to be found

 ******************************************************************************/
int
find_keyword(FILE *fp, char *name)

{

	int  i,j,c;
	int  rc = TRUE;
	char string[strlen(name)*2];
	int  l;
	char sep[]={' ','\n',':',',',';','=','\t','\0'};

	rewind(fp);
	l  = strlen(name);
	i  = 0;

	while (rc != EOF) {

		rc=fgetc(fp);
		if ( rc != EOF ) {

			string[i++] = rc;
			string[i]   = '\0';

			if ( strstr(string,name) != NULL) {
				// wait for one more character to judge whether this string
				// has the correct end delimiters
				if (strchr(sep,string[i-1]) != NULL) {
					// now check for preceeding delimiter

					if (i-l-2 < 0) // this means "name" was the first string in file
						return TRUE;
					else if (strchr(sep,string[i-l-2]) != NULL) //otherwise check delim
						return TRUE;
				}
			}

			if (i >= 2*l-1) {
				strcpy(string,&(string[i-l]));
				i = strlen(string);
			}

		}

	}

	return FALSE;

}

