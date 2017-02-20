/*
 * Subset of SL structs needed
 *
 * SL0.h
 *
 *  Created on: Feb 6, 2017
 *      Author: okoc
 */

#ifndef SL_H_
#define SL_H_

/*! defines for the preference and config files */
#define CONFIG   "config/"
#define PREFS    "prefs/"

#define TRUE 1
#define FALSE 0
#define START_INDEX 1
#define N_CART 3
#define N_QUAT 4
#define N_ENDEFFS 1
#define N_LINKS 6
#define _X_ (0+START_INDEX)
#define _Y_ (1+START_INDEX)
#define _Z_ (2+START_INDEX)
#define _Q0_ 1
#define _Q1_ 2
#define _Q2_ 3
#define _Q3_ 4
#define N_DOFS 7

/*! defines that are used to parse the config and prefs files */
#define MIN_THETA    1
#define MAX_THETA    2
#define THETA_OFFSET 3

// dimensions of the robot
#define ZSFE 0.346              //!< z height of SAA axis above ground
#define ZHR  0.505              //!< length of upper arm until 4.5cm before elbow link
#define YEB  0.045              //!< elbow y offset
#define ZEB  0.045              //!< elbow z offset
#define YWR -0.045              //!< elbow y offset (back to forewarm)
#define ZWR  0.045              //!< elbow z offset (back to forearm)
#define ZWFE 0.255              //!< forearm length (minus 4.5cm)

typedef double*  Vector;
typedef double** Matrix;

typedef struct { /*!< Link parameters */
  double   m;             /*!< Mass */
  double   mcm[N_CART+1]; /*!< Center of mass multiplied with the mass */
  double   inertia[N_CART+1][N_CART+1];  /*!< Moment of inertia */
  double   vis;           /*!< viscous friction term */
  double   coul;          /*!< coulomb friction */
  double   stiff;         /*!< spring stiffness */
  double   cons;          /*!< constant term */
} SL_link;

typedef struct { /*!< joint space state for each DOF */
  double   th;   /*!< theta */
  double   thd;  /*!< theta-dot */
  double   thdd; /*!< theta-dot-dot */
  double   ufb;  /*!< feedback portion of command */
  double   u;    /*!< torque command */
  double   load; /*!< sensed torque */
} SL_Jstate;

typedef struct { /*!< desired values for controller */
  double   th;   /*!< theta */
  double   thd;  /*!< theta-dot */
  double   thdd; /*!< theta-dot-dot */
  double   uff;  /*!< feedforward torque command */
  double   uex;  /*!< externally imposed torque */
} SL_DJstate;

typedef struct { /*!< Cartesian state */
  double   x[N_CART+1];    /*!< Position [x,y,z] */
  double   xd[N_CART+1];   /*!< Velocity */
  double   xdd[N_CART+1];  /*!< Acceleration */
} SL_Cstate;

typedef struct { /*!< Quaternion orientation */
  double   q[N_QUAT+1];    /*!< Position [q0,q1,q2,q3] */
  double   qd[N_QUAT+1];   /*!< Velocity */
  double   qdd[N_QUAT+1];  /*!< Acceleration */
  double   ad[N_CART+1];   /*!< Angular Velocity [alpha,beta,gamma] */
  double   add[N_CART+1];  /*!< Angular Acceleration */
} SL_quat;

typedef struct { /*!< desired state for optimization */
  double   th;   /*!< desired theta */
  double   w;    /*!< feedforward command */
} SL_OJstate;

typedef struct { /*!< end-effector parameters */
  double   m;             /*!< Mass */
  double   mcm[N_CART+1]; /*!< mass times Center of mass */
  double   x[N_CART+1];   /*!< end position of endeffector in local coordinates*/
  double   a[N_CART+1];   /*!< orientation of the tool in Euler Angles (x-y-z) */
  int      c[2*N_CART+1]; /*!< constraint status of the endeffector */
  double   cf[N_CART+1];  /*!< contact force in world coordinates at the endeffector */
  double   ct[N_CART+1];  /*!< contact torques in world coordinates at the endeffector */
} SL_endeff;

typedef struct { /*!< Vision Blob */
  char       status;
  SL_Cstate  blob;
} SL_VisionBlob;

#endif /* SL0_H_ */
