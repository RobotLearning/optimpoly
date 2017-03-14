/*
 * constants.h
 *
 *  Created on: Feb 2, 2017
 *      Author: okoc
 */

#ifndef CONSTANTS_H_
#define CONSTANTS_H_

// indices for simplifying code comprehension
#define X 0
#define Y 1
#define Z 2
#define W 3 // for quaternion
#define DX 3
#define DY 4
#define DZ 5
#define NCART 3
#define NDOF 7

static const double VHPY = -0.3; // in front of the robot
static const double dt = 0.002; // 500 Hz robot operation

/************************************* constants  ***********************/

/* Contact Coefficients */
/* coefficient of restitution for the table (i.e. rebound z-velocity multiplier) */
static const double CRT = 0.88;//0.91;cannon//0.8; human //0.88 das hatten wir mal experimentell bestimmt, 0.92 fittet aber besser in die beobachteten Daten
/* coefficient of table contact model on Y-direction */
static const double CFTY = 0.72;//0.89;cannon//0.78; human //0.6 is eigentlich das richtige, aber da immer ein wenig spin drin sein wird, passt das am besten in die beobachteten Daten
/* coefficient of table contact model on X-direction */
static const double CFTX = 0.68; //0.74
/* coefficent of restitution for racket */
static const double CRR = 0.78;//0.9;cannon//0.78;human//0.78;

/* Air drag coefficient */
static const double Cdrag = 0.1414;

/* for simulating different gravities */
static const double gravity = -9.802;

static const double Clift = 0.001; // coefficient of lift for the magnus force


#endif /* INCLUDE_CONSTANTS_H_ */
