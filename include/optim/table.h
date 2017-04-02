/*
 * table.h
 *
 *  Created on: Feb 1, 2017
 *      Author: okoc
 */

#ifndef TABLE_H_
#define TABLE_H_

/* Table Variables */
// where is the origin? dist to tables edge on the robot side
static const double dist_to_table = -1.15; // -0.8; //-3.50; //-3.34;//-3.6;//-3.34//-2.74 // 0.86
static const double table_height = -0.76; // this changes +/- 1 cm across the table, measured from table top
static const double table_length = 2.76; //more like 273.8
static const double net_height   = 0.144; // 0.1525;
static const double net_overhang = 0.1525; // seems inaccurate
static const double net_thickness = 0.01;
static const double table_width  = 1.525;
static const double table_thickness = 0.056; // do we need this?
static const double net_restitution = 0.05;
static const double table_center = 0.0;

/* Table Tennis Ball Variables */
static const double ball_radius  = 0.02;
static const double ball_mass    = 0.0027;
static const double ball_contact_damping = 0; // filled in SimBall
static const double ball_contact_spring = 0;  // filled in SimBall

/* Table Tennis Racket Radius */
static const double racket_radius = 0.076; // shorter axis about 15.2 cm, longer approx 15.5 - 15.6

/* Stand Variables */
static const double stand_height  = -1.16; //-0.95;
static const double radius_bottom = 0.1;
static const double radius_top    = 0.02;
static const double stand_x       = 0.02;  //-0.85
static const double stand_y       = -0.59;
static const double stand_z	    = 0.9;

/* Floor */
static const double floor_level = -1.71;

#endif /* TABLE_H_ */
