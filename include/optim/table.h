/**
 * @file table.h
 *
 * @brief Constants related to the Table Tennis Table.
 *
 * All table contants related to table tennis simulation and ball prediction
 * are stored here.
 *
 *  Created on: Feb 1, 2017
 *      Author: okoc
 */

#pragma once
/* Table Variables */
// where is the origin? dist to tables edge on the robot side
static const double dist_to_table = -1.15;
static const double table_height =
    -0.76; //! this changes +/- 1 cm across the table, measured from table top
static const double table_length = 2.76;   //! more like 273.8
static const double net_height = 0.144;    // 0.1525;
static const double net_overhang = 0.1525; // seems inaccurate
static const double net_thickness = 0.01;
static const double table_width = 1.525;
static const double table_thickness = 0.056; // do we need this?
static const double net_restitution = 0.05;
static const double table_center = 0.0;

/* Stand Variables */
static const double stand_height = -1.16; //-0.95;
static const double radius_bottom = 0.1;
static const double radius_top = 0.02;
static const double stand_x = 0.02; //-0.85
static const double stand_y = -0.59;
static const double stand_z = 0.9;

/* Floor */
static const double floor_level = -1.67;
