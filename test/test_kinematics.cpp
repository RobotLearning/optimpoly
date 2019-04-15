/*
 * kinematics.cpp
 *
 * Testing the kinematics inside the optimization
 *
 *  Created on: Feb 20, 2017
 *      Author: okoc
 */

#include "kinematics.h"
#include "kinematics.hpp"
#include "optim.h"
#include "player.hpp"
#include "gtest/gtest.h"
#include <armadillo>
#include <iostream>

using namespace arma;
using namespace optim;

#define OPTIM_DIM 2 * NDOF + 1

TEST(KinematicsTests, CompareRacketStateGivenQAndQdWithMatlab0) {
  // Comparing both C and C++ versions of calc_racket_state function

  static double q[NDOF], qdot[NDOF], pos[NCART], vel[NCART], normal[NCART];
  for (int i = 0; i < NDOF; i++) {
    q[i] = 1.0;
    qdot[i] = 0.0;
  }
  vec7 q0_cpp = ones<vec>(NDOF);
  joint q_cpp;
  q_cpp.q = q0_cpp;
  player::racket robot_racket;

  // C version
  calc_racket_state(q, qdot, pos, vel, normal);
  // C++ version
  calc_racket_state(q_cpp, robot_racket);

  vec3 x_c(pos);
  vec3 xdot_c(vel);
  vec3 n_c(normal);

  vec3 x_MATLAB = {-0.8309, 0.0861, -0.6672};
  vec3 xdot_MATLAB = zeros<vec>(3);
  vec3 n_MATLAB = {0.3857, 0.0174, -0.9225};

  EXPECT_TRUE(approx_equal(x_c, x_MATLAB, "absdiff", 0.002));
  EXPECT_TRUE(approx_equal(xdot_c, xdot_MATLAB, "absdiff", 0.002));
  EXPECT_TRUE(approx_equal(n_c, n_MATLAB, "absdiff", 0.002));
  EXPECT_TRUE(approx_equal(robot_racket.pos, x_c, "absdiff", 0.002));
  EXPECT_TRUE(approx_equal(robot_racket.vel, xdot_c, "absdiff", 0.002));
  EXPECT_TRUE(approx_equal(robot_racket.normal, n_c, "absdiff", 0.002));
}
