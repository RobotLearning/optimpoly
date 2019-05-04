/**
 * @file vhp_optim.cpp
 *
 * @brief Fixing a virtual hitting plane to find qf_, qfdot_ joint angles
 * and velocities.
 *
 *  Created on: Mar 5, 2017
 *      Author: okoc
 */
#include "constants.h"
#include "kinematics.h"
#include "math.h"
#include "optim.h"
#include "stdlib.h"
#include "utils.h"
#include <armadillo>
#include <assert.h>

namespace optim {

/*
 * Penalize (unweighted) squared distance (in joint space) to joint limit
 * averages
 *
 * Adds also joint velocity penalties
 *
 */
static double penalize_dist_to_limits(unsigned n, const double *x, double *grad,
                                      void *my_func_params);

/*
 * Constant function for simple inverse kinematics
 */
// static double const_costfunc(unsigned n,
//                                const double *x,
//                                double *grad,
//                               void *my_func_params);

/*
 * This is the constraint that makes sure we hit the ball
 */
static void kinematics_eq_constr(unsigned m, double *result, unsigned n,
                                 const double *x, double *grad,
                                 void *my_function_data);

HittingPlane::HittingPlane(const vec7 &qrest, double lb[], double ub[]) {

  double tol_eq[EQ_CONSTR_DIM];
  const_vec(EQ_CONSTR_DIM, 1e-2, tol_eq);
  // set tolerances equal to second argument

  // LN = does not require gradients //
  opt_ = nlopt_create(NLOPT_LN_COBYLA, OPTIM_DIM_);
  nlopt_set_xtol_rel(opt_, 1e-2);
  nlopt_set_lower_bounds(opt_, lb);
  nlopt_set_upper_bounds(opt_, ub);
  nlopt_set_min_objective(opt_, penalize_dist_to_limits, this);
  nlopt_add_equality_mconstraint(opt_, EQ_CONSTR_DIM, kinematics_eq_constr, this,
                                 tol_eq);

  for (int i = 0; i < NDOF; i++) {
    qrest_[i] = qrest(i);
    ub_[i] = ub[i];
    lb_[i] = lb[i];
    limit_avg[i] = (ub[i] + lb[i]) / 2.0;
  }
}

void HittingPlane::fix_hitting_time(double time_pred) {
  if (time_pred > 0.05)
    T_ = time_pred;
}

void HittingPlane::init_last_soln(double x[2 * NDOF]) const {

  // initialize first dof entries to q0
  for (int i = 0; i < NDOF; i++) {
    x[i] = qf_[i];
    x[i + NDOF] = qfdot_[i];
  }
}

void HittingPlane::init_rest_soln(double x[2 * NDOF]) const {

  // initialize first dof entries to q0
  for (int i = 0; i < NDOF; i++) {
    x[i] = qrest_[i];
    x[i + NDOF] = 0.0;
  }
}

void HittingPlane::finalize_soln(const double x[2 * NDOF],
                                 double time_elapsed) {

  if (T_ > 0.05) {
    // initialize first dof entries to q0
    for (int i = 0; i < NDOF; i++) {
      qf_[i] = x[i];
      qfdot_[i] = x[i + NDOF];
    }
    if (detach_) {
      T_ -= (time_elapsed / 1e3);
    }
    update_ = true;
  }
}

double HittingPlane::test_soln(const double x[]) const {

  double x_[2 * NDOF + 1];
  for (int i = 0; i < 2 * NDOF; i++)
    x_[i] = x[i];
  x_[2 * NDOF] = T_;

  // give info on constraint violation
  double *grad = 0;
  double kin_violation[EQ_CONSTR_DIM];
  double lim_violation[INEQ_CONSTR_DIM]; // joint limit violations on strike and
                                         // return
  kinematics_eq_constr(EQ_CONSTR_DIM, kin_violation, 2 * NDOF, x, grad,
                       (void *)this);
  joint_limits_ineq_constr(INEQ_CONSTR_DIM, lim_violation, 2 * NDOF, x_, grad,
                           (void *)this);
  // double cost = costfunc(OPTIM_DIM_, x, grad, coparams);

  if (verbose_) {
    // give info on solution vector
    print_optim_vec(x_);
    // printf("f = %.2f\n",cost);
    printf("Position constraint violation: [%.2f %.2f %.2f]\n",
           kin_violation[0], kin_violation[1], kin_violation[2]);
    printf("Velocity constraint violation: [%.2f %.2f %.2f]\n",
           kin_violation[3], kin_violation[4], kin_violation[5]);
    printf("Normal constraint violation: [%.2f %.2f %.2f]\n", kin_violation[6],
           kin_violation[7], kin_violation[8]);
    for (int i = 0; i < INEQ_CONSTR_DIM; i++) {
      if (lim_violation[i] > 0.0)
        printf("Joint limit violated by %.2f on joint %d\n", lim_violation[i],
               i % NDOF + 1);
    }
  }

  return fmax(max_abs_array(kin_violation, EQ_CONSTR_DIM),
              max_array(lim_violation, INEQ_CONSTR_DIM));
}

static double penalize_dist_to_limits(unsigned n, const double *x, double *grad,
                                      void *my_func_params) {

  double cost = 0.0;
  HittingPlane *vhp = (HittingPlane *)my_func_params;

  if (grad) {
    assert(n == 2 * NDOF); // to turn off unused-parameter warning
    for (int i = 0; i < NDOF; i++) {
      grad[i] = 2 * (x[i] - vhp->limit_avg[i]);
      grad[i + NDOF] = 2 * x[i + NDOF];
    }
  }

  for (int i = 0; i < NDOF; i++) {
    cost += pow(x[i] - vhp->limit_avg[i], 2);
    cost += pow(x[i + NDOF], 2);
  }
  return cost;
}

// static double const_costfunc(unsigned n,
//                             const double *x,
//                             double *grad,
//                             void *my_func_params) {
//
//	return 1.0;
//}

static void kinematics_eq_constr(unsigned m, double *result, unsigned n,
                                 const double *x, double *grad,
                                 void *my_function_data) {

  static double pos[NCART];
  static double qfdot_[NDOF];
  static double vel[NCART];
  static double normal[NCART];
  static double qf_[NDOF];

  HittingPlane *vhp = (HittingPlane *)my_function_data;

  // extract state information from optimization variables
  for (int i = 0; i < NDOF; i++) {
    qf_[i] = x[i];
    qfdot_[i] = x[i + NDOF];
  }

  // compute the actual racket pos,vel and normal
  calc_racket_state(qf_, qfdot_, pos, vel, normal);

  // deviations from the desired racket frame
  for (int i = 0; i < NCART; i++) {
    result[i] = pos[i] - vhp->param_des_->racket_pos(i);
    result[i + NCART] = vel[i] - vhp->param_des_->racket_vel(i);
    result[i + 2 * NCART] = normal[i] - vhp->param_des_->racket_normal(i);
  }

  if (grad) { // to turn off unused-param warning
    assert(m == 3 * NCART);
    assert(n == 2 * NDOF);
    // TODO
  }
}

} // namespace optim
