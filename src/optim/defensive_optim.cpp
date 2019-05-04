/**
 * @file lazyoptim.cpp
 *
 * @brief NLOPT polynomial optimization functions for DEFENSIVE PLAYER are
 * stored here.
 *
 *
 *  Created on: Sep 11, 2016
 *      Author: okan
 */

#include "constants.h"
#include "kinematics.h"
#include "lookup.h"
#include "math.h"
#include "optim.h"
#include "stdlib.h"
#include "table.h"
#include "utils.h"
#include <armadillo>
#include <iostream>

#define INEQ_HIT_CONSTR_DIM 3
#define INEQ_LAND_CONSTR_DIM 8 // 11
#define INEQ_JOINT_CONSTR_DIM 2 * NDOF + 2 * NDOF

using namespace const_tt;

namespace optim {

/*
 * Calculates the cost function for table tennis Lazy Player (LP)
 * to find spline (3rd order strike+return) polynomials
 */
static double costfunc(unsigned n, const double *x_, double *grad,
                       void *my_func_params);

/*
 * This is the constraint that makes sure we LAND the ball
 */
static void land_ineq_constr(unsigned m, double *result, unsigned n,
                             const double *x_, double *grad,
                             void *my_func_params);

/*
 * This is the constraint that makes sure we only TOUCH the ball
 *
 */
static void hit_ineq_constr(unsigned m, double *result, unsigned n,
                            const double *x_, double *grad,
                            void *my_func_params);

/*
 * Update the incoming ball velocity with outgoing ball velocity using MIRROR
 * LAW
 *
 * The racket contact model in vector form is O = I + (1 + eps_R)*N*N'*(V - I)
 * where I is the incoming ball velocity
 *       N is the racket normal
 *       V is the racket velocity
 *       eps_R is the coefficient of restitution of the racket
 *
 * The outgoing ball velocity is post-multiplied by some constants \mu < 1
 * to account for ball drag
 *
 *
 */
static void racket_contact_model(const double *racketVel,
                                 const double *racketNormal,
                                 const std::vector<double> &mult_vel_,
                                 double *ballVel);

/*
 * First order hold to interpolate linearly at time T_
 * between racket pos,vel,normal entries
 *
 * IF T_ is nan, racket variables are assigned to zero-element of
 * relevant racket entries
 *
 */
static void interp_ball(const optim_des *params, const double T_,
                        double *ballpos, double *ballvel);

DefensiveOptim::DefensiveOptim(const vec7 &qrest, double lbIn[], double ubIn[],
                               bool land, bool lookup)
    : Optim(), mult_vel_(3, 0.0), land_(land), x_last_{0.0}, t_land_(-1.0),
      t_net_(-1.0), x_land_{0.0}, x_net_{0.0}, dist_b2r_norm_(1.0),
      dist_b2r_proj_(1.0), penalty_loc_(4, 0.0) {

  mult_vel_[0] = 0.9;
  mult_vel_[1] = 0.8;
  mult_vel_[2] = 0.83;
  penalty_loc_[1] = 0.23;
  penalty_loc_[3] = -3.22;
  lookup_ = lookup;
  player::load_lookup_table(lookup_table_);
  const_vec(NDOF, 1.0, w_.R_strike);
  w_.R_net = 1e1;
  w_.R_hit = 1e3;
  w_.R_land = 1e1;

  for (int i = 0; i < NDOF; i++) {
    qrest_[i] = qrest(i);
  }
  for (int i = 0; i < OPTIM_DIM_; i++) {
    ub_[i] = ubIn[i];
    lb_[i] = lbIn[i];
  }

  if (land_) {
    set_land_constr();
  } else {
    set_hit_constr();
  }
}

void DefensiveOptim::init_last_soln(double x[]) const {

  // initialize first dof entries to q0_
  for (int i = 0; i < NDOF; i++) {
    x[i] = qf_[i];
    x[i + NDOF] = qfdot_[i];
  }
  x[2 * NDOF] = T_;
  // cout << "Initialization from T_ = " << T_ << endl;
}

void DefensiveOptim::init_rest_soln(double x[]) const {

  // initialize first dof entries to q0_
  for (int i = 0; i < NDOF; i++) {
    x[i] = qrest_[i];
    x[i + NDOF] = 0.0;
  }
  x[2 * NDOF] = 0.5;
}

void DefensiveOptim::set_weights(const std::vector<double> &weights) {

  w_.R_net = weights[1];
  w_.R_hit = weights[0];
  w_.R_land = weights[2];
}

void DefensiveOptim::set_velocity_multipliers(const std::vector<double> &mult) {
  mult_vel_ = mult;
}

void DefensiveOptim::set_penalty_loc(const std::vector<double> &loc) {
  penalty_loc_ = loc;
}

void DefensiveOptim::set_land_constr() {

  double max_opt_time = 0.05;
  double tol_ineq_land[INEQ_LAND_CONSTR_DIM];
  double tol_ineq_joint[INEQ_JOINT_CONSTR_DIM];
  const_vec(INEQ_LAND_CONSTR_DIM, 1e-2, tol_ineq_land);
  const_vec(INEQ_JOINT_CONSTR_DIM, 1e-3, tol_ineq_joint);

  // LN = does not require gradients //
  opt_ = nlopt_create(NLOPT_AUGLAG, 2 * NDOF + 1);
  nlopt_opt local_opt = nlopt_create(NLOPT_LD_VAR2, 2 * NDOF + 1);
  nlopt_set_xtol_rel(local_opt, 1e-2);
  nlopt_set_lower_bounds(local_opt, lb_);
  nlopt_set_upper_bounds(local_opt, ub_);
  // nlopt_set_vector_storage(local_opt, 20);
  nlopt_set_local_optimizer(opt_, local_opt);
  nlopt_set_min_objective(opt_, costfunc, this);
  nlopt_add_inequality_mconstraint(opt_, INEQ_LAND_CONSTR_DIM, land_ineq_constr,
                                   this, tol_ineq_land);
  nlopt_add_inequality_mconstraint(opt_, INEQ_JOINT_CONSTR_DIM,
                                   joint_limits_ineq_constr, this,
                                   tol_ineq_joint);
  nlopt_set_lower_bounds(opt_, lb_);
  nlopt_set_upper_bounds(opt_, ub_);
  nlopt_set_xtol_rel(opt_, 1e-2);
  nlopt_set_maxtime(opt_, max_opt_time);
  nlopt_destroy(local_opt);
}

void DefensiveOptim::set_hit_constr() {

  double tol_ineq_hit[INEQ_HIT_CONSTR_DIM];
  double tol_ineq_joint[INEQ_JOINT_CONSTR_DIM];
  const_vec(INEQ_HIT_CONSTR_DIM, 1e-2, tol_ineq_hit);
  const_vec(INEQ_JOINT_CONSTR_DIM, 1e-3, tol_ineq_joint);

  // LN = does not require gradients //
  opt_ = nlopt_create(NLOPT_LD_SLSQP, 2 * NDOF + 1);
  nlopt_set_min_objective(opt_, costfunc, this);
  nlopt_add_inequality_mconstraint(opt_, INEQ_HIT_CONSTR_DIM, hit_ineq_constr,
                                   this, tol_ineq_hit);
  // nlopt_add_inequality_mconstraint(opt_, INEQ_JOINT_CONSTR_DIM,
  //		joint_limits_ineq_constr, this, tol_ineq_joint);
  nlopt_set_lower_bounds(opt_, lb_);
  nlopt_set_upper_bounds(opt_, ub_);
  nlopt_set_xtol_rel(opt_, 1e-2);
}

void DefensiveOptim::finalize_soln(const double x_[], double time_elapsed) {

  if (x_[2 * NDOF] > fmax(time_elapsed / 1e3, 0.05)) {
    // initialize first dof entries to q0_
    for (int i = 0; i < NDOF; i++) {
      qf_[i] = x_[i];
      qfdot_[i] = x_[i + NDOF];
    }
    T_ = x_[2 * NDOF];
    if (detach_)
      T_ -= (time_elapsed / 1e3);
    update_ = true;
  }
  // trigger_optim();
}

void DefensiveOptim::calc_times(
    const double x_[]) { // ball projected to racket plane

  static double qf_[NDOF];
  static double qfdot_[NDOF];
  static double vel[NCART];
  static double normal[NCART];
  static double pos[NCART];
  static double ballpos[NCART];
  static double ballvel[NCART];
  static double g = -9.8;
  static double net_y = dist_to_table - table_length / 2.0;
  static double table_z = floor_level - table_height;
  double discr = 0.0;
  double d;

  if (!vec_is_equal(OPTIM_DIM_, x_, x_last_)) {
    // extract state information from optimization variables
    for (int i = 0; i < NDOF; i++) {
      qf_[i] = x_[i];
      qfdot_[i] = x_[i + NDOF];
    }
    calc_racket_state(qf_, qfdot_, pos, vel, normal);
    interp_ball(param_des_, x_[2 * NDOF], ballpos, ballvel);
    // calculate deviation of ball to racket - hitting constraints
    calc_hit_distance(ballpos, pos, normal);
    racket_contact_model(vel, normal, mult_vel_, ballvel);
    t_net_ = (net_y - ballpos[Y]) / ballvel[Y];
    x_net_[Z] = ballpos[Z] + t_net_ * ballvel[Z] + 0.5 * g * t_net_ * t_net_;
    x_net_[X] = ballpos[X] + t_net_ * ballvel[X];

    // calculate ball planned landing
    d = ballpos[Z] - table_z;
    if (pow(ballvel[Z], 2) > 2 * g * d) {
      discr = sqrt(pow(ballvel[Z], 2) - 2 * g * d);
    }
    t_land_ = fmax((-ballvel[Z] - discr) / g, (-ballvel[Z] + discr) / g);
    t_net_ = (net_y - ballpos[Y]) / ballvel[Y];
    x_land_[X] = ballpos[X] + t_land_ * ballvel[X];
    x_land_[Y] = ballpos[Y] + t_land_ * ballvel[Y];

    make_equal(OPTIM_DIM_, x_, x_last_);
  }
}

double DefensiveOptim::calc_punishment() {

  double Jhit = 0;
  double Jland = 0;
  if (land_) {
    Jland = pow(x_net_[X] - penalty_loc_[0], 2) * w_.R_net +
            pow(x_net_[Z] - penalty_loc_[1], 2) * w_.R_net +
            pow(x_land_[X] - penalty_loc_[2], 2) * w_.R_land +
            pow(x_land_[Y] - penalty_loc_[3], 2) * w_.R_land;
  }
  Jhit = w_.R_hit * pow(dist_b2r_proj_, 2); // punish for hitting properly
  return Jhit + Jland;
}

void DefensiveOptim::calc_hit_distance(const double ball_pos[],
                                       const double racket_pos[],
                                       const double racket_normal[]) {

  double e[NCART];
  for (int i = 0; i < NCART; i++) {
    e[i] = ball_pos[i] - racket_pos[i];
  }
  dist_b2r_norm_ = inner_prod(NCART, racket_normal, e);
  dist_b2r_proj_ =
      sqrt(inner_prod(NCART, e, e) - dist_b2r_norm_ * dist_b2r_norm_);
}

double DefensiveOptim::test_soln(const double x_[]) const {

  double max_viol;
  // give info on constraint violation
  static double table_xmax = table_width / 2.0;
  static double table_ymax = dist_to_table - table_length;
  static double wall_z = 1.0;
  static double net_y = dist_to_table - table_length / 2.0;
  static double net_z = floor_level - table_height + net_height;
  static int count = 0;
  double *grad = 0;
  static double max_acc_violation; // at hitting time
  static double land_violation[INEQ_LAND_CONSTR_DIM];
  static double lim_violation[INEQ_JOINT_CONSTR_DIM]; // joint limit violations
                                                      // on strike and return
  joint_limits_ineq_constr(INEQ_JOINT_CONSTR_DIM, lim_violation, OPTIM_DIM_, x_,
                           grad, (void *)this);
  land_ineq_constr(INEQ_LAND_CONSTR_DIM, land_violation, OPTIM_DIM_, x_, grad,
                   (void *)this);
  double cost = costfunc(OPTIM_DIM_, x_, grad, (void *)this);

  if (verbose_) {
    // give info on solution vector
    printf("Optim count: %d\n", (++count));
    print_optim_vec(x_);
    printf("f = %.2f\n", cost);
    printf("Hitting constraints (b2r):\n");
    printf("Distance along normal: %.2f\n", dist_b2r_norm_);
    printf("Distance along racket: %.2f\n", dist_b2r_proj_);
    printf("Landing constraints:\n");
    printf("NetTime: %f\n", t_net_);
    printf("LandTime: %f\n", t_land_);
    printf("Below wall by : %f\n", wall_z - x_net_[Z]);
    printf("Above net by: %f\n", x_net_[Z] - net_z);
    printf("X: between table limits by [%f, %f]\n", table_xmax - x_land_[X],
           x_land_[X] + table_xmax);
    printf("Y: between table limits by [%f, %f]\n", net_y - x_land_[Y],
           x_land_[Y] - table_ymax);
    for (int i = 0; i < INEQ_JOINT_CONSTR_DIM; i++) {
      if (lim_violation[i] > 0.0) {
        printf("Joint limit violated by %.2f on joint %d\n", lim_violation[i],
               i % NDOF + 1);
      }
    }
  }
  max_acc_violation = calc_max_acc_violation(x_, q0_, q0dot_);
  max_viol =
      fmax(max_array(lim_violation, INEQ_JOINT_CONSTR_DIM), max_acc_violation);
  if (land_)
    max_viol = fmax(max_viol, max_array(land_violation, INEQ_LAND_CONSTR_DIM));
  return max_viol;
}

static double costfunc(unsigned n, const double *x_, double *grad,
                       void *my_func_params) {

  static double J1, Jland;
  static double a1[NDOF], a2[NDOF];
  double T_ = x_[2 * NDOF];

  DefensiveOptim *opt_ = (DefensiveOptim *)my_func_params;
  double *q0_ = opt_->q0_;
  double *q0dot_ = opt_->q0dot_;
  weights w_ = opt_->w_;

  // calculate the polynomial coeffs which are used in the cost calculation
  calc_strike_poly_coeff(q0_, q0dot_, x_, a1, a2);

  // calculate the landing time
  opt_->calc_times(x_);

  // calculate the punishments;
  Jland = opt_->calc_punishment();

  J1 = T_ * (3 * T_ * T_ * inner_w_prod(NDOF, w_.R_strike, a1, a1) +
             3 * T_ * inner_w_prod(NDOF, w_.R_strike, a1, a2) +
             inner_w_prod(NDOF, w_.R_strike, a2, a2));

  if (grad) {
    static double h = 1e-6;
    static double val_plus, val_minus;
    static double xx[2 * NDOF + 1];
    for (unsigned i = 0; i < n; i++)
      xx[i] = x_[i];
    for (unsigned i = 0; i < n; i++) {
      xx[i] += h;
      val_plus = costfunc(n, xx, NULL, my_func_params);
      xx[i] -= 2 * h;
      val_minus = costfunc(n, xx, NULL, my_func_params);
      grad[i] = (val_plus - val_minus) / (2 * h);
      xx[i] += h;
    }
  }

  // std::cout << J1 << "\t" << Jhit << "\t" << Jland << std::endl;
  return J1 + Jland;
}

static void land_ineq_constr(unsigned m, double *result, unsigned n,
                             const double *x_, double *grad,
                             void *my_func_params) {

  static double table_xmax = table_width / 2.0;
  // static double table_ymax = dist_to_table - table_length;
  static double wall_z = 1.0;
  // static double net_y = dist_to_table - table_length/2.0;
  static double net_z = floor_level - table_height + net_height;

  DefensiveOptim *opt_ = (DefensiveOptim *)my_func_params;
  opt_->calc_times(x_);

  result[0] = -opt_->dist_b2r_norm_;
  result[1] = opt_->dist_b2r_norm_ - ball_radius;
  result[2] = opt_->dist_b2r_proj_ - racket_radius;
  result[3] = opt_->x_net_[Z] - wall_z;
  result[4] = -opt_->x_net_[Z] + net_z;
  result[5] = opt_->x_net_[X] - table_xmax;
  result[6] = -opt_->x_net_[X] - table_xmax;
  result[7] = -opt_->t_net_;

  if (grad) {
    static double h = 1e-6;
    static double res_plus[INEQ_LAND_CONSTR_DIM],
        res_minus[INEQ_LAND_CONSTR_DIM];
    static double xx[2 * NDOF + 1];
    for (unsigned i = 0; i < n; i++)
      xx[i] = x_[i];
    for (unsigned i = 0; i < n; i++) {
      xx[i] += h;
      land_ineq_constr(m, res_plus, n, xx, NULL, my_func_params);
      xx[i] -= 2 * h;
      land_ineq_constr(m, res_minus, n, xx, NULL, my_func_params);
      xx[i] += h;
      for (unsigned j = 0; j < m; j++)
        grad[j * n + i] = (res_plus[j] - res_minus[j]) / (2 * h);
    }
  }
}

static void hit_ineq_constr(unsigned m, double *result, unsigned n,
                            const double *x_, double *grad,
                            void *my_func_params) {

  static double qf_[NDOF];
  static double qfdot_[NDOF];
  static double vel[NCART];
  static double normal[NCART];
  static double pos[NCART];
  static double ballpos[NCART];
  static double ballvel[NCART];

  DefensiveOptim *opt_ = (DefensiveOptim *)my_func_params;

  for (unsigned i = 0; i < NDOF; i++) {
    qf_[i] = x_[i];
    qfdot_[i] = x_[i + NDOF];
  }
  interp_ball(opt_->param_des_, x_[2 * NDOF], ballpos, ballvel);
  calc_racket_state(qf_, qfdot_, pos, vel, normal);
  // calculate deviation of ball to racket - hitting constraints
  opt_->calc_hit_distance(ballpos, pos, normal);

  if (grad) {
    static double h = 1e-6;
    static double res_plus[INEQ_HIT_CONSTR_DIM], res_minus[INEQ_HIT_CONSTR_DIM];
    static double xx[2 * NDOF + 1];
    for (unsigned i = 0; i < n; i++)
      xx[i] = x_[i];
    for (unsigned i = 0; i < n; i++) {
      xx[i] += h;
      hit_ineq_constr(m, res_plus, n, xx, NULL, my_func_params);
      xx[i] -= 2 * h;
      hit_ineq_constr(m, res_minus, n, xx, NULL, my_func_params);
      xx[i] += h;
      for (unsigned j = 0; j < m; j++)
        grad[j * n + i] = (res_plus[j] - res_minus[j]) / (2 * h);
    }
  }

  result[0] = -opt_->dist_b2r_norm_;
  result[1] = opt_->dist_b2r_norm_ - ball_radius;
  result[2] = opt_->dist_b2r_proj_ - racket_radius;
}

static void racket_contact_model(const double *racketVel,
                                 const double *racketNormal,
                                 const std::vector<double> &mult_vel_,
                                 double *ballVel) {

  static const double racket_param = 0.78;
  static double diffVel[NCART];
  static double normalMultSpeed[NCART];
  double speed;

  for (int i = 0; i < NCART; i++)
    diffVel[i] = racketVel[i] - ballVel[i];

  speed = (1 + racket_param) * inner_prod(NCART, racketNormal, diffVel);

  for (int i = 0; i < NCART; i++) {
    normalMultSpeed[i] = speed * racketNormal[i];
  }
  vec_plus(NCART, normalMultSpeed, ballVel);
  for (int i = 0; i < NCART; i++) {
    ballVel[i] *= mult_vel_[i];
  }
}

static void interp_ball(const optim_des *data, const double T_, double *ballpos,
                        double *ballvel) {

  const double dt = data->dt;
  if (std::isnan(T_)) {
    printf("Warning: T_ value is nan!\n");
    for (int i = 0; i < NCART; i++) {
      ballpos[i] = data->ball_pos(i, 0);
      ballvel[i] = data->ball_vel(i, 0);
    }
  } else {
    const unsigned Nmax = data->Nmax;
    unsigned N = (int)(T_ / dt);
    double Tdiff = T_ - N * dt;
    for (int i = 0; i < NCART; i++) {
      if (N < Nmax - 1) {
        ballpos[i] =
            data->ball_pos(i, N) +
            (Tdiff / dt) * (data->ball_pos(i, N + 1) - data->ball_pos(i, N));
        ballvel[i] =
            data->ball_vel(i, N) +
            (Tdiff / dt) * (data->ball_vel(i, N + 1) - data->ball_vel(i, N));
      } else {
        ballpos[i] = data->ball_pos(i, Nmax - 1);
        ballvel[i] = data->ball_vel(i, Nmax - 1);
      }
    }
  }
}

} // namespace optim
