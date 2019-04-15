/**
 * @file tabletennis.h
 *
 * @brief Header for table tennis ball prediction functions.
 *
 * Mostly taken from table_tennis_common.h file.
 *
 *  Created on: Feb 1, 2017
 *      Author: okoc
 */
#pragma once
#include "constants.h"
#include "table.h"
#include <armadillo>

using arma::mat;
using arma::vec;
using arma::vec3;
using arma::vec6;
using arma::zeros;

namespace player {

/**
 * @brief Racket positions, velocities and normal
 * used to predict racket contact.
 */
struct racket {
  vec3 pos = zeros<vec>(3);
  vec3 vel = zeros<vec>(3);
  vec3 normal = zeros<vec>(3);
};

/**
 * @brief Information about the status of the game.
 *
 * This information is useful in SIMULATION. For instance
 * when we're running tests for performance of optimizers, these
 * are useful to collect statistics.
 */
struct status {
  bool hit = false;         //!< robot hit the incoming ball
  bool has_bounced = false; //!< the ball has bounced on robot court (without
                            //!< checking for legal bounce)
  bool legal_bounce =
      false; //!< the ball has bounced only once on the robot court
  bool has_landed = false;  //!< the ball has landed on the other side (without
                            //!< checking for legal land)
  bool legal_serve = false; //!< the ball was hit, landed first on robot court
                            //!< and then on opponent's court
  bool legal_land = false;  //!< ball has landed legally on the other side
                            //!< (bounced once only on robot court)
  bool touched_ground =
      false; //!< the ball touched the ground level (vel. zeroed)
};

/**
 * @brief Ball parameters used to predict future ball path
 * and to calculate desired racket parameters.
 */
struct ball_params {

  /* Contact Coefficients */
  double CRT = 0.88; //!< coefficient of restitution for the table (i.e. rebound
                     //!< z-velocity multiplier)
  double CFTY = 0.72; //!< coefficient of table contact model on Y-direction
  double CFTX = 0.68; //!< coefficient of table contact model on X-direction
  double CRR = 0.78;  //!< coefficent of restitution for racket

  double Cdrag = 0.1414;       //!< Air drag coefficient
  double gravity = -9.802;     //!< gravity
  double Clift = 0.001;        //!< coefficient of lift for the magnus force
  double mu = 0.10;            //!< dynamic coefficient of friction
  double init_topspin = -50.0; //!< initial topspin amount
};

/**
 * @brief Table Tennis ball prediction methods
 *
 * Spin and spinless ball flight models can be used the predict the
 * next ball state. Checking for contact, moreover, enables to predict the
 * future of the table tennis trial, given racket and the table.
 */
class TableTennis {

private:
  bool CHECK_CONTACTS = true; // turn off for simplified debugging
  bool SPIN_MODE = false;     // turn on prediction with a spin model
  bool VERBOSE = false;
  status stats;
  ball_params params; // ball prediction parameters
  vec3 ball_pos = zeros<vec>(3);
  vec3 ball_vel = zeros<vec>(3);
  vec3 ball_spin =
      zeros<vec>(3); // ball angular velocity = 0 if spin mode is turned OFF

  /** @brief Initialize constant angular velocity (a.k.a. spin) for the spinning
   * ball. */
  void init_topspin(const double val = -50);

  /**
   * @brief Spinning nonlinear ball flight model.
   *
   * Flight model including including airdrag and Magnus force (for spin).
   *
   * @return Ball accelerations as a 3-vector.
   *
   */
  vec3 flight_model() const;

  /**
   * @brief Spin free nonlinear flight model including airdrag
   *
   * Cdrag is a parameter used to add drag force.
   *
   * @return Ball accelerations as 3-vector.
   */
  vec3 drag_flight_model() const;

  /**
   * @brief Simple contact model that uses spin if spin mode is turned on.
   *
   * Linear contact model that updates the outgoing rebound velocities only.
   * Checks if the contact is of roll or slide type. Based on a spin model
   * from a Japanese table tennis paper.
   *
   * Note: Spin is not changed!
   * Coeff of restitution and friction used.
   *
   */
  vec3 table_contact_model(const vec3 &ball_vel_in) const;

  /**
   * @brief FIRST ORDER Symplectic Euler integration for dt seconds.
   *
   * First integrating the accelerations to velocities by dt
   * Then integrating the velocities to positions by dt
   * These (pos and vel) are kept in the ball candidate vector
   *
   * We dont directly update the ball states, this is done in
   * integrate_ball_state. This way we can check for contacts
   * in between.
   *
   * @param dt Prediction horizon.
   * @param ball_next_pos Next position using calculated ball accelerations.
   * @param ball_next_vel Next velocity using calculated ball accelerations.
   */
  void symplectic_euler(const double dt, vec3 &ball_cand_pos,
                        vec3 &ball_cand_vel) const;

  /**
   * @brief Checks if a contact will occur.
   *
   * Checking against table, net, racket, ground, and possibly
   * a simulated human opponent!
   *
   * @param robot_racket Racket centre positions,velocities and normal of the
   * robot
   * @param ball_cand_pos Balls next candidate positions (after symplectic
   * int.).
   * @param ball_cand_vel Balls next candidate vels. (after symplectic int.).
   */
  void check_contact(const racket &robot_racket, vec3 &ball_cand_pos,
                     vec3 &ball_cand_vel); // calls the contact functions below

  /**
   * @brief Checks for legal bounce on robot court
   * and legal landing on the opponents court
   *
   * If the ball bounced only once before being hit on the robot court
   * then it is a LEGAL_BOUNCE.
   *
   * If there was already a hit and ball hasn't bounced before,
   * then the bounce location is checked and if it is on the opponent's court
   * and if LEGAL_BOUNCE is TRUE then it is a LAND.
   */
  void check_legal_bounce(const vec3 &ball_cand_pos, const vec3 &ball_cand_vel);

  /**
   * @brief Checks for legal LAND on opponents court.
   *
   * Called when ball's z-location is BELOW the table.
   * If there was already a hit and ball hasn't LANDed before
   * then if bounce location is on the opponent's court it is a LEGAL LAND.
   *
   */
  void check_legal_land(const vec3 &ball_cand_pos, const vec3 &ball_cand_vel);

  /**
   * @brief Condition to determine if ball hits the table.
   *
   * Useful for prediction including a rebound model
   * Useful also in KF/EKF filtering.
   * If contact is detected, then a table contact model (with constant
   * parameters) will update the next candidate ball velocities.
   *
   * @param ball_cand_pos Next candidate ball positions. Used to check contact
   * only.
   * @param ball_cand_vel Next candidate ball velocities. Updated if contact
   * happens.
   */
  void check_ball_table_contact(const vec3 &ball_cand_pos, vec3 &ball_cand_vel);

  /**
   * @brief Checks contact with net.
   *
   * Curious way to check contact with net:
   * if the net's distance to integrated y-state and distance to current y-state
   * signs do not match, it means that the ball is in contact with the net.
   * Then the ball candidate velocities are updated according to a (simplistic)
   * net contact model.
   *
   * @param ball_cand_pos Next candidate ball positions. Updated if contact
   * happens.
   * @param ball_cand_vel Next candidate ball velocities. Updated if contact
   * happens.
   */
  void check_ball_net_contact(vec3 &ball_cand_pos, vec3 &ball_cand_vel) const;

  /**
   * @brief  Checks contact with racket.
   *
   * If there is contact, then the predicted state will be transformed according
   * to a racket-ball contact model. If racket is going to hit the ball, hit
   * data member is set to TRUE so that we do not hit it the next time.
   *
   * @param robot_racket Racket center pos,vel and normals of the robot
   * @param ball_cand_pos Next candidate ball positions. Used to check contact
   * only.
   * @param ball_cand_vel Next candidate ball velocities. Updated if contact
   * happens.
   */
  void check_ball_racket_contact(const racket &robot, const vec3 &ball_cand_pos,
                                 vec3 &ball_cand_vel);

  /**
   *
   * @brief Checks contact with ground and zeros the velocities.
   *
   * Checking contact with ground. Zeros the velocities and
   * hardsets the next candidate positions to ground level!
   *
   * @param ball_cand_pos Next candidate ball positions.
   * If contact occurs, z-position is set to floor level.
   * @param ball_cand_vel Next candidate ball velocities.
   * If contact occurs, velocities are set to zero.
   */
  void check_ball_ground_contact(vec3 &ball_cand_pos, vec3 &ball_cand_vel);

public:
  /**
   * @brief Initialize ball variables to zero.
   *
   * Initializes all ball variables (positions,velocities,spin) to ZERO.
   *
   * @param spin_flag Turn ON for spin modelling.
   * @param verbosity Turn ON for printing events (bounce, hit, etc.)
   * @param check_contacts Turn OFF for disabling contact checking during
   * integration (default ON).
   */
  TableTennis(bool spin = false, bool verbose = false,
              bool check_contacts = true);

  /**
   * @brief Initialize ball pos and velocity to input vector of size 6
   * and initialize ball spin to zero.
   * @param ball_state Initial ball state (pos and vel).
   * @param spin_flag Turn ON for spin modelling.
   * @param verbosity Turn ON for printing events (bounce, hit, etc.)
   */
  TableTennis(const vec6 &ball_state, bool spin = false, bool verbose = false);

  /** @brief Turn off contact checking so integrate ball_state will not check
   * contacts! */
  void turn_off_contact_checking();

  /**
   * @brief Checks for legal landing (on the opponents court).
   *
   * Function that returns true if ball has landed legally and bounced only
   * once on OPPONENTS court.
   * Useful for generating statistics, and on robot vs. robot mode.
   *
   * @return TRUE if ball has landed legally.
   *
   */
  bool has_legally_landed() const;

  /** @brief Checks for legal bounce of the ball (on the robot court).*/
  bool has_legally_bounced() const;

  /** @brief Checks for a hit by the racket during serve
   *
   * TODO: should check for a bounce on robot court + bounce on opponent court
   */
  bool was_legally_served() const;

  /** @brief Checks for contact with ground */
  bool touched_ground() const;

  /**
   * @brief Reset statistics of the game.
   *
   * Statistics are stored for counting return successes in SIM mode.
   */
  void reset_stats();

  /**
   * @brief Load ball prediction and other SIM parameters from a CONFIG file
   * @param file_name_relative Relative file name (base is table-tennis)
   *
   */
  void load_params(const std::string &file_name_relative);

  /** @return Ball position as a 3-vector. */
  vec3 get_ball_position() const;

  /** @return Ball velocity as a 3-vector. */
  vec3 get_ball_velocity() const;

  /** @return Return ball state as a 6-vector. */
  vec6 get_ball_state() const;

  /** @return Set ball state as a 6-vector. */
  void set_ball_state(const vec6 &ball_state);

  /** @brief Set topspin equal to argument (revolutions/sec) */
  void set_topspin(const double val);

  /**
   * @brief Reset the simulated ball state.
   *
   * Set the ball-gun somewhere behind the table and launches a ball,
   * i.e., initializes the ball state accordingly.
   * Adds noise to the initial ball launch velocity.
   * DOES NOT MODIFY SPIN (from before)!
   * Method is used for testing purposes (see Unit Tests).
   *
   * @param std Standard deviation of the initial ball pos and vel distribution.
   * @param ballgun_side Position the ballgun: 0 = LEFT, 1 = CENTER (DEFAULT), 2
   * = RIGHT.
   *
   */
  void set_ball_gun(double std, int ballgun_side = 1);

  /**
   * @brief FOURTH ORDER Symplectic integration for dt seconds.
   *
   * Unlike Symplectic Euler, this is only used for accurate and fast racket
   * dynamics computations so it integrates already the positions and velocities
   */
  void symplectic_int_fourth(const double dt);

  /**
   * @brief Integrate the ball state dt seconds later.
   * Checking contacts with environment, i.e. net, table, ground.
   * Does not check the racket!!
   * Modified: July-August 2016
   *
   * @param dt Prediction horizon.
   *
   */
  void integrate_ball_state(const double dt);

  /**
   *
   * @brief Main function to integrate the ball state (for dt seconds).
   *
   * This function is used to predict the ball state and typically is called
   * many times. Can be used to predict the landing point
   * (if there is a strike of course).
   *
   * Modified: February 2017
   *
   * Integrate the ball state dt seconds later.
   * Checking contacts with environment, i.e. racket, net, table, ground.
   *
   * @param robot_racket Racket of the robot for checking a strike
   * @param dt Prediction horizon.
   *
   * Takes around 1mu sec to run
   *
   */
  void integrate_ball_state(const racket &robot, const double dt);

  /**
   * @brief Calculate desired racket normal assuming mirror law
   *
   * Does not use any state of the table tennis class, only the parameters.
   *
   * @param v_in Incoming ball velocity
   * @param v_out Outgoing ball velocity (desired)
   * @param normal Desired normals of the racket calculated (output)
   */
  void calc_des_racket_normal(const mat &v_in, const mat &v_out,
                              mat &normal) const;

  /**
   * @brief Computes the desired outgoing velocity of the ball after possible
   * contact
   *
   * To return to a desired landing position at a desired landing time on the
   * opponents court, we calculate the desired outgoing velocities.
   * Does not use any state of the table tennis class, only the parameters.
   *
   * @param ball_land_des Desired landing position of the ball
   * @param time_land_des Time it should take for the ball to land on opponents
   * court
   * @param balls_predicted The incoming balls that are predicted (for a fixed
   * time horizon)
   * @param balls_out_vel Outgoing velocities on the predicted ball locations
   * (output)
   */
  void calc_des_ball_out_vel(const arma::vec2 &ball_land_des,
                             const double time_land_des, const bool hack,
                             const mat &balls_predicted,
                             mat &balls_out_vel) const;

  /**
   * @brief Calculate desired racket velocity given ball incoming and outgoing
   * velocities
   *
   * Assuming a mirror law.
   * Assumes no desired spin, i.e. racket velocity along the racket will be set
   * to zero
   *
   * @param vel_ball_in Incoming ball velocities already predicted
   * @param vel_ball_out Outgoing desired ball velocities already calculated
   * @param racket_normal Desired racket normals already calculated
   * @param racket_vel Desired racket velocities (output)
   */
  void calc_des_racket_vel(const mat &vel_ball_in, const mat &vel_ball_out,
                           const mat &racket_normal, mat &racket_vel) const;
};

/**
 * @brief Function that integrates a table tennis ball for an outside filter.
 *
 * Function exposes the table tennis integration to filters, e.g. an EKF.
 * They can use then to apply predict() using the this function pointer.
 *
 * Warning: spin is turned OFF!
 * Prediction with a spin model assumes that spin is kept constant
 * as changes to spin are not saved!
 *
 *
 * @param xnow Consists of current ball position and velocity.
 * @param dt Prediction horizon.
 * @param fp Function parameters, not used.
 * @return Next ball positions and velocities.
 */
vec calc_next_ball(const vec &xnow, const double dt, const void *fp);
/**
 * @brief Function that integrates a table tennis ball for an outside filter.
 * including predicting a potential racket contact.
 *
 * Overloaded function exposes the table tennis integration to filters, e.g. an
 * EKF. They can use then to apply predict() using the this function pointer.
 * This version also checks for robot's racket to predict next ball state.
 *
 * Warning: spin is turned OFF!
 *
 * @param robot Interactions with the robot will be checked in integration.
 * @param xnow Consists of current ball position and velocity.
 * @param dt Prediction horizon.
 * @return Next ball positions and velocities.
 */
vec calc_next_ball(const racket &robot, const vec &xnow, const double dt);

/**
 * @brief Function that integrates a table tennis ball for an outside filter.
 *
 * Function exposes the table tennis integration to filters, e.g. an EKF.
 * They can use then to apply predict() using the this function pointer.
 *
 * Warning: spin is turned ON!
 *
 * @param xnow Consists of current ball position and velocity.
 * @param dt Prediction horizon.
 * @param fp Function parameters are in this case the topspin value.
 * @return Next ball positions and velocities.
 */
vec calc_spin_ball(const vec &xnow, const double dt, const void *fp);

/**
 * @brief Predict ball state FORWARDS till net
 *
 * Predict the ball forwards till the net so that we can look up the
 * corresponding robot joint parameters qf, qfdot, T around net.
 *
 * Used for lookup table based trajectory generation.
 * TODO: it would NOT work backwards
 * @param ball_est
 */
void predict_till_net(vec6 &ball_est);

} // namespace player
