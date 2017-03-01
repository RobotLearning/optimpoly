/*
 * Table tennis common includes all the
 * functions used for modelling table tennis interactions between
 * ball, table, racket (and robots, including possible opponents).
 *
 * Most of them are used in simulation mode, such as
 * the simulate_ball function.
 *
 * We support in simulation opponent modelling, simple spin modelling, etc.
 *
 *
 */

#include <armadillo>
#include "tabletennis.h"

/*
 *
 * Initialize ball pos and velocity to input vector of size 6
 * and initialize ball spin to zero
 *
 * TODO: initialized racket variables to inf!
 *
 */
TableTennis::TableTennis(const vec6 & ball_state, bool spin_flag, bool verbosity) {

	SPIN_MODE = spin_flag;
	VERBOSE = verbosity;
	LAND = false;
	HIT = false;
	ball_pos = ball_state(span(X,Z));
	ball_vel = ball_state(span(DX,DZ));
	ball_spin = zeros<vec>(3);

	//init_topspin();
}

/*
 *
 * Initialize ball variables to zero and init racket variables to inf!
 *
 */
TableTennis::TableTennis(bool spin_flag, bool verbosity) {

	SPIN_MODE = spin_flag;
	VERBOSE = verbosity;
	LAND = false;
	HIT = false;
	ball_pos = zeros<vec>(3);
	ball_vel = zeros<vec>(3);
	ball_spin = zeros<vec>(3);

	//init_topspin();
}

/*
 * Initialize constant angular velocity (a.k.a. spin) for the spinning ball
 */
void TableTennis::init_topspin() {

	if (SPIN_MODE) {
		ball_spin(X) = -50*2*datum::pi; // constant 3000 rpm topsin
		// others are zero
	}
	else {
		// do nothing, they are all zero - spinless model
	}
}

/*
 *
 * Reset the simulated ball state.
 * Set the ball-gun somewhere behind the table.
 * Adds noise to the initial ball launch velocity.
 * Unless you set the seed to a reasonable value - Okan. *
 *
 */
void TableTennis::set_ball_state(double std) {

	vec3 ballgun;
	ballgun << table_center + 0.4 << endr
			<< dist_to_table - table_length - 0.2 << endr
			<< floor_level - table_height + 0.15 << endr;

	vec3 rand_ball_pos = ballgun + std * randn<vec>(3);
	vec3 good_ball_vel;
	good_ball_vel << -1.08 << endr << 4.80 << endr << 3.84 << endr;
	vec3 rand_ball_vel = good_ball_vel + std * randn<vec>(3);

	this->ball_pos = rand_ball_pos;
	this->ball_vel = rand_ball_vel;
	this->ball_spin = zeros<vec>(3);
	this->LAND = false;
	this->HIT = false;
}

/*
 * Return the ball position as a 3-vector
 */
vec3 TableTennis::get_ball_position() const {

	return this->ball_pos;
}

/*
 * Return the ball velocity as a 3-vector
 */
vec3 TableTennis::get_ball_velocity() const {

	return this->ball_vel;
}

/*
 * Modified: February 2017
 *
 * Integrate the ball state dt seconds later.
 * Checking contacts with environment, i.e. racket, net, table, ground.
 * TODO: implementing Symplectic Euler, implement RK4 also!
 *
 * Takes around 1mu sec to run
 *
 */
void TableTennis::integrate_ball_state(const racket & robot_racket,
		                               const double dt) {

	// Symplectic Euler for No-Contact-Situation (Flight model)
	vec3 ball_acc = flight_model();
	vec3 ball_cand_pos, ball_cand_vel;
	symplectic_euler(ball_acc, ball_cand_pos, ball_cand_vel, dt);

	if (CHECK_CONTACTS) {
		check_contact(robot_racket,ball_cand_pos,ball_cand_vel);
	}

	// Pass the computed ball variables to ball pos and vel
	ball_pos = ball_cand_pos;
	ball_vel = ball_cand_vel;
}

/*
 * Modified: July-August 2016
 *
 * Integrate the ball state dt seconds later.
 * Checking contacts with environment, i.e. net, table, ground.
 * Does not check the racket!!
 *
 */
void TableTennis::integrate_ball_state(double dt) {

	// Symplectic Euler for No-Contact-Situation (Flight model)
	vec3 ball_acc = flight_model();
	vec3 ball_cand_pos, ball_cand_vel;
	symplectic_euler(ball_acc, ball_cand_pos, ball_cand_vel, dt);

	if (CHECK_CONTACTS) {
		check_ball_table_contact(ball_cand_pos,ball_cand_vel);
		check_ball_net_contact(ball_cand_pos,ball_cand_vel);
		check_ball_ground_contact(ball_cand_pos,ball_cand_vel);
	}

	// Pass the computed ball variables to ball pos and vel
	ball_pos = ball_cand_pos;
	ball_vel = ball_cand_vel;
}

/*
 * Spinning nonlinear ball flight model
 * [including airdrag and Magnus force]
 *
 * Returns ball accelerations as a 3-vector (ARMA class)
 *
 */
vec3 TableTennis::flight_model() const {

	vec3 ball_acc = drag_flight_model();
	// add Magnus force
	vec3 magnus = cross(ball_spin,ball_vel); // acceleration due to magnus force

	magnus *= Clift;
	ball_acc += magnus;
	return ball_acc;

}

/*
 * Spin free nonlinear flight model including airdrag (Cdrag as parameter)
 */
vec3 TableTennis::drag_flight_model() const {

	double velBall = norm(ball_vel);
	vec3 ball_acc;
	ball_acc(X) = -ball_vel(X) * Cdrag * velBall;
	ball_acc(Y) = -ball_vel(Y) * Cdrag * velBall;
	ball_acc(Z) = gravity - ball_vel(Z) * Cdrag * velBall;

	return ball_acc;
}

/*
 * First integrating the accelerations to velocities by dt
 * Then integrating the velocities to positions by dt
 * These (pos and vel) are kept in the ball candidate vector
 *
 * We dont directly update the ball states, this is done in
 * integrate_ball_state. This way we can check for contacts
 * in between.
 *
 */
void TableTennis::symplectic_euler(const vec3 & ball_acc,
		vec3 & ball_next_pos, vec3 & ball_next_vel, const double dt) const {

	// ball candidate velocities
	ball_next_vel(X) = ball_vel(X) + ball_acc(X) * dt;
	ball_next_vel(Y) = ball_vel(Y) + ball_acc(Y) * dt;
	ball_next_vel(Z) = ball_vel(Z) + ball_acc(Z) * dt;

	// ball candidate positions
	ball_next_pos(X) = ball_pos(X) + ball_next_vel(X) * dt;
	ball_next_pos(Y) = ball_pos(Y) + ball_next_vel(Y) * dt;
	ball_next_pos(Z) = ball_pos(Z) + ball_next_vel(Z) * dt;
}

/*
 * Check if a contact will occur,
 *
 * Checking against table, net, racket, ground, and possibly
 * a simulated human opponent!
 *
 */
void TableTennis::check_contact(const racket & robot_racket,
		                        vec3 & ball_cand_pos,
		                        vec3 & ball_cand_vel) {

	// Check contact to table
	check_ball_table_contact(ball_cand_pos,ball_cand_vel);
	// Check contact with net
	check_ball_net_contact(ball_cand_pos,ball_cand_vel);
	// Check contact with racket
	check_ball_racket_contact(robot_racket,ball_cand_pos,ball_cand_vel);
	// Check if it hits the ground...
	check_ball_ground_contact(ball_cand_pos,ball_cand_vel);
}

/*
 * Condition to determine if ball hits the table
 * Useful for prediction including a rebound model
 * Useful also in KF/EKF filtering.
 *
 *
 */
void TableTennis::check_ball_table_contact(const vec3 & ball_cand_pos, vec3 & ball_cand_vel) {

	static const double contact_table_level = floor_level - table_height + ball_radius;
	static const double table_human_end = dist_to_table - table_length;

	// check to see if ball is over the table
	if ((ball_cand_pos(Y) > table_human_end) && (ball_cand_pos(Y) < dist_to_table) &&
			(fabs(ball_cand_pos(X) - table_center) <= table_width/2.0)) {
		// check if the ball hits the table coming from above
		if ((ball_cand_pos(Z) <= contact_table_level) && (ball_cand_vel(Z) < 0.0)) {
			LAND = check_landing(ball_cand_pos(Y),HIT,VERBOSE);
			table_contact_model(ball_spin,ball_cand_vel,SPIN_MODE);
		}
	}
}

/*
 * Checks contact with net.
 * Curious way to check it: if the net's distance to integrated y-state and distance to current y-state
 * signs do not match, it means that the ball is in contact with the net
 */
void TableTennis::check_ball_net_contact(vec3 & ball_cand_pos, vec3 & ball_cand_vel) const {


	static const double contact_table_level = floor_level - table_height + ball_radius;
	static const double net_dist_robot = dist_to_table - 0.5 * table_length;
	static double dist_state_net;
	static double dist_cand_net;

	// Check contact with net
	if ((ball_cand_pos(Z) <= contact_table_level + net_height)
			&& (fabs(ball_cand_pos(X)) <= table_width/2.0 + net_overhang)) {
		dist_state_net = net_dist_robot - ball_pos(Y);
		dist_cand_net = net_dist_robot - ball_cand_pos(Y);

		// If on the other side of the net after integration
		// apply super simplistic model for contact with net
		if (dist_state_net >= 0.0 && dist_cand_net < 0.0) {
			if (VERBOSE)
				std::cout << "Touches the net!" << std::endl;
			// Reflect to Front
			ball_cand_vel(Y) *= -net_restitution;
			ball_cand_pos(Y) = net_dist_robot + (0.5 * net_thickness + ball_radius);
		}
		else if (dist_state_net < 0.0 && dist_cand_net >= 0.0){
			if (VERBOSE)
				std::cout << "Touches the net!" << std::endl;
			// Reflect to Back
			ball_cand_vel(Y) *= -net_restitution;
			ball_cand_pos(Y) = net_dist_robot - (0.5 * net_thickness + ball_radius);
		}
	}
}

/*
 * Checks contact with racket.
 * If there is contact, then the predicted state will be transformed according to a racket-ball contact model.
 * If racket is going to hit the ball, static hit variable is set to TRUE so that
 * we do not hit it the next time.
 *
 *
 */
void TableTennis::check_ball_racket_contact(const racket & robot_racket,
		                                    const vec3 & ball_cand_pos,
		                                    vec3 & ball_cand_vel) {

	vec3 racket2ball = robot_racket.pos - ball_cand_pos;
	double normal_dist_ball2racket = dot(robot_racket.normal,racket2ball);
	double parallel_dist_ball2racket = sqrt(fabs(dot(racket2ball,racket2ball)
							- normal_dist_ball2racket*normal_dist_ball2racket));

	// check for contact with racket
	if ((parallel_dist_ball2racket < racket_radius) &&
			(fabs(normal_dist_ball2racket) < ball_radius) && !HIT) {
		HIT = true;
		if (VERBOSE)
			std::cout << "Contact with racket!" << std::endl;
		// Reflect to back
		racket_contact_model(robot_racket.vel, robot_racket.normal, ball_cand_vel);
	}
}

/*
 * Checks contact with ground and zeros the velocities
 *
 * Also hardsets the candidate positions to ground level!
 */
void TableTennis::check_ball_ground_contact(vec3 & ball_cand_pos, vec3 & ball_cand_vel) const {

	static double ball_last_z_pos;
	if (ball_cand_pos(Z) <= floor_level) {
		if (VERBOSE && (ball_pos(Z) != ball_last_z_pos)) // we dont want to print all the time
			std::cout << "Contact with ground Zeroing the velocities!" << std::endl;
		// zero the velocities
		ball_cand_vel = zeros<vec>(3);
		ball_cand_pos(Z) = floor_level;
	}
	ball_last_z_pos = ball_pos(Z);

}

/*
 * Returns TRUE if ball has landed legally
 *
 */
bool TableTennis::has_landed() const {

	return this->LAND;
}

/*
 * Checks for legal landing on the opponents court
 * if there was already a hit then the bounce location is checked.
 *
 */
bool check_landing(const double ball_y, const bool hit, const bool verbose) {

	bool land = false;
	if (verbose) {
		if (ball_y < dist_to_table - table_length/2.0)
			std::cout << "Bounces on opponents court!" << std::endl;
		else
			std::cout << "Bounces on robot court!" << std::endl;
	}
	if (hit) {
		if (ball_y < dist_to_table - table_length/2.0) {
			land = true;
			if (verbose) {
				std::cout << "Legal land!" << std::endl;
			}
		}
		else {
			land = false;
			if (verbose) {
				std::cout << "Illegal ball! Lost a point..." << std::endl;
			}
		}
	}
	return land;
}

/*
 * Forms the rotation matrix that corresponds to the quaternion
 *
 *
 */
mat33 quat2mat(const vec4 & q) {

	mat33 R;
	R(X,X) = 2*q(X)*q(X) - 1 + 2*q(Y)*q(Y);
	R(X,Y) = 2*q(Y)*q(Z) - 2*q(X)*q(W);
	R(X,Z) = 2*q(Y)*q(W) + 2*q(X)*q(Z);
	R(Y,X) = 2*q(Y)*q(Z) + 2*q(X)*q(W);
	R(Y,Y) = 2*q(X)*q(X) - 1 + 2*q(Z)*q(Z);
	R(Y,Z) = 2*q(Z)*q(W) - 2*q(X)*q(Y);
	R(Z,X) = 2*q(Y)*q(W) - 2*q(X)*q(Z);
	R(Z,Y) = 2*q(Z)*q(W) + 2*q(X)*q(Y);
	R(Z,Z) = 2*q(X)*q(X) - 1 + 2*q(W)*q(W);
	return R;
}

/*
 * Update the incoming ball velocity with outgoing ball velocity using MIRROR LAW
 *
 * The racket contact model in vector form is O = I + (1 + eps_R)*N*N'*(V - I)
 * where I is the incoming ball velocity
 *       N is the racket normal
 *       V is the racket velocity
 *       eps_R is the coefficient of restitution of the racket
 *
 * TODO: add a spinning contact model
 *
 */
void racket_contact_model(const vec3 & racket_vel, const vec3 & racket_normal, vec3 & ball_vel) {

	double speed = (1 + CRR) * dot(racket_normal, racket_vel - ball_vel);
	ball_vel += speed * racket_normal;
}

/*
 * Simple contact model that uses spin if spin mode is turned on.
 * Assuming roll contact instead of slide contact in rebound calculations for simplicity.
 * Spin vector is also modified.
 * Coeff of restitution and friction used.
 */
void table_contact_model(vec3 & ball_spin, vec3 & ball_vel,
		                 bool spin_flag) {

	if (spin_flag) { // if spin mode is on ballvec is not a null pointer
		// compute effect on spin
		ball_spin(X) -= (3*(1-CFTX)/(2*ball_radius))*ball_vel(Y) + (3*(1-CFTX)/2)*ball_spin(X);
		ball_spin(Y) += (3*(1-CFTY)/(2*ball_radius))*ball_vel(X) - (3*(1-CFTY)/2)*ball_spin(Y);
		// in z-direction spin is preserved
		// compute effect on velocity
		ball_vel(Z) = -CRT * ball_vel(Z);
		ball_vel(Y) = CFTY * ball_vel(Y) - (1-CFTY) * ball_radius * ball_spin(X);
		ball_vel(X) = CFTX * ball_vel(X) + (1-CFTX) * ball_radius * ball_spin(Y);
	}
	else {
		// reflect ball velocity
		ball_vel(Z) = -CRT * ball_vel(Z);
		ball_vel(Y) = CFTY * ball_vel(Y);
		ball_vel(X) = CFTX * ball_vel(X);
	}
}

/*
 * Friend function that exposes table tennis ball integration function
 * to an outside filter (e.g. EKF) except for racket contact
 *
 * Warning: spin is turned off!
 *
 * xnow consists of current ball position and velocity
 *
 *
 */
vec calc_next_ball(const vec & xnow, double dt) {

	TableTennis tennis = TableTennis(xnow);
	tennis.integrate_ball_state(dt);
	static vec6 out = zeros<vec>(6);
	out(span(X,Z)) = tennis.ball_pos;
	out(span(DX,DZ)) = tennis.ball_vel;
	return out;
}

/*
 * Friend function that exposes table tennis ball integration function
 * to an outside filter (e.g. EKF) including potential racket contact
 *
 * Warning: spin is turned off!
 *
 * xnow consists of current ball position and velocity
 *
 *
 */
vec calc_next_ball(const racket & robot, const vec & xnow, double dt) {

	TableTennis tennis = TableTennis(xnow);
	tennis.integrate_ball_state(robot,dt);
	static vec6 out = zeros<vec>(6);
	out(span(X,Z)) = tennis.ball_pos;
	out(span(DX,DZ)) = tennis.ball_vel;
	return out;
}
