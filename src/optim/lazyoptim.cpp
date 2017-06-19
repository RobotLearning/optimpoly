/**
 * @file lazyoptim.cpp
 *
 * @brief NLOPT polynomial optimization functions for LAZY PLAYER are stored here.
 *
 *
 *  Created on: Sep 11, 2016
 *      Author: okan
 */

#include <armadillo>
#include <iostream>
#include "constants.h"
#include "table.h"
#include "utils.h"
#include "stdlib.h"
#include "math.h"
#include "kinematics.h"
#include "optim.h"

#define INEQ_HIT_CONSTR_DIM 3
#define INEQ_LAND_CONSTR_DIM 8 //11
#define INEQ_JOINT_CONSTR_DIM 2*NDOF + 2*NDOF

static double costfunc(unsigned n, const double *x, double *grad, void *my_func_params);
static double punish_land_robot(const double *xland,
								const double *xnet,
								const double Rland,
								const double Rnet);
static void land_ineq_constr(unsigned m, double *result, unsigned n,
		                     const double *x, double *grad,
		                     void *my_func_params);
static void hit_ineq_constr(unsigned m, double *result, unsigned n, const double *x, double *grad,
		                     void *my_func_params);
static void racket_contact_model(const double* racketVel,
		                         const double* racketNormal, double* ballVel);
static void interp_ball(const optim_des *params, const double T,
		                double *ballpos, double *ballvel);


LazyOptim::LazyOptim(double qrest_[NDOF], double lb_[], double ub_[], bool land_, bool grad_)
                          : FocusedOptim() { //FocusedOptim(qrest_, lb_, ub_) {

	GRAD_BASED_OPT = grad_;
	param_des = nullptr;
	const_vec(NDOF,1.0,w.R_strike);
	w.R_net = 1e2;
	w.R_hit = 1e2;
	w.R_land = 1e2;

	for (int i = 0; i < NDOF; i++) {
		qrest[i] = qrest_[i];
		ub[i] = ub_[i];
		lb[i] = lb_[i];
	}

	if (land_) {
		land = true;
		set_land_constr();
	}
	else {
		land = false;
		set_hit_constr();
	}
}

/**
 * @brief Setting LANDING constraints, i.e., not just hitting the ball.
 *
 * Here we consider the full landing constraints where the ball has
 * to pass over the net and land on the opponent's court.
 * We can study the simpler optimization problem of hitting the ball
 * by switching to pure hitting constraints.
 */
void LazyOptim::set_land_constr() {

	double tol_ineq_land[INEQ_LAND_CONSTR_DIM];
	double tol_ineq_joint[INEQ_JOINT_CONSTR_DIM];
	const_vec(INEQ_LAND_CONSTR_DIM,1e-2,tol_ineq_land);
	const_vec(INEQ_JOINT_CONSTR_DIM,1e-3,tol_ineq_joint);

	// LN = does not require gradients //
	opt = nlopt_create(NLOPT_LN_COBYLA, 2*NDOF+1);
	nlopt_set_min_objective(opt, costfunc, this);
	nlopt_add_inequality_mconstraint(opt, INEQ_LAND_CONSTR_DIM,
			land_ineq_constr, this, tol_ineq_land);
	nlopt_add_inequality_mconstraint(opt, INEQ_JOINT_CONSTR_DIM,
			joint_limits_ineq_constr, this, tol_ineq_joint);
	nlopt_set_xtol_rel(opt, 1e-2);
}

/**
 * @brief Switch to only HITTING constraints, not landing.
 *
 * We can study the simpler optimization problem of hitting the ball
 * by switching to pure hitting constraints.
 */
void LazyOptim::set_hit_constr() {

	double tol_ineq_hit[INEQ_HIT_CONSTR_DIM];
	double tol_ineq_joint[INEQ_JOINT_CONSTR_DIM];
	const_vec(INEQ_HIT_CONSTR_DIM,1e-2,tol_ineq_hit);
	const_vec(INEQ_JOINT_CONSTR_DIM,1e-3,tol_ineq_joint);

	if (GRAD_BASED_OPT) {
		printf("Using Gradient Based Optimizer!\n");
		opt = nlopt_create(NLOPT_LD_SLSQP, OPTIM_DIM);
		/*nlopt_opt local_opt = nlopt_create(NLOPT_LD_SLSQP, OPTIM_DIM);
		nlopt_set_xtol_rel(local_opt, 1e-2);
		nlopt_set_local_optimizer(opt,local_opt);*/
		param_des = new optim_des();
		generate_tape();
	}
	else {	// LN = does not require gradients //
		param_des = nullptr;
		opt = nlopt_create(NLOPT_LN_COBYLA, OPTIM_DIM);
	}

	nlopt_set_min_objective(opt, costfunc, this);
	nlopt_add_inequality_mconstraint(opt, INEQ_HIT_CONSTR_DIM,
			hit_ineq_constr, this, tol_ineq_hit);
	//nlopt_add_inequality_mconstraint(opt, INEQ_JOINT_CONSTR_DIM,
	//		joint_limits_ineq_constr, this, tol_ineq_joint);
	nlopt_set_xtol_rel(opt, 1e-2);
}


void LazyOptim::finalize_soln(const double x[], double time_elapsed) {

	if (x[2*NDOF] > fmax(time_elapsed/1e3,0.05)) {
		// initialize first dof entries to q0
		for (int i = 0; i < NDOF; i++) {
			qf[i] = x[i];
			qfdot[i] = x[i+NDOF];
		}
		T = x[2*NDOF];
		if (detach)
			T -= (time_elapsed/1e3);
		update = true;
	}
	//trigger_optim();
}


/**
 * @brief Calculate the net hitting time and the landing time
 * Assuming that there was an interaction (hitting) between racket and ball.
 *
 * Since we call this function multiple times within one optimization step
 * we store the optimization variable x,
 * landTime and netTime variables and return these when the optimization
 * variable is the same (meaning optimization algorithm did not update it yet).
 *
 *
 */
void LazyOptim::calc_times(const double x[]) { // ball projected to racket plane

	static double qf[NDOF];
	static double qfdot[NDOF];
	static double vel[NCART];
	static double normal[NCART];
	static double pos[NCART];
	static double ballpos[NCART];
	static double ballvel[NCART];
	static double g = -9.8;
	static double net_y = dist_to_table - table_length/2.0;
	static double table_z = floor_level - table_height;
	double discr = 0.0;
	double d;

	if (!vec_is_equal(OPTIM_DIM,x,x_last)) {
		// extract state information from optimization variables
		for (int i = 0; i < NDOF; i++) {
			qf[i] = x[i];
			qfdot[i] = x[i + NDOF];
		}
		calc_racket_state(qf, qfdot, pos, vel, normal);
		interp_ball(this->param_des, x[2*NDOF], ballpos, ballvel);
		// calculate deviation of ball to racket - hitting constraints
		calc_hit_distance(ballpos,pos,normal);
		racket_contact_model(vel, normal, ballvel);
		t_net = (net_y - ballpos[Y])/ballvel[Y];
		x_net[Z] = ballpos[Z] + t_net * ballvel[Z] + 0.5*g*t_net*t_net;
		x_net[X] = ballpos[X] + t_net * ballvel[X];

		// calculate ball planned landing
		d = ballpos[Z] - table_z;
		if (sqr(ballvel[Z]) > 2*g*d) {
			discr = sqrt(sqr(ballvel[Z]) - 2*g*d);
		}
		t_land = fmax((-ballvel[Z] - discr)/g,(-ballvel[Z] + discr)/g);
		t_net = (net_y - ballpos[Y])/ballvel[Y];
		x_land[X] = ballpos[X] + t_land * ballvel[X];
		x_land[Y] = ballpos[Y] + t_land * ballvel[Y];

		make_equal(OPTIM_DIM,x,x_last);
	}
}


/**
 * @brief Calculate deviation and projection norms of ball to racket.
 *
 * Calculates normal and projected distance from ball to racket
 * For satisfying hitting constraints.
 *
 */
void LazyOptim::calc_hit_distance(const double ball_pos[],
		                          const double racket_pos[],
								  const double racket_normal[]) {

	double e[NCART];
	for (int i = 0; i < NCART; i++) {
		e[i] = ball_pos[i] - racket_pos[i];
	}
	dist_b2r_norm = inner_prod(NCART,racket_normal,e);
	dist_b2r_proj = sqrt(inner_prod(NCART,e,e) - dist_b2r_norm*dist_b2r_norm);

}

/*
 * Initialize ADOLC operations
 * Objective has easily computed gradients so no need to auto-diff
 */
void LazyOptim::generate_tape() {

	int n = OPTIM_DIM;
	int m = INEQ_HIT_CONSTR_DIM;

	adouble *x_auto = new adouble[n];
	adouble f_auto = 1;
	adouble *g_auto = new adouble[m];
	double *x = new double[n];
	double dummy;
	jac = new double*[m];
	for (int i = 0; i < m; i++)
		jac[i] = new double[n];

	init_rest_soln(x);

	trace_on(1);
	for (int i = 0; i < n; i++)
		x_auto[i] <<= x[i];
	f_auto = costfunc(n,x_auto,nullptr,this);
	f_auto >>= dummy;
	trace_off();

	trace_on(2);
	for (int i = 0; i < n; i++)
		x_auto[i] <<= x[i];

	kinematics_eq_constr(m,g_auto,n,x_auto,nullptr,this);
	for(int i = 0; i < m; i++)
		g_auto[i] >>= dummy;
	trace_off();

	delete[] x_auto;
	delete[] x;
	delete[] g_auto;
}

double LazyOptim::test_soln(const double x[]) const {

	double max_viol;
	// give info on constraint violation
	static double table_xmax = table_width/2.0;
	static double table_ymax = dist_to_table - table_length;
	static double wall_z = 1.0;
	static double net_y = dist_to_table - table_length/2.0;
	static double net_z = floor_level - table_height + net_height;
	static int count = 0;
	double *grad = 0;
	static double land_violation[INEQ_LAND_CONSTR_DIM];
	static double lim_violation[INEQ_JOINT_CONSTR_DIM]; // joint limit violations on strike and return
	joint_limits_ineq_constr(INEQ_JOINT_CONSTR_DIM, lim_violation, OPTIM_DIM, x, grad, (void*)this);
	land_ineq_constr(INEQ_LAND_CONSTR_DIM, land_violation, OPTIM_DIM, x, grad, (void*)this);
	double cost = costfunc(OPTIM_DIM, x, grad, (void*)this);

	if (verbose) {
		// give info on solution vector
		printf("Optim count: %d\n", (++count));
		print_optim_vec(x);
		printf("f = %.2f\n",cost);
		printf("Hitting constraints (b2r):\n");
		printf("Distance along normal: %.2f\n",dist_b2r_norm);
		printf("Distance along racket: %.2f\n",dist_b2r_proj);
		printf("Landing constraints:\n");
		printf("NetTime: %f\n", t_net);
	    printf("LandTime: %f\n", t_land);
		printf("Below wall by : %f\n", wall_z - x_net[Z]);
		printf("Above net by: %f\n", x_net[Z] - net_z);
		printf("X: between table limits by [%f, %f]\n", table_xmax - x_land[X], x_land[X] + table_xmax);
		printf("Y: between table limits by [%f, %f]\n", net_y - x_land[Y], x_land[Y] - table_ymax);
		for (int i = 0; i < INEQ_JOINT_CONSTR_DIM; i++) {
			if (lim_violation[i] > 0.0) {
				printf("Joint limit violated by %.2f on joint %d\n", lim_violation[i], i%NDOF + 1);
			}
		}
	}

	max_viol = max_array(lim_violation,INEQ_JOINT_CONSTR_DIM);
	if (land)
		max_viol = fmax(max_viol,max_array(land_violation,INEQ_LAND_CONSTR_DIM));
	return max_viol;
}

/*
 * Calculates the cost function for table tennis Lazy Player (LP)
 * to find spline (3rd order strike+return) polynomials
 */
static double costfunc(unsigned n, const double *x, double *grad, void *my_func_params) {

	static double J1, Jhit, Jland;
	static double a1[NDOF], a2[NDOF];
	double T = x[2*NDOF];

	LazyOptim* opt = (LazyOptim*) my_func_params;
	double *q0 = opt->q0;
	double *q0dot = opt->q0dot;
	weights w = opt->w;

	// calculate the polynomial coeffs which are used in the cost calculation
	calc_strike_poly_coeff(q0,q0dot,x,a1,a2);

	// calculate the landing time
	opt->calc_times(x);

	J1 = T * (3*T*T*inner_winv_prod(NDOF,w.R_strike,a1,a1) +
			3*T*inner_winv_prod(NDOF,w.R_strike,a1,a2) +
			inner_winv_prod(NDOF,w.R_strike,a2,a2));

	Jhit = w.R_hit * opt->dist_b2r_proj;

	if (opt->land)
		Jland = punish_land_robot(opt->x_land,opt->x_net,w.R_land, w.R_net);

	//std::cout << J1 << "\t" << Jhit << "\t" << Jland << std::endl;

	return J1 + Jhit + Jland;
}

/*
 * Calculates the cost function for table tennis Lazy Player (LP)
 * to find spline (3rd order strike+return) polynomials
 */
static adouble costfunc(unsigned n, const adouble *x, double *grad, void *my_func_params) {

	static adouble J1, Jhit, Jland;
	static adouble a1[NDOF], a2[NDOF];
	adouble T = x[2*NDOF];

	LazyOptim* opt = (LazyOptim*) my_func_params;
	double *q0 = opt->q0;
	double *q0dot = opt->q0dot;
	weights w = opt->w;

	// calculate the polynomial coeffs which are used in the cost calculation
	calc_strike_poly_coeff(q0,q0dot,x,a1,a2);

	// calculate the landing time
	opt->calc_times(x);

	J1 = T * (3*T*T*inner_winv_prod(NDOF,w.R_strike,a1,a1) +
			3*T*inner_winv_prod(NDOF,w.R_strike,a1,a2) +
			inner_winv_prod(NDOF,w.R_strike,a2,a2));

	Jhit = w.R_hit * opt->dist_b2r_proj;

	if (opt->land)
		Jland = punish_land_robot(opt->x_land,opt->x_net,w.R_land, w.R_net);

	//std::cout << J1 << "\t" << Jhit << "\t" << Jland << std::endl;

	return J1 + Jhit + Jland;
}

/*
 * Punish the robot sufficiently as to induce proper landing behaviour
 *
 * The desired landing location chosen is the center of the table
 */
static double punish_land_robot(const double *xland,
								const double *xnet,
								const double Rland,
								const double Rnet) {
	// desired landing locations
	static double x_des_net = 0.0;
	static double z_des_net = floor_level - table_height + net_height + 1.0;

	return sqr(xnet[Z] - z_des_net)*Rnet + sqr(xnet[X] - x_des_net)*Rnet;

}

/*
 * This is the constraint that makes sure we land the ball
 *
 */
static void land_ineq_constr(unsigned m, double *result, unsigned n, const double *x, double *grad,
		                     void *my_func_params) {

	static double table_xmax = table_width/2.0;
	static double table_ymax = dist_to_table - table_length;
	static double wall_z = 1.0;
	static double net_y = dist_to_table - table_length/2.0;
	static double net_z = floor_level - table_height + net_height;

	LazyOptim* opt = (LazyOptim*)my_func_params;
	opt->calc_times(x);

	result[0] = -opt->dist_b2r_norm;
	result[1] = opt->dist_b2r_norm - ball_radius;
	result[2] = opt->dist_b2r_proj - racket_radius;
	result[3] = opt->x_net[Z] - wall_z;
	result[4] = -opt->x_net[Z] + net_z;
	result[5] = opt->x_net[X] - table_xmax;
	result[6] = -opt->x_net[X] - table_xmax;
	result[7] = -opt->t_net;
}

/*
 * This is the constraint that makes sure we only TOUCH the ball
 *
 */
static void hit_ineq_constr(unsigned m, double *result, unsigned n, const double *x, double *grad,
		                     void *my_func_params) {

	static double qf[NDOF];
	static double qfdot[NDOF];
	static double vel[NCART];
	static double normal[NCART];
	static double pos[NCART];
	static double ballpos[NCART];
	static double ballvel[NCART];

	LazyOptim* opt = (LazyOptim*)my_func_params;

	for (int i = 0; i < NDOF; i++) {
		qf[i] = x[i];
		qfdot[i] = x[i + NDOF];
	}
	interp_ball(opt->param_des, x[2*NDOF], ballpos, ballvel);
	calc_racket_state(qf, qfdot, pos, vel, normal);
	// calculate deviation of ball to racket - hitting constraints
	opt->calc_hit_distance(ballpos,pos,normal);

	result[0] = -opt->dist_b2r_norm;
	result[1] = opt->dist_b2r_norm - ball_radius;
	result[2] = opt->dist_b2r_proj - racket_radius;
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
 * The outgoing ball velocity is post-multiplied by some constants \mu < 1
 * to account for ball drag
 *
 *
 */
static void racket_contact_model(const double* racketVel,
		const double* racketNormal, double* ballVel) {

	static const double mult_vel[NCART] = {0.9,0.8,0.83};
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
	vec_plus(NCART,normalMultSpeed,ballVel);
	for (int i = 0; i < NCART; i++ ) {
		ballVel[i] *= mult_vel[i];
	}
}

/*
 * First order hold to interpolate linearly at time T
 * between racket pos,vel,normal entries
 *
 * IF T is nan, racket variables are assigned to zero-element of
 * relevant racket entries
 *
 */
static void interp_ball(const optim_des* data, const double T, double *ballpos, double *ballvel) {

    const double dt = data->dt;
	if (std::isnan(T)) {
		printf("Warning: T value is nan!\n");
		for(int i = 0; i < NCART; i++) {
			ballpos[i] = data->ball_pos(i,0);
			ballvel[i] = data->ball_vel(i,0);
		}
	}
	else {
		const int Nmax = data->Nmax;
		unsigned N = (int) (T/dt);
		double Tdiff = T - N*dt;
		for (int i = 0; i < NCART; i++) {
			if (N < Nmax - 1) {
				ballpos[i] = data->ball_pos(i,N) +
						(Tdiff/dt) * (data->ball_pos(i,N+1) - data->ball_pos(i,N));
				ballvel[i] = data->ball_vel(i,N) +
						(Tdiff/dt) * (data->ball_vel(i,N+1) - data->ball_vel(i,N));
			}
			else {
				ballpos[i] = data->ball_pos(i,Nmax-1);
				ballvel[i] = data->ball_vel(i,Nmax-1);
			}
		}
	}
}
