/**
 * @file lazyoptim.cpp
 *
 * @brief NLOPT polynomial optimization functions for LAZY PLAYER are stored here.
 *
 * Lazy Player makes use of FIXED PLAYER to seed good trajectories.
 *
 *  Created on: Sep 11, 2016
 *      Author: okan
 */

#include "constants.h"
#include "table.h"
#include "utils.h"
#include "stdlib.h"
#include "math.h"
#include "kinematics.h"
#include "optim.h"

#define INEQ_LAND_CONSTR_DIM 11
#define INEQ_JOINT_CONSTR_DIM 2*NDOF + 2*NDOF

static bool firsttime[3]; //! TODO: remove this global var!

static double nlopt_optim_lazy(lazy_data *data, optim *params);
static void print_input_structs(lazy_data* data, optim* params);
static double costfunc(unsigned n, const double *x, double *grad, void *my_func_params);
static double punish_land_robot(const double *xland,
								const double xnet,
								const double *Rland,
								const double Rnet);
static double punish_wait_robot(const lazy_data* data,
								double *Rwait,
								double *qrest);
static void land_ineq_constr(unsigned m, double *result, unsigned n, const double *x, double *grad,
		                     void *my_func_params);
static void joint_limits_ineq_constr(unsigned m, double *result,
		unsigned n, const double *x, double *grad, void *my_func_params);
static void calc_strike_poly_coeff(const double *q0, const double *q0dot, const double *x,
		                    double *a1, double *a2);
static void calc_return_poly_coeff(const double *q0,
		                           const double *x, const double T,
		                           const double *Rwait,
								   double *a1, double *a2);
static void calc_strike_extrema_cand(const double *a1, const double *a2, const double T,
		                      const double *q0, const double *q0dot,
		                      double *joint_max_cand, double *joint_min_cand);
static void calc_return_extrema_cand(double *a1, double *a2, const double *x, double landTime,
		                     double *joint_max_cand, double *joint_min_cand);
static double test_optim(double *x, lazy_data* data, bool verbose);
static void modify_ball_outgoing_vel(double* ballVel);
static void calc_times(const lazy_data* data,
		               const double *x,
					   double *netTime,
					   double *landTime,
					   double *xnet,
					   double *xland,
					   double *ball_norm, // ball to racket normal distance
					   double *ball_proj);
static void calc_rest_robot(const lazy_data* data,
		                        const double *x,
		                        const double *Rret,
								const double *a3,
								const double *a2,
								const double landTime,
								double *qrest);
static void finalize_soln(const double* x,
		lazy_data* data, optim * params, double time_elapsed);
static void set_penalty_matrices(weights * pen);
static void racket_contact_model(double* racketVel, double* racketNormal, double* ballVel);
static void interp_ball(double** ballpred, const double T,
		                const double dt, const int Nmax,
						double *ballpos, double *ballvel);
static void init_soln(const optim * params, double x[OPTIM_DIM]);
static void init_right_posture(double* qwait);


/**
 * @brief Launch LAZY (also known as DEFENSIVE) PLAYER trajectory generation.
 *
 * Multi-threading entry point for the NLOPT optimization.
 * The optimization problem is solved online using COBYLA (see NLOPT).
 * As a trick, we first launch the FOCUSED PLAYER and then adapt the
 * trajectories using LAZY PLAYER optimization (only inequalities here).
 *
 * @param ballpred Ball prediction matrix (as double pointer).
 * @param coparams Co-optimization parameters held fixed during optimization.
 * @param racketdata Predicted racket position,velocity and normals.
 * @param params Optimization parameters updated if the solution found is FEASIBLE.
 * @return maximum violations (Zero if feasible).
 */
double nlopt_optim_lazy_run(double** ballpred,
		              coptim *coparams,
	                  racketdes *racketdata,
		              optim *params) {

	static double qwait[NDOF];
	init_right_posture(qwait);
	lazy_data data = {racketdata,coparams,ballpred,qwait,
			          racketdata->dt,racketdata->Nmax};
	static int count = 0;
	printf("==========================================\n");
	printf("Running NLOPT\n");
	//initTime = get_time();
	//print_input_structs(&data,params);
	nlopt_optim_fixed_run(data.coparams,data.racketdata,params);
	/*double x[OPTIM_DIM];
	for (int i = 0; i < NDOF; i++) {
		x[i] = params->qf[i];
		x[i+NDOF] = params->qfdot[i];
	}
	x[2*NDOF] = params->T;
	double maxviol = test_optim(x,&data,true);*/
	double maxviol = nlopt_optim_lazy(&data,params);
	printf("Optim count: %d\n", (++count));
	printf("==========================================\n");
	return maxviol;
}

/*
 * Print input structs to give info about the arguments
 * For debugging purposes useful
 *
 */
static void print_input_structs(lazy_data* data, optim* params) {

	for (int i = 0; i < NDOF; i++) {
		printf("qwait[%d] = %f\n", i, data->qwait[i]);
		printf("q0[%d] = %f\n", i, data->coparams->q0[i]);
		printf("q0dot[%d] = %f\n", i, data->coparams->q0dot[i]);
		printf("lb[%d] = %f\n", i, data->coparams->lb[i]);
		printf("ub[%d] = %f\n", i, data->coparams->ub[i]);
	}
	printf("Tret = %f\n", data->coparams->time2return);
	for (int i = 0; i < NDOF; i++) {
		printf("qf[%d] = %f\n", i, params->qf[i]);
		printf("qfdot[%d] = %f\n", i, params->qfdot[i]);
	}
	printf("Thit = %f\n", params->T);

	print_mat_size("ball = ", data->ballpred, 2*NCART, 5);
	print_mat_size("pos = ", data->racketdata->pos, NCART, 5);
	print_mat_size("vel = ", data->racketdata->vel, NCART, 5);
	print_mat_size("normal = ", data->racketdata->normal, NCART, 5);

}

/*
 * NLOPT optimization routine for table tennis LAZY PLAYER (LP)
 *
 * If constraints are violated, it will not modify the lookup values (input)
 *
 * TODO: try different optimization routines
 *
 */
static double nlopt_optim_lazy(lazy_data *data, optim *params) {

	// firsttime checking
	firsttime[0] = firsttime[1] = firsttime[2] = true;

	//print_input_structs(coparams, racketdata, params);
	params->update = false;
	params->running = true;
	nlopt_opt opt;
	double x[OPTIM_DIM];
	double tol_ineq_land[INEQ_LAND_CONSTR_DIM];
	double tol_ineq_joint[INEQ_JOINT_CONSTR_DIM];
	const_vec(INEQ_LAND_CONSTR_DIM,1e-3,tol_ineq_land);
	const_vec(INEQ_JOINT_CONSTR_DIM,1e-3,tol_ineq_joint);
	init_soln(params,x); //parameters are the initial joint positions q0*/
	// set tolerances equal to second argument //

	// LN = does not require gradients //
	opt = nlopt_create(NLOPT_LN_COBYLA, OPTIM_DIM);
	nlopt_set_lower_bounds(opt, data->coparams->lb);
	nlopt_set_upper_bounds(opt, data->coparams->ub);
	nlopt_set_min_objective(opt, costfunc, data);
	nlopt_add_inequality_mconstraint(opt, INEQ_LAND_CONSTR_DIM,
			land_ineq_constr, data, tol_ineq_land);
	nlopt_add_inequality_mconstraint(opt, INEQ_JOINT_CONSTR_DIM,
			joint_limits_ineq_constr, data, tol_ineq_joint);
	nlopt_set_xtol_rel(opt, 1e-2);

	double init_time = get_time();
	double past_time = 0.0;
	double minf; // the minimum objective value, upon return
	int res; // error code
	double max_violation;

	if ((res = nlopt_optimize(opt, x, &minf)) < 0) {
		if (params->verbose)
			printf("NLOPT failed with exit code %d!\n", res);
	    past_time = (get_time() - init_time)/1e3;
	    max_violation = test_optim(x,data,params->update);
	}
	else {
		past_time = (get_time() - init_time)/1e3;
		if (params->verbose) {
			printf("NLOPT success with exit code %d!\n", res);
			printf("NLOPT took %f ms\n", past_time);
			printf("Found minimum at f = %0.10g\n", minf);
		}
	    max_violation = test_optim(x,data,params->update);
	    if (max_violation < 1e-2)
	    	finalize_soln(x,data,params,past_time);
	}
	params->running = false;
	//nlopt_destroy(opt);
	return max_violation;
}

/*
 * Calculates the cost function for table tennis Lazy Player (LP)
 * to find spline (3rd order strike+return) polynomials
 */
static double costfunc(unsigned n, const double *x, double *grad, void *my_func_params) {

	static double balldist;
	static double ballproj[NCART];
	static double qrest[NDOF];
	static double *q0dot; // initial joint velocity
	static double *q0; // initial joint pos
	static lazy_data * data;
	static weights *w = (weights*)malloc(sizeof(weights));
	static double J1, J2, Jhit, Jland, Jwait;
	static double a1[NDOF], a2[NDOF], b1[NDOF], b2[NDOF];
	static double T2, netTime, landTime, xnet;
	static double xland[2];
	double T = x[2*NDOF];

	if (firsttime[0]) {
		firsttime[0] = false;
		set_penalty_matrices(w);
		data = (lazy_data*)my_func_params;
		q0 = data->coparams->q0;
		q0dot = data->coparams->q0dot;
	}

	// calculate the polynomial coeffs which are used in the cost calculation
	calc_strike_poly_coeff(q0,q0dot,x,a1,a2);

	// calculate the landing time
	calc_times(data, x, &netTime, &landTime, &xnet, xland, &balldist, ballproj);

	// calculate the follow up coeffs which are added to the cost
	calc_return_poly_coeff(q0,x,landTime,w->R_wait,b1,b2);

	// calculate resting state
	calc_rest_robot(data,x,w->R_return,b1,b2,landTime,qrest);

	J1 = T * (3*T*T*inner_winv_prod(NDOF,w->R_strike,a1,a1) +
			3*T*inner_winv_prod(NDOF,w->R_strike,a1,a2) +
			inner_winv_prod(NDOF,w->R_strike,a2,a2));

	T2 = landTime;
	J2 = T2 * (3*T2*T2*inner_winv_prod(NDOF,w->R_return,b1,b1) +
			3*T2*inner_winv_prod(NDOF,w->R_return,b1,b2) +
			inner_winv_prod(NDOF,w->R_return,b2,b2));

	Jhit = inner_w_prod(NCART,w->R_hit,ballproj,ballproj);
	Jland = punish_land_robot(xland,xnet,w->R_land, w->R_net);
	Jwait = punish_wait_robot(data,w->R_wait,qrest);

	return J1 + J2 + Jhit + Jland + Jwait;
}

/*
 * Punish the robot sufficiently as to induce proper landing behaviour
 *
 * The desired landing location chosen is the center of the table
 */
static double punish_land_robot(const double *xland,
								const double xnet,
								const double *Rland,
								const double Rnet) {
	// desired landing locations
	static double x_des_land = 0.0;
	static double y_des_land = dist_to_table - 3*table_length/4;
	static double z_des_net = floor_level - table_height + net_height + 0.5;

	return sqr(xland[X] - x_des_land)*Rland[X] +
			sqr(xland[Y] - y_des_land)*Rland[Y] +
			sqr(xnet - z_des_net)*Rnet;

}

/*
 * Punish the deviations of resting joints (as a function of optim params) from a chosen init_joint_state
 *
 *
 */
static double punish_wait_robot(const lazy_data* data,
								double *Rwait,
								double *qrest) {
	static double qdiff[NDOF];

	for (int i = 0; i < NDOF; i++) {
		qdiff[i] = qrest[i] - data->qwait[i];
	}
	return inner_w_prod(NDOF,Rwait,qdiff,qdiff);
}

/*
 * Punish the deviations of resting joints (as a function of optim params) from a chosen init_joint_state
 *
 *
 */
static void calc_rest_robot(const lazy_data* data,
		                        const double *x,
		                        const double *Rret,
								const double *a3,
								const double *a2,
								const double landTime,
								double *qrest) {

	for (int i = 0; i < NDOF; i++) {
		qrest[i] = (x[i] + x[i+NDOF] * landTime +
				    a2[i] * pow(landTime,2) +
					a3[i] * pow(landTime,3)) / (2*Rret[i]);
	}
}

/*
 * This is the constraint that makes sure we land the ball
 *
 */
static void land_ineq_constr(unsigned m, double *result, unsigned n, const double *x, double *grad,
		                     void *my_func_params) {

	static lazy_data *data;
	static double xnet;
	static double xland[2];
	static double balldist;
	static double ballproj[NCART];
	static double netTime;
	static double landTime;
	static double table_xmax = table_width/2.0;
	static double table_ymax = dist_to_table - table_length;
	static double wall_z = 1.0;
	static double net_y = dist_to_table - table_length/2.0;
	static double net_z = floor_level - table_height + net_height;

	/* initialization of static variables */
	if (firsttime[1]) {
		firsttime[1] = false;
		data = (lazy_data*)my_func_params;
	}

	calc_times(data, x, &netTime, &landTime, &xnet, xland, &balldist, ballproj);

	result[0] = -balldist;
	result[1] = balldist - ball_radius;
	result[2] = inner_prod(NCART,ballproj,ballproj) - sqr(racket_radius);
	result[3] = -netTime;
	result[4] = xnet - wall_z;
	result[5] = -xnet + net_z;
	result[6] = netTime - landTime;
	result[7] = xland[X] - table_xmax;
	result[8] = -xland[X] - table_xmax;
	result[9] = xland[Y] - net_y;
	result[10] = -xland[Y] + table_ymax;
}

/*
 * This is the inequality constraint that makes sure we never exceed the
 * joint limits during the striking and returning motion
 *
 *
 */
static void joint_limits_ineq_constr(unsigned m, double *result,
		unsigned n, const double *x, double *grad, void *my_func_params) {

	static lazy_data* data;
	static weights *w = (weights *)malloc(sizeof(weights));
	static double a1[NDOF], a2[NDOF], a1ret[NDOF], a2ret[NDOF]; // coefficients for the polynomials
	static double joint_strike_max_cand[NDOF];
	static double joint_strike_min_cand[NDOF];
	static double joint_return_max_cand[NDOF];
	static double joint_return_min_cand[NDOF];
	static double *q0, *q0dot, *qwait;
	static double *ub, *lb;
	static double xnet, xland[2];
	static double balldist;
	static double ballproj[NCART];
	static double netTime, landTime;

	if (firsttime[2]) {
		firsttime[2] = false;
		set_penalty_matrices(w);
		data = (lazy_data*)my_func_params;
		q0 = data->coparams->q0;
		q0dot = data->coparams->q0dot;
		qwait = data->qwait;
		ub = data->coparams->ub;
		lb = data->coparams->lb;
	}

	calc_times(data, x, &netTime, &landTime, &xnet, xland, &balldist, ballproj);

	// calculate the polynomial coeffs which are used for checking joint limits
	calc_strike_poly_coeff(q0,q0dot,x,a1,a2);
	calc_return_poly_coeff(qwait,x,landTime,w->R_wait,a1ret,a2ret);
	// calculate the candidate extrema both for strike and return
	calc_strike_extrema_cand(a1,a2,x[2*NDOF],q0,q0dot,
			joint_strike_max_cand,joint_strike_min_cand);
	calc_return_extrema_cand(a1ret,a2ret,x,landTime,joint_return_max_cand,joint_return_min_cand);

	/* deviations from joint min and max */
	for (int i = 0; i < NDOF; i++) {
		result[i] = joint_strike_max_cand[i] - ub[i];
		result[i+NDOF] = lb[i] - joint_strike_min_cand[i];
		result[i+2*NDOF] = joint_return_max_cand[i] - ub[i];
		result[i+3*NDOF] = lb[i] - joint_return_min_cand[i];
		//printf("%f %f %f %f\n", result[i],result[i+DOF],result[i+2*DOF],result[i+3*DOF]);
	}
}


/*
 * Calculate the polynomial coefficients from the optimized variables qf,qfdot,T
 * p(t) = a1*t^3 + a2*t^2 + a3*t + a4
 *
 *
 */
static void calc_strike_poly_coeff(const double *q0, const double *q0dot, const double *x,
		                    double *a1, double *a2) {

	double T = x[2*NDOF];

	for (int i = 0; i < NDOF; i++) {
		a1[i] = (2/pow(T,3))*(q0[i]-x[i]) + (1/(T*T))*(q0dot[i] + x[i+NDOF]);
		a2[i] = (3/(T*T))*(x[i]-q0[i]) - (1/T)*(x[i+NDOF] + 2*q0dot[i]);
	}
	return;
}

/*
 * Calculate the returning polynomial coefficients from the optimized variables qf,qfdot
 * and time to return constant T
 * p(t) = a1*t^3 + a2*t^2 + a3*t + a4
 *
 * TODO: is this correct?
 *
 */
static void calc_return_poly_coeff(const double *qwait,
		                           const double *x, const double T,
		                           const double *Rwait,
								   double *a1, double *a2) {
	double weights[NDOF];
	for (int i = 0; i < NDOF; i++) {
		weights[i] = Rwait[i] / (6 + pow(T,3)*Rwait[i]/2);
		a1[i] = weights[i] * (x[i+NDOF]*T + 2*x[i] - 2*qwait[i]);
		a2[i] = -x[i+NDOF]/(2*T) - (3/2)*T*a1[i];
	}
}

/*
 * Calculate the extrema candidates for each joint (2*7 candidates in total)
 * For the striking polynomial
 * Clamp to [0,T]
 *
 */
static void calc_strike_extrema_cand(const double *a1, const double *a2, const double T,
		                      const double *q0, const double *q0dot,
		                      double *joint_max_cand, double *joint_min_cand) {

	static double cand1, cand2;

	for (int i = 0; i < NDOF; i++) {
		cand1 = fmin(T,fmax(0,(-a2[i] + sqrt(a2[i]*a2[i] - 3*a1[i]*q0dot[i]))/(3*a1[i])));
		cand2 =  fmin(T,fmax(0,(-a2[i] - sqrt(a2[i]*a2[i] - 3*a1[i]*q0dot[i]))/(3*a1[i])));
		cand1 = a1[i]*pow(cand1,3) + a2[i]*pow(cand1,2) + q0dot[i]*cand1 + q0[i];
		cand2 = a1[i]*pow(cand2,3) + a2[i]*pow(cand2,2) + q0dot[i]*cand2 + q0[i];
		joint_max_cand[i] = fmax(cand1,cand2);
		joint_min_cand[i] = fmin(cand1,cand2);
	}
}

/*
 * Calculate the extrema candidates for each joint (2*7 candidates in total)
 * For the return polynomial
 * Clamp to [0,TIME2RETURN]
 *
 */
static void calc_return_extrema_cand(double *a1, double *a2, const double *x, double landTime,
		                     double *joint_max_cand, double *joint_min_cand) {

	int i;
	static double cand1, cand2;
	double Tret = landTime;

	for (i = 0; i < NDOF; i++) {
		// find the extrema time candidates
		cand1 = fmin(Tret, fmax(0,(-a2[i] + sqrt(a2[i]*a2[i] - 3*a1[i]*x[i+NDOF]))/(3*a1[i])));
		cand2 = fmin(Tret, fmax(0,(-a2[i] - sqrt(a2[i]*a2[i] - 3*a1[i]*x[i+NDOF]))/(3*a1[i])));
		// find the joint extrema values at those candidate points
		cand1 = a1[i]*pow(cand1,3) + a2[i]*pow(cand1,2) + x[i+NDOF]*cand1 + x[i];
		cand2 = a1[i]*pow(cand2,3) + a2[i]*pow(cand2,2) + x[i+NDOF]*cand2 + x[i];
		joint_max_cand[i] = fmax(cand1,cand2);
		joint_min_cand[i] = fmin(cand1,cand2);
	}
}

/*
 * Calculate the net hitting time and the landing time
 * Assuming that there was an interaction (hitting) between racket and ball.
 *
 * Since we call this function multiple times within one optimization step
 * we store the optimization variable x,
 * landTime and netTime variables and return these when the optimization
 * variable is the same (meaning optimization algorithm did not update it yet).
 *
 * TODO: simplify code! If landtime is not real, should we set it to 1.0?
 *
 */
static void calc_times(const lazy_data* data,
		               const double *x,
					   double *netTime,
					   double *landTime,
					   double *xnet,
					   double *xland,
					   double *ball_norm, // ball to racket normal distance
					   double *ball_proj) { // ball projected to racket plane

	static double x_[OPTIM_DIM];
	static double b2rn_;
	static double b2rp_[NCART];
	static double xnet_;
	static double xland_[2];
	static double landTime_;
	static double netTime_; // local static variables for saving state

	static double qf[NDOF];
	static double qfdot[NDOF];
	static double vel[NCART];
	static double normal[NCART];
	static double pos[NCART];
	static double ballpos[NCART];
	static double ballvel[NCART];
	static double table_z = floor_level - table_height + ball_radius;
	static double g = fabs(gravity);
	static double net_y = dist_to_table - table_length/2;
	double distBall2TableZ;

	if (vec_is_equal(OPTIM_DIM,x,x_)) {
		*netTime = netTime_;
		*landTime = landTime_;
		*xnet = xnet_;
		make_equal(2,xland_,xland);
		*ball_norm = b2rn_;
		make_equal(NCART,b2rp_,ball_proj);
	}
	else {
		// extract state information from optimization variables
		for (int i = 0; i < NDOF; i++) {
			qf[i] = x[i];
			qfdot[i] = x[i + NDOF];
		}
		calc_racket_state(qf, qfdot, pos, vel, normal);
		interp_ball(data->ballpred, x[2*NDOF], data->dt, data->Nmax, ballpos, ballvel);
		racket_contact_model(vel, normal, ballvel);
		modify_ball_outgoing_vel(ballvel);
		distBall2TableZ = ballpos[Z] - table_z;

		if (sqr(ballvel[3]) > -2*g*distBall2TableZ) {
			landTime_ = *landTime = fmax(ballvel[Z] +
					sqrt(sqr(ballvel[Z]) + 2*g*distBall2TableZ),
				ballvel[Z] - sqrt(sqr(ballvel[Z]) + 2*g*distBall2TableZ)) / g;
		}
		else {
			// landTime is not real!
			landTime_ = *landTime = 1.0;
		}
		netTime_ = *netTime = (net_y - ballpos[Y])/ballvel[Y];

		xnet_ = *xnet = ballpos[Z] + netTime_ * ballvel[Z] +
				        0.5*gravity*netTime_*netTime_;
		xland_[0] = xland[0] = ballpos[X] + landTime_ * ballvel[X];
		xland_[1] = xland[1] = ballpos[Y] + landTime_ * ballvel[Y];

		// calculate deviation of ball to racket - hitting constraints
		make_equal(NCART,ballpos,ball_proj);
		vec_minus(NCART,pos,ball_proj);
		b2rn_ = *ball_norm = inner_prod(NCART,normal,ball_proj);
		for (int i = 0; i < NCART; i++) {
			normal[i] *= b2rn_;
		}
		vec_minus(NCART,normal,ball_proj); // we get (e - nn'e) where e is ballpos - racketpos
		make_equal(NCART,ball_proj,b2rp_);
		make_equal(OPTIM_DIM,x,x_);
	}
}

/*
 * Debug by testing the cost function value and
 * the constraint violation of the solution vector x
 *
 * Returns FALSE if joint inequality constraints are violated!
 *
 *
 */
static double test_optim(double *x, lazy_data* data, bool verbose) {

	// give info on constraint violation
	double *grad = 0;
	static double land_violation[INEQ_LAND_CONSTR_DIM];
	static double lim_violation[INEQ_JOINT_CONSTR_DIM]; // joint limit violations on strike and return
	joint_limits_ineq_constr(INEQ_JOINT_CONSTR_DIM, lim_violation, OPTIM_DIM, x, grad, data);
	land_ineq_constr(INEQ_LAND_CONSTR_DIM, land_violation, OPTIM_DIM, x, grad, data);
	double cost = costfunc(OPTIM_DIM, x, grad, data);

	if (verbose) {
		// give info on solution vector
		print_optim_vec(x);
		printf("f = %.2f\n",cost);
		printf("Constraint info:\n");
		printf("Ball to racket normal distance: %.2f\n",land_violation[0]);
		printf("Ball distance projected to racket: %.2f\n", land_violation[2]);
		printf("netTime: %f\n", -land_violation[3]);
	    printf("landTime - netTime: %f\n", -land_violation[6]);
		printf("wall_z - zNet: %f\n", -land_violation[4]);
		printf("zNet - net_z: %f\n", -land_violation[5]);
		printf("table_xmax - xLand: %f\n", -land_violation[6]);
		printf("table_xmax + xLand: %f\n", -land_violation[8]);
		printf("net_y - yLand: %f\n", -land_violation[9]);
	    printf("yLand - table_ymax: %f\n", -land_violation[10]);
		for (int i = 0; i < INEQ_JOINT_CONSTR_DIM; i++) {
			if (lim_violation[i] > 0.0) {
				printf("Joint limit violated by %.2f on joint %d\n", lim_violation[i], i % NDOF + 1);
			}
		}
	}

	return fmax(max_array(lim_violation,INEQ_JOINT_CONSTR_DIM),
			max_array(land_violation,INEQ_LAND_CONSTR_DIM));
}

/*
 * Modify ball outgoing velocity by multiplying with a constant vector (less than one)
 *
 * Necessary since the landing and net hitting locations are calculated using projectile motion
 *
 */
static void modify_ball_outgoing_vel(double* ballVel) {

	static double multiplier_x, multiplier_y, multiplier_z;
	static bool firsttime = true;

	if (firsttime) {
		firsttime = false;
		multiplier_x = 0.9;
		multiplier_y = 0.8;
		multiplier_z = 0.83;
	}

	ballVel[X] *= multiplier_x;
	ballVel[Y] *= multiplier_y;
	ballVel[Z] *= multiplier_z;
}

/*
 * Initialize the penalty matrices used to punish the robot
 *
 * Penalty structure needs to be initialized before
 *
 */
static void set_penalty_matrices(weights * pen) {

	double* R1 = (double*)calloc(NDOF,sizeof(double));
	double* R2 = (double*)calloc(NDOF,sizeof(double));
	double* Rwait = (double*)calloc(NDOF,sizeof(double));
	double* Rhit = (double*)calloc(NDOF,sizeof(double));
	double* Rland = (double*)calloc(2,sizeof(double));
	double Rnet = 1.0;

	const_vec(NDOF,1.0,R1);
	const_vec(NDOF,1.0,R2);
	const_vec(NDOF,1.0,Rwait);
	const_vec(NCART,1.0,Rhit);
	const_vec(2,1.0,Rland);

	pen->R_hit = Rhit;
	pen->R_land = Rland;
	pen->R_return = R2;
	pen->R_strike = R1;
	pen->R_wait = Rwait;
	pen->R_net = Rnet;
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
 *
 */
static void racket_contact_model(double* racketVel, double* racketNormal, double* ballVel) {

	static double diffVel[NCART];
	static double normalMultSpeed[NCART];
	double speed;

	for (int i = 0; i < NCART; i++)
		diffVel[i] = racketVel[i] - ballVel[i];

	speed = (1 + CRR) * inner_prod(NCART, racketNormal, diffVel);

	for (int i = 0; i < NCART; i++) {
		normalMultSpeed[i] = speed * racketNormal[i];
	}
	vec_plus(NCART,normalMultSpeed,ballVel);
}

/*
 * First order hold to interpolate linearly at time T
 * between racket pos,vel,normal entries
 *
 * IF T is nan, racket variables are assigned to zero-element of
 * relevant racket entries
 *
 */
static void interp_ball(double** ballpred, const double T,
		                const double dt, const int Nmax,
				        double *ballpos, double *ballvel) {
	int N = (int) (T/dt);
	double Tdiff = T - N*dt;

	for (int i = 0; i < NCART; i++) {
		if (N < Nmax) {
			ballpos[i] = ballpred[i][N] +
					(Tdiff/dt) * (ballpred[i][N+1] - ballpred[i][N]);
			ballvel[i] = ballpred[i+NCART][N] +
					(Tdiff/dt) * (ballpred[i+NCART][N+1] - ballpred[i+NCART][N]);
		}
		else {
			ballpos[i] = ballpred[i][N];
			ballvel[i] = ballpred[i+NCART][N];
		}
	}
}

/*
 * Estimate an initial solution for NLOPT
 * 2*dof + 1 dimensional problem
 *
 * The closer to the optimum it is the faster alg should converge
 */
static void init_soln(const optim * params, double x[OPTIM_DIM]) {

	// initialize first dof entries to q0
	for (int i = 0; i < NDOF; i++) {
		x[i] = params->qf[i];
		x[i+NDOF] = params->qfdot[i];
	}
	x[2*NDOF] = params->T;
}

/*
 * Initialize robot posture on the right size of the robot
 */
static void init_right_posture(double* qwait) {

	qwait[0] = 1.0;
	qwait[1] = -0.2;
	qwait[2] = -0.1;
	qwait[3] = 1.8;
	qwait[4] = -1.57;
	qwait[5] = 0.1;
	qwait[6] = 0.3;
}

/*
 * Finalize the solution and update target optimization
 * parameter structure and hitTime values
 *
 *
 */
static void finalize_soln(const double* x,
		lazy_data* data, optim * params, double time_elapsed) {

	static weights *w = (weights*)malloc(sizeof(weights));
	static double netTime, landTime, xnet, balldist, xland[2], ballproj[3];
	static double *qrest = (double*)calloc(NDOF,sizeof(double));
	static double b1[NDOF], b2[NDOF];

	set_penalty_matrices(w);
	calc_times(data, x, &netTime, &landTime, &xnet, xland, &balldist, ballproj);

	calc_return_poly_coeff(data->qwait,x,landTime,w->R_wait,b1,b2);
	calc_rest_robot(data,x,w->R_wait,b1,b2,landTime,qrest);

	// initialize first dof entries to q0
	for (int i = 0; i < NDOF; i++) {
		params->qf[i] = x[i];
		params->qfdot[i] = x[i+NDOF];
	}
	params->T = x[2*NDOF] - (time_elapsed/1e3);
	data->coparams->time2return = landTime;
	data->coparams->qrest = qrest;
	params->update = true;
}
