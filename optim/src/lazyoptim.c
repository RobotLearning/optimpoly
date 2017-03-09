/*
 * lazyoptim.c
 *
 * NLOPT polynomial optimization functions for LAZY PLAYER are stored here.
 *
 * Lazy Player make use of them to generate trajectories.
 *
 *  Created on: Sep 11, 2016
 *      Author: okan
 */

#include "constants.h"
#include "table.h"
#include "utils.h"
#include "kinematics.h"
#include "stdlib.h"
#include "math.h"
#include "optim.h"

#define EQ_HIT_CONSTR_DIM 1
#define INEQ_LAND_CONSTR_DIM 9
#define INEQ_JOINT_CONSTR_DIM 2*NDOF + 2*NDOF

static double nlopt_optim_lazy(ball_data *data, optim *params);
static void print_input_structs(ball_data* data, optim* params);
static double costfunc(unsigned n, const double *x, double *grad, void *my_func_params);
static double punish_hit_robot(const ball_data * data,
		                       const double *x,
							   const double *Rhit);
static double punish_land_robot(const ball_data * data,
		                        const double *x,
								const double *Rland,
								const double Rnet);
static double punish_wait_robot(const ball_data* data,
		                        const double *x,
		                        const double *Rwait,
								double landTime);
static void land_ineq_constr(unsigned m, double *result, unsigned n, const double *x, double *grad,
		                     void *my_func_params);
static void hit_eq_constr(unsigned m, double *result, unsigned n, const double *x, double *grad,
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
static double test_optim(double *x, ball_data* data, int verbose);
static void modify_ball_outgoing_vel(double* ballVel);
static void calc_times(const ball_data* data, const double *x, double *netTime, double *landTime);
static void finalize_soln(const double* x, optim * params, double time_elapsed);
static void set_penalty_matrices(weights * pen);
static void racket_contact_model(double* racketVel, double* racketNormal, double* ballVel);
static void interp_ball(double** ballpred, const double T,
		                const double dt, const int Nmax,
						double *ballpos, double *ballvel);
static void init_soln(const optim * params, double x[OPTIM_DIM]);


/*
 * Multi-threading entry point for the NLOPT optimization
 */
double optim_lazy_run(double** ballpred,
		              coptim *coparams,
	                  racketdes *racketdata,
		              optim *params) {

	ball_data data = {racketdata,coparams,ballpred,racketdata->dt,racketdata->Nmax};
	static int count = 0;
	printf("==========================================\n");
	printf("Running NLOPT\n");
	//initTime = get_time();
	//print_input_structs(data,params);
	nlopt_optim_fixed_run(data.coparams,data.racketdata,params);
	/*double x[OPTIM_DIM];
	for (int i = 0; i < NDOF; i++) {
		x[i] = params->qf[i];
		x[i+NDOF] = params->qfdot[i];
	}
	x[2*NDOF] = params->T;
	double maxviol = test_optim(x,&data,TRUE);*/
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
static void print_input_structs(ball_data* data, optim* params) {

	for (int i = 0; i < NDOF; i++) {
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
 */
static double nlopt_optim_lazy(ball_data *data, optim *params) {

	//print_input_structs(coparams, racketdata, params);
	params->update = FALSE;
	params->running = TRUE;
	nlopt_opt opt;
	double x[OPTIM_DIM];
	double tol_eq_hit[EQ_HIT_CONSTR_DIM];
	double tol_ineq_land[INEQ_LAND_CONSTR_DIM];
	double tol_ineq_joint[INEQ_JOINT_CONSTR_DIM];
	const_vec(EQ_CONSTR_DIM,1e-2,tol_eq_hit);
	const_vec(INEQ_LAND_CONSTR_DIM,1e-3,tol_ineq_land);
	const_vec(INEQ_JOINT_CONSTR_DIM,1e-3,tol_ineq_joint);
	init_soln(params,x); //parameters are the initial joint positions q0*/
	// set tolerances equal to second argument //

	// LN = does not require gradients //
	opt = nlopt_create(NLOPT_LN_COBYLA, OPTIM_DIM);
	nlopt_set_xtol_rel(opt, 1e-2);
	nlopt_set_lower_bounds(opt, data->coparams->lb);
	nlopt_set_upper_bounds(opt, data->coparams->ub);
	nlopt_set_min_objective(opt, costfunc, data);
	nlopt_add_equality_mconstraint(opt, EQ_HIT_CONSTR_DIM, hit_eq_constr, data, tol_eq_hit);
	nlopt_add_inequality_mconstraint(opt, INEQ_LAND_CONSTR_DIM,
			land_ineq_constr, data, tol_ineq_land);
	nlopt_add_inequality_mconstraint(opt, INEQ_JOINT_CONSTR_DIM,
			joint_limits_ineq_constr, data, tol_ineq_joint);
	nlopt_set_xtol_rel(opt, 1e-2);

	double init_time = get_time();
	double past_time = 0.0;
	double minf; // the minimum objective value, upon return //
	int res; // error code
	double max_violation;

	if ((res = nlopt_optimize(opt, x, &minf)) < 0) {
	    printf("NLOPT failed with exit code %d!\n", res);
	    past_time = (get_time() - init_time)/1e3;
	    max_violation = test_optim(x,data,TRUE);
	}
	else {
		past_time = (get_time() - init_time)/1e3;
		printf("NLOPT success with exit code %d!\n", res);
		printf("NLOPT took %f ms\n", past_time);
	    printf("Found minimum at f = %0.10g\n", minf);
	    max_violation = test_optim(x,data,TRUE);
	    if (max_violation < 1e-2)
	    	finalize_soln(x,params,past_time);
	}
	params->running = FALSE;
	//nlopt_destroy(opt);
	return max_violation;
}

/*
 * Calculates the cost function for table tennis Lazy Player (LP)
 * to find spline (3rd order strike+return) polynomials
 */
static double costfunc(unsigned n, const double *x, double *grad, void *my_func_params) {

	static double *q0dot; // initial joint velocity
	static double *q0; // initial joint pos
	static int firsttime = TRUE;
	static ball_data * data;
	static weights *w = (weights*)malloc(sizeof(weights));
	static double J1, J2, Jhit, Jland, Jwait;
	static double a1[NDOF];
	static double a2[NDOF];
	static double b1[NDOF];
	static double b2[NDOF];
	static double T2, netTime, landTime;
	double T = x[2*NDOF];

	if (firsttime) {
		firsttime = FALSE;
		set_penalty_matrices(w);
		data = (ball_data*)my_func_params;
		q0 = data->coparams->q0;
		q0dot = data->coparams->q0dot;
	}

	// calculate the polynomial coeffs which are used in the cost calculation
	calc_strike_poly_coeff(q0,q0dot,x,a1,a2);

	// calculate the landing time
	calc_times(data, x, &netTime, &landTime);

	// calculate the follow up coeffs which are added to the cost
	calc_return_poly_coeff(q0,x,landTime,w->R_wait,b1,b2);

	J1 = T * (3*T*T*inner_winv_prod(NDOF,w->R_strike,a1,a1) +
			3*T*inner_winv_prod(NDOF,w->R_strike,a1,a2) +
			inner_winv_prod(NDOF,w->R_strike,a2,a2));

	T2 = landTime;
	J2 = T2 * (3*T2*T2*inner_winv_prod(NDOF,w->R_return,b1,b1) +
			3*T2*inner_winv_prod(NDOF,w->R_return,b1,b2) +
			inner_winv_prod(NDOF,w->R_return,b2,b2));

	Jhit = punish_hit_robot(data,x,w->R_hit);
	Jland = punish_land_robot(data,x,w->R_land, w->R_net);
	Jwait = punish_wait_robot(data,x,w->R_wait,landTime);

	return J1 + Jhit + Jland + Jwait + J2;
}

/*
 * Punish the robot to encourage hitting
 *
 * We punish the square of the deviations of the ball projected to the racket to the racket center
 */
static double punish_hit_robot(const ball_data * data,
		                       const double *x,
							   const double *Rhit) {
	static double pos[NCART];
	static double vel[NCART];
	static double normal[NCART];
	static double ball_pos[NCART];
	static double ball_vel[NCART];
	static double qf[NDOF];
	static double qfdot[NDOF];
	double ball2rob_along_n;
	double T = x[2*NDOF];

	// extract state information from optimization variables
	for (int i = 0; i < NDOF; i++) {
		qf[i] = x[i];
		qfdot[i] = x[i + NDOF];
	}

	interp_ball(data->ballpred, T, data->dt, data->Nmax, ball_pos, ball_vel);
	// calculate deviation of ball to racket - hitting constraints
	calc_racket_state(qf,qfdot,pos,vel,normal);
	vec_minus(NCART,pos,ball_pos);
	ball2rob_along_n = inner_prod(NCART,normal,ball_pos);
	for (int i = 0; i < NCART; i++) {
		normal[i] *= ball2rob_along_n;
	}
	vec_minus(NCART,normal,ball_pos);

	return inner_w_prod(NCART,Rhit,ball_pos,ball_pos);
}

/*
 * Punish the robot sufficiently as to induce proper landing behaviour
 *
 * The desired landing location chosen is the center of the table
 */
static double punish_land_robot(const ball_data * data,
		                        const double *x,
								const double *Rland,
								const double Rnet) {
	static double pos[NCART];
	static double vel[NCART];
	static double normal[NCART];
	static int firsttime = TRUE;
	static double ball_pos[NCART];
	static double ball_vel[NCART];
	static double qf[NDOF];
	static double qfdot[NDOF];
	static double netTime;
	static double landTime;
	// desired landing locations
	static double x_des_land;
	static double y_des_land;
	static double z_des_net;
	double z_net;
	double x_land;
	double y_land;
	double T = x[2*NDOF];

	/* initialization of static variables */
	if (firsttime) {
		firsttime = FALSE;
		x_des_land = 0.0;
		y_des_land = dist_to_table - 3*table_length/4;
		z_des_net = floor_level - table_height + net_height + 0.5;
	}

	// extract state information from optimization variables
	for (int i = 0; i < NDOF; i++) {
		qf[i] = x[i];
		qfdot[i] = x[i+NDOF];
	}

	interp_ball(data->ballpred, T, data->dt, data->Nmax, ball_pos, ball_vel);
	// calculate deviation of ball to racket - hitting constraints
	calc_racket_state(qf,qfdot,pos,vel,normal);
	calc_times(data, x, &netTime, &landTime);

	racket_contact_model(vel, normal, ball_vel);
	modify_ball_outgoing_vel(ball_vel);
	z_net = ball_pos[Z] + netTime * ball_vel[Z] + 0.5*gravity*netTime*netTime;
	x_land = ball_pos[X] + landTime * ball_vel[X];
	y_land = ball_pos[Y] + landTime * ball_vel[Y];

	return sqr(x_land - x_des_land)*Rland[X] +
			sqr(y_land - y_des_land)*Rland[Y] +
			sqr(z_net - z_des_net)*Rnet;

}

/*
 * Punish the deviations of resting joints (as a function of optim params) from a chosen init_joint_state
 */
static double punish_wait_robot(const ball_data* data,
		                        const double *x,
		                        const double *Rwait,
								double landTime) {

	int i;
	static int firsttime = TRUE;
	static double* qwait;
	static double qdiff[NDOF];
	static double qrest[NDOF];

	if (firsttime) {
		firsttime = FALSE;
		qwait = data->coparams->qrest;
	}

	for (i = 0; i < NDOF; i++) {
		qrest[i] = x[i] + x[i+NDOF] * landTime/2;
		qdiff[i] = qrest[i] - qwait[i];
	}

	return inner_w_prod(NDOF,Rwait,qdiff,qdiff);
}

/*
 * This is the constraint that makes sure we land the ball
 *
 */
static void land_ineq_constr(unsigned m, double *result, unsigned n, const double *x, double *grad,
		                     void *my_func_params) {

	static double qf[NDOF];
	static double qfdot[NDOF];
	static ball_data* data;
	static double racketvel[NCART];
	static double racketnormal[NCART];
	static double racketpos[NCART];
	static int firsttime = TRUE;
	static double netTime;
	static double landTime;
	static double table_xmax = table_width/2.0;
	static double table_ymax = dist_to_table - table_length;
	static double wall_z = 1.0;
	static double net_y = dist_to_table - table_length/2.0;
	static double net_z = floor_level - table_height + net_height;
	static double ballpos[NCART];
	static double ballvel[NCART];
	double ball2rob;
	double zNet;
	double xLand;
	double yLand;
	double T = x[2*NDOF];

	/* initialization of static variables */
	if (firsttime) {
		firsttime = FALSE;
		data = (ball_data*)my_func_params;
	}

	// extract state information from optimization variables
	for (int i = 0; i < NDOF; i++) {
		qf[i] = x[i];
		qfdot[i] = x[i + NDOF];
	}

	interp_ball(data->ballpred, T, data->dt, data->Nmax, ballpos, ballvel);
	calc_times(data, x, &netTime, &landTime);
	// calculate landing and net positions
	calc_racket_state(qf,qfdot,racketpos,racketvel,racketnormal);
	racket_contact_model(racketvel, racketnormal, ballvel);
	modify_ball_outgoing_vel(ballvel);
	zNet = ballpos[Z] + netTime * ballvel[Z] + 0.5*gravity*netTime*netTime;
	xLand = ballpos[X] + landTime * ballvel[X];
	yLand = ballpos[Y] + landTime * ballvel[Y];

	// calculate deviation of ball to racket - hitting constraints
	vec_minus(NCART,racketpos,ballpos);
	ball2rob = inner_prod(NCART,racketnormal,ballpos);
	for (int i = 0; i < NCART; i++) {
		racketnormal[i] *= ball2rob;
	}
	vec_minus(NCART,racketnormal,ballpos); // we get (e - nn'e) where e is ballpos - racketpos

	result[0] = inner_prod(NCART,ballpos,ballpos) - sqr(racket_radius);
	result[1] = -netTime;
	result[2] = zNet - wall_z;
	result[3] = -zNet + net_z;
	result[4] = netTime - landTime;
	result[5] = xLand - table_xmax;
	result[6] = -xLand - table_xmax;
	result[7] = yLand - net_y;
	result[8] = -yLand + table_ymax;

}

/*
 * Single equality constraint to make sure we hit the ball
 */
static void hit_eq_constr(unsigned m, double *result, unsigned n, const double *x, double *grad,
		                  void *my_func_params) {

	static int firsttime = TRUE;
	static double qf[NDOF];
	static double qfdot[NDOF];
	static double pos[NCART];
	static double vel[NCART];
	static double normal[NCART];
	static double ballpos[NCART];
	static double ballvel[NCART];
	static ball_data* data;
	double ball2rob;
	double T = x[2*NDOF];

	/* initialization of static variables */
	if (firsttime) {
		firsttime = FALSE;
		data = (ball_data*)my_func_params;
	}

	// extract state information from optimization variables
	for (int i = 0; i < NDOF; i++) {
		qf[i] = x[i];
		qfdot[i] = x[i + NDOF];
	}

	interp_ball(data->ballpred, T, data->dt, data->Nmax, ballpos, ballvel);
	// calculate deviation of ball to racket - hitting constraints
	calc_racket_state(qf,qfdot,pos,vel,normal);
	vec_minus(NCART,pos,ballpos);
	ball2rob = inner_prod(NCART,normal,ballpos);

	result[0] = ball2rob;
}

/*
 * This is the inequality constraint that makes sure we never exceed the
 * joint limits during the striking and returning motion
 *
 * FIXME: this is broken!
 *
 */
static void joint_limits_ineq_constr(unsigned m, double *result,
		unsigned n, const double *x, double *grad, void *my_func_params) {

	static int firsttime = TRUE;
	static ball_data* data;
	static weights *w = (weights *)malloc(sizeof(weights));
	static double a1[NDOF];
	static double a2[NDOF];
	static double a1ret[NDOF]; // coefficients for the returning polynomials
	static double a2ret[NDOF];
	static double joint_strike_max_cand[NDOF];
	static double joint_strike_min_cand[NDOF];
	static double joint_return_max_cand[NDOF];
	static double joint_return_min_cand[NDOF];
	static double *q0;
	static double *q0dot;
	static double *qrest;
	static double *ub;
	static double *lb;
	static double Tret;
	double netTime, landTime = 0.0;

	if (firsttime) {
		firsttime = FALSE;
		set_penalty_matrices(w);
		data = (ball_data*)my_func_params;
		q0 = data->coparams->q0;
		q0dot = data->coparams->q0dot;
		qrest = data->coparams->qrest;
		ub = data->coparams->ub;
		lb = data->coparams->lb;
	}

	calc_times(data, x, &netTime, &landTime);

	// calculate the polynomial coeffs which are used for checking joint limits
	calc_strike_poly_coeff(q0,q0dot,x,a1,a2);
	calc_return_poly_coeff(qrest,x,landTime,w->R_wait,a1ret,a2ret);
	// calculate the candidate extrema both for strike and return
	calc_strike_extrema_cand(a1,a2,x[2*NDOF],q0,q0dot,
			joint_strike_max_cand,joint_strike_min_cand);
	calc_return_extrema_cand(a1ret,a2ret,x,Tret,joint_return_max_cand,joint_return_min_cand);

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
 */
static void calc_return_poly_coeff(const double *q0,
		                           const double *x, const double T,
		                           const double *Rwait,
								   double *a1, double *a2) {
	double weights[NDOF];
	for (int i = 0; i < NDOF; i++) {
		weights[i] = Rwait[i] / (6 + pow(T,3)*Rwait[i]);
		a1[i] = weights[i] * (x[i+NDOF]*T + 2*x[i] - 2*q0[i]);
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
 * Debug by testing the cost function value and
 * the constraint violation of the solution vector x
 *
 * Returns FALSE if joint inequality constraints are violated!
 *
 */
static double test_optim(double *x, ball_data* data, int verbose) {

	// give info on constraint violation
	double *grad = FALSE;
	static double hit_violation[EQ_HIT_CONSTR_DIM];
	static double land_violation[INEQ_LAND_CONSTR_DIM];
	static double lim_violation[INEQ_JOINT_CONSTR_DIM]; // joint limit violations on strike and return
	joint_limits_ineq_constr(INEQ_JOINT_CONSTR_DIM, lim_violation, OPTIM_DIM, x, grad, data);
	hit_eq_constr(EQ_HIT_CONSTR_DIM, hit_violation, OPTIM_DIM, x, grad, data);
	land_ineq_constr(INEQ_LAND_CONSTR_DIM, land_violation, OPTIM_DIM, x, grad, data);
	double cost = costfunc(OPTIM_DIM, x, grad, data);

	if (verbose) {
		// give info on solution vector
		print_optim_vec(x);
		printf("f = %.2f\n",cost);
		printf("Hitting violations:\n");
		printf("Ball to racket normal distance: %.2f\n",hit_violation[0]);
		printf("Ball distance projected to racket: %.2f\n", land_violation[0]);
		//printf("Ball along racket: %.2f\n", hit_violation[1]);
		printf("Landing info:\n");
		printf("netTime: %f\n", -land_violation[1]);
	    printf("landTime - netTime: %f\n", -land_violation[4]);
		printf("wall_z - zNet: %f\n", -land_violation[2]);
		printf("zNet - net_z: %f\n", -land_violation[3]);
		printf("table_xmax - xLand: %f\n", -land_violation[5]);
		printf("table_xmax + xLand: %f\n", -land_violation[6]);
		printf("net_y - yLand: %f\n", -land_violation[7]);
	    printf("yLand - table_ymax: %f\n", -land_violation[8]);
	}

	for (int i = 0; i < INEQ_JOINT_CONSTR_DIM; i++) {
		if (lim_violation[i] > 0.0) {
			printf("Joint limit violated by %.2f on joint %d\n", lim_violation[i], i % NDOF + 1);
		}
	}

	return fmax(fmax(max_abs_array(hit_violation,EQ_HIT_CONSTR_DIM),
		    max_array(lim_violation,INEQ_JOINT_CONSTR_DIM)),
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
	static int firsttime = TRUE;

	if (firsttime) {
		firsttime = FALSE;
		multiplier_x = 0.9;
		multiplier_y = 0.8;
		multiplier_z = 0.83;
	}

	ballVel[X] *= multiplier_x;
	ballVel[Y] *= multiplier_y;
	ballVel[Z] *= multiplier_z;

}

/*
 * Calculate the net hitting time and the landing time
 * Assuming that there was an interaction (hitting) between racket and ball
 *
 */
static void calc_times(const ball_data* data, const double *x, double *netTime, double *landTime) {

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

	// extract state information from optimization variables
	for (int i = 0; i < NDOF; i++) {
		qf[i] = x[i];
		qfdot[i] = x[i + NDOF];
	}

	calc_racket_state(qf, qfdot,pos, vel, normal);
	interp_ball(data->ballpred, x[2*NDOF], data->dt, data->Nmax, ballpos, ballvel);
	racket_contact_model(vel, normal, ballvel);
	modify_ball_outgoing_vel(ballvel);
	distBall2TableZ = ballpos[Z] - table_z;

	if (sqr(ballvel[3]) > -2*g*distBall2TableZ) {
		*landTime = fmax(ballvel[Z] + sqrt(sqr(ballvel[Z]) + 2*g*distBall2TableZ),
		                ballvel[Z] - sqrt(sqr(ballvel[Z]) + 2*g*distBall2TableZ)) / g;
	}
	else {
		// landTime is not real!
		*landTime = 1.0;
	}

	*netTime = (net_y - ballpos[Y])/ballvel[Y];
}

/*
 * Finalize the solution and update target SL structure and hitTime value
 */
static void finalize_soln(const double* x, optim * params, double time_elapsed) {

	// initialize first dof entries to q0
	for (int i = 0; i < NDOF; i++) {
		params->qf[i] = x[i];
		params->qfdot[i] = x[i+NDOF];
	}
	params->T = x[2*NDOF] - (time_elapsed/1e3);
	params->update = TRUE;
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
	const_vec(NCART,10.0,Rhit);
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

	if (isnan(T)) {
		printf("Warning: T value is nan!\n");

		for(int i = 0; i < NCART; i++) {
			ballpos[i] = ballpred[i][0];
			ballvel[i] = ballpred[i+NCART][0];
		}
	}
	else {
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
