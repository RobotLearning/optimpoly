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
#include "utils.h"
#include "kinematics.h"
#include "stdlib.h"
#include "math.h"
#include "optim.h"

#define EQ_HIT_CONSTR_DIM 1
#define INEQ_LAND_CONSTR_DIM 9
#define INEQ_JOINT_CONSTR_DIM 2*NDOF + 2*NDOF

typedef struct {
	coptim *coparams;
	double **ballpred;
} lazy_data;

/*
 * Multi-threading entry point for the NLOPT optimization
 */
double optim_lazy_run(coptim *coparams,
	                     racketdes *racketdata,
		                 optim *params) {

	static int count = 0;
	printf("==========================================\n");
	printf("Running NLOPT\n");
	//initTime = get_time();
	nlopt_optim_fixed_run(coparams,racketdata,params);
	lazy_data data = {coparams,racketdata->pos};
	double maxviol = nlopt_optim_lazy(&data,params);
	printf("Optim count: %d\n", (++count));
	printf("==========================================\n");
	return maxviol;
}

/*
 * NLOPT optimization routine for table tennis LAZY PLAYER (LP)
 *
 * If constraints are violated, it will not modify the lookup values (input)
 */
double nlopt_optim_lazy(lazy_data *data, optim *params) {

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
	nlopt_add_equality_mconstraint(opt, EQ_HIT_CONSTR_DIM, hit_eq_constr, NULL, tol_eq_hit);
	nlopt_add_inequality_mconstraint(opt, INEQ_LAND_CONSTR_DIM,
			land_ineq_constr, data, tol_ineq_land);
	nlopt_add_inequality_mconstraint(opt, INEQ_JOINT_CONSTR_DIM,
			joint_lim_ineq_constr, data->coparams, tol_ineq_joint);
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
	check_optim_result(res);
	//nlopt_destroy(opt);
	return max_violation;
}

/*
 * Calculates the cost function for table tennis Lazy Player (LP)
 * to find spline (3rd order strike+return) polynomials
 */
static double costfunc(unsigned n, double *x, double *grad, void *my_func_params) {

	static double **ballpred;
	static double *q0dot; // initial joint velocity
	static double *q0; // initial joint pos
	static int firsttime = TRUE;
	static weights * w;
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
		lazy_data* data = (lazy_data*)my_func_params;
		q0 = data->coparams->q0;
		q0dot = data->coparams->q0dot;
		ballpred = data->ballpred;
	}

	// calculate the polynomial coeffs which are used in the cost calculation
	calc_strike_poly_coeff(q0,q0dot,x,a1,a2);

	// calculate the landing time
	calc_times(x, &netTime, &landTime);

	// calculate the follow up coeffs which are added to the cost
	calc_return_poly_coeff(q0,x,landTime,w->R_wait,b1,b2);

	J1 = T * (3*T*T*inner_winv_prod(a1,a1,w->R_strike,NDOF) +
			3*T*inner_winv_prod(a1,a2,w->R_strike,NDOF) +
			inner_winv_prod(a2,a2,w->R_strike,NDOF));

	T2 = landTime;
	J2 = T2 * (3*T2*T2*inner_winv_prod(b1,b1,w->R_return,NDOF) +
			3*T2*inner_winv_prod(b1,b2,w->R_return,NDOF) +
			inner_winv_prod(b2,b2,w->R_return,NDOF));

	Jhit = punish_hit_robot(ballpred,x,w->R_hit);
	Jland = punish_land_robot(ballpred,x,w->R_land, w->R_net);
	Jwait = punish_wait_robot(x,w->R_wait,landTime);

	return J1 + Jhit + Jland + Jwait; //J2
}

/*
 * Punish the robot to encourage hitting
 *
 * We punish the square of the deviations of the ball projected to the racket to the racket center
 */
static double punish_hit_robot(const double** ballpred,
		                       const double *x,
							   const double *Rhit) {

	int i;
	static Vector  xfdot;
	static Vector  normal;
	static Vector  pos;
	static int firsttime = TRUE;
	static SL_Cstate ballPredict;
	static Vector vecBallAlongRacket;

	double distB2RAlongN;
	double T = x[2*NDOF];

	/* initialization of static variables */
	if (firsttime) {
		firsttime = FALSE;
		xfdot = my_vector(1,2*NCART);
		normal = my_vector(1,NCART);
		pos = my_vector(1,NCART);
		vecBallAlongRacket = my_vector(1,NCART);
		bzero((char *)&ballPredict, sizeof(ballPredict));
	}

	interp_ball(&ballPredict, T);
	for (i = 1; i <= NCART; i++) {
		vecBallAlongRacket[i] = ballPredict.x[i];
	}

	// calculate deviation of ball to racket - hitting constraints
	calc_racket_state(x,pos,xfdot,normal);
	vec_sub(vecBallAlongRacket,pos,vecBallAlongRacket);
	distB2RAlongN = vec_mult_inner(normal,vecBallAlongRacket);
	vec_mult_scalar(normal,distB2RAlongN,normal);
	vec_sub(vecBallAlongRacket,normal,vecBallAlongRacket);

	return inner_w_prod(vecBallAlongRacket,vecBallAlongRacket,Rhit,NCART);

}

/*
 * First order hold to interpolate linearly at time T
 * between racket pos,vel,normal entries
 *
 * IF T is nan, racket variables are assigned to zero-element of
 * relevant racket entries
 *
 */
static void interp_ball(const double** ballpred, const double T,
		                     const double dt, const int Nmax,
							 double *ballpos) {

	if (isnan(T)) {
		printf("Warning: T value is nan!\n");

		for(int i = 0; i < NCART; i++) {
			ballpos[i] = ballpred[i][0];
		}
	}
	else {
		int N = (int) (T/dt);
		double Tdiff = T - N*dt;

		for (int i = 0; i < NCART; i++) {
			if (N < Nmax) {
				ballpos[i] = ballpred[i][N] +
			(Tdiff/dt) * (ballpred[i][N+1] - ballpred[i][N]);
			}
			else {
				ballpos[i] = ballpred[i][N];
			}
		}
	}
}

/*
 * Punish the robot sufficiently as to induce proper landing behaviour
 *
 * The desired landing location chosen is the center of the table
 */
static double punish_land_robot(const double **ballpred,
		                        const double *x,
								const double *Rland,
								const double Rnet) {

	int i;
	static Vector  racketVel;
	static Vector  xfdot;
	static Vector  normal;
	static Vector  pos;
	static int firsttime = TRUE;
	static double netTime;
	static double landTime;
	static double g;
	static SL_Cstate ballPredict;
	static Vector ballVel;
	// desired landing locations
	static double xDesLand;
	static double yDesLand;
	static double zDesNet;

	double zNet;
	double xLand;
	double yLand;
	double T = x[2*NDOF];

	/* initialization of static variables */
	if (firsttime) {
		firsttime = FALSE;
		racketVel = my_vector(1,NCART);
		xfdot = my_vector(1,2*NCART);
		normal = my_vector(1,NCART);
		pos = my_vector(1,NCART);
		ballVel = my_vector(1,NCART);
		bzero((char *)&ballPredict, sizeof(ballPredict));
		xDesLand = 0.0;
		yDesLand = dist_to_table - 3*table_length/4;
		zDesNet = floor_level - table_height + net_height + 0.5;
	}
	// calculate deviation of ball to racket - hitting constraints
	calc_racket_state(x,pos,xfdot,normal);
	interp_ball(&ballPredict, T);
	for (i = 1; i <= NCART; i++) {
		racketVel[i] = xfdot[i];
		ballVel[i] = ballPredict.xd[i];
	}
	calc_times(x, &netTime, &landTime);
	racket_contact_model(racketVel, normal, ballVel);
	modify_ball_outgoing_vel(ballVel);
	zNet = ballPredict.x[Z] + netTime * ballVel[Z] + 0.5*gravity*netTime*netTime;
	xLand = ballPredict.x[X] + landTime * ballVel[X];
	yLand = ballPredict.x[Y] + landTime * ballVel[Y];

	return sqr(xLand - xDesLand)*Rland[0] +
			sqr(yLand - yDesLand)*Rland[1] +
			sqr(zNet - zDesNet)*Rnet[0];

}

/*
 * Punish the deviations of resting joints (as a function of optim params) from a chosen init_joint_state
 */
static double punish_wait_robot(double *x, double *Rwait, double landTime) {

	int i;
	static int firsttime = TRUE;
	static double qwait[NDOF];
	static double qdiff[NDOF];
	static double qrest[NDOF];

	if (firsttime) {
		firsttime = FALSE;
		for (i = 0; i < NDOF; i++) {
			qwait[i] = init_joint_state[i+1].th;
		}
	}

	for (i = 0; i < NDOF; i++) {
		qrest[i] = x[i] + x[i+NDOF] * landTime/2;
		qdiff[i] = qrest[i] - qwait[i];
	}

	return inner_w_prod(qdiff,qdiff,Rwait,NDOF);

}

/*
 * This is the constraint that makes sure we land the ball
 *
 */
static void land_ineq_constr(unsigned m, double *result, unsigned n, double *x, double *grad, static void *data) {

	int i;
	static Vector  racketVel;
	static Vector  xfdot;
	static Vector  normal;
	static Vector  pos;
	static int firsttime = TRUE;
	static double netTime;
	static double landTime;
	static double table_xmax;
	static double table_ymax;
	static double g;
	static double wall_z;
	static double net_y;
	static double net_z;
	static SL_Cstate ballPredict;
	static Vector vecBallAlongRacket;
	static Vector ballVel;

	double distB2RAlongN;
	double zNet;
	double xLand;
	double yLand;
	double T = x[2*NDOF];

	/* initialization of static variables */
	if (firsttime) {
		firsttime = FALSE;
		racketVel = my_vector(1,NCART);
		xfdot = my_vector(1,2*NCART);
		normal = my_vector(1,NCART);
		pos = my_vector(1,NCART);
		ballVel = my_vector(1,NCART);
		table_xmax = table_width/2.0;
		table_ymax = dist_to_table - table_length;
		wall_z = 1.0;
		net_y = dist_to_table - table_length/2;
		net_z = floor_level - table_height + net_height;
		vecBallAlongRacket = my_vector(1,NCART);
		bzero((char *)&ballPredict, sizeof(ballPredict));
	}

	interp_ball(&ballPredict, T);
	// calculate deviation of ball to racket - hitting constraints
	calc_racket_state(x,pos,xfdot,normal);

	for (i = 1; i <= NCART; i++) {
		vecBallAlongRacket[i] = ballPredict.x[i];
		ballVel[i] = ballPredict.xd[i];
		racketVel[i] = xfdot[i];
	}

	// calculate landing constraints
	racket_contact_model(racketVel, normal, ballVel);
	modify_ball_outgoing_vel(ballVel);
	calc_times(x, &netTime, &landTime);
	vec_sub(vecBallAlongRacket,pos,vecBallAlongRacket);
	distB2RAlongN = vec_mult_inner(normal,vecBallAlongRacket);
	vec_mult_scalar(normal,distB2RAlongN,normal);
	vec_sub(vecBallAlongRacket,normal,vecBallAlongRacket);

	zNet = ballPredict.x[Z] + netTime * ballVel[Z] + 0.5*gravity*netTime*netTime;
	xLand = ballPredict.x[X] + landTime * ballVel[X];
	yLand = ballPredict.x[Y] + landTime * ballVel[Y];

	result[0] = vec_mult_inner(vecBallAlongRacket,vecBallAlongRacket) - sqr(racket_radius);
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
static void hit_eq_constr(unsigned m, double *result, unsigned n, double *x, double *grad, static void *data) {

	int i;
	static Vector vecBallAlongRacket;
	static Vector  racketVel;
	static Vector  xfdot;
	static Vector  normal;
	static Vector  pos;
	static SL_Cstate ballPredict;
	static int firsttime = TRUE;
	double distB2RAlongN = 0;

	if (firsttime) {
		firsttime = FALSE;
		racketVel = my_vector(1,NCART);
		xfdot = my_vector(1,2*NCART);
		normal = my_vector(1,NCART);
		pos = my_vector(1,NCART);
		vecBallAlongRacket = my_vector(1,NCART);
		bzero((char *)&ballPredict, sizeof(ballPredict));
	}

	interp_ball(&ballPredict, x[2*NDOF]);
	for (i = 1; i <= NCART; i++) {
		vecBallAlongRacket[i] = ballPredict.x[i];
		racketVel[i] = xfdot[i];
	}

	// calculate deviation of ball to racket - hitting constraints
	calc_racket_state(x,pos,xfdot,normal);
	vec_sub(vecBallAlongRacket,pos,vecBallAlongRacket);
	distB2RAlongN = vec_mult_inner(normal,vecBallAlongRacket);
	vec_mult_scalar(normal,distB2RAlongN,normal);
	vec_sub(vecBallAlongRacket,normal,vecBallAlongRacket);

	result[0] = distB2RAlongN - ball_radius;
	//result[1] = vec_mult_inner(vecBallAlongRacket,vecBallAlongRacket);
}

/*
 * This is the inequality constraint that makes sure we never exceed the
 * joint limits during the striking and returning motion
 *
 */
static void joint_limits_ineq_constr(unsigned m, double *result,
		unsigned n, const double *x, double *grad, void *my_func_params) {

	static double a1[NDOF];
	static double a2[NDOF];
	static double a1ret[NDOF]; // coefficients for the returning polynomials
	static double a2ret[NDOF];
	static double qdot_rest[NDOF];
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
	static int firsttime = TRUE;

	if (firsttime) {
		firsttime = FALSE;
		coptim* optim_data = (coptim*)my_func_params;
		q0 = optim_data->q0;
		q0dot = optim_data->q0dot;
		qrest = optim_data->qrest;
		ub = optim_data->ub;
		lb = optim_data->lb;
		Tret = optim_data->time2return;
	}

	// calculate the polynomial coeffs which are used for checking joint limits
	calc_strike_poly_coeff(q0,q0dot,x,a1,a2);
	calc_return_poly_coeff(qrest,qdot_rest,x,Tret,a1ret,a2ret);
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
	double T = x[2*NDOF];

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
static double test_optim(double *x, double *params, int verbose) {

	// give info on constraint violation
	double *grad = FALSE;
	static double hit_violation[EQ_HIT_CONSTR_DIM];
	static double land_violation[INEQ_LAND_CONSTR_DIM];
	static double lim_violation[INEQ_JOINT_CONSTR_DIM]; // joint limit violations on strike and return
	joint_lim_ineq_constr(INEQ_JOINT_CONSTR_DIM, lim_violation, OPTIM_DIM, x, grad, NULL);
	hit_eq_constr(EQ_HIT_CONSTR_DIM, hit_violation, OPTIM_DIM, x, grad, NULL);
	land_ineq_constr(INEQ_LAND_CONSTR_DIM, land_violation, OPTIM_DIM, x, grad, NULL);
	double cost = costfunc(OPTIM_DIM, x, grad, params);

	if (verbose) {
		// give info on solution vector
		print_optim_vec(x);
		printf("Hitting violations:\n");
		printf("Ball to racket normal distance: %.2f\n",hit_violation[0]);
		printf("Ball distance projected to racket: %.2f\n", land_violation[0]);
		//printf("Ball along racket: %.2f\n", hit_violation[1]);
		printf("Landing info:\n");
		printf("netTime: %f\n", -land_violation[1]);
	    printf("landTime - netTime: %f\n", -land_violation[4]);
		/*printf("wall_z - zNet: %f\n", -land_violation[2]);
		printf("zNet - net_z: %f\n", -land_violation[3]);
		printf("table_xmax - xLand: %f\n", -land_violation[5]);
		printf("table_xmax + xLand: %f\n", -land_violation[6]);
		printf("net_y - yLand: %f\n", -land_violation[7]);
	    printf("yLand - table_ymax: %f\n", -land_violation[8]);*/
	}

	for (int i = 0; i < EQ_HIT_CONSTR_DIM; i++) {
		if (hit_violation[i] > 0.0) {
			printf("Hitting equality constraints violated!\n");
			break;
		}
	}
	for (int i = 0; i < INEQ_LAND_CONSTR_DIM; i++) {
		if (land_violation[i] > 0.0) {
			printf("Landing inequality constraints violated!\n");
			break;
		}
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
static void modify_ball_outgoing_vel(Vector ballVel) {

	static double multiplier_x, multiplier_y, multiplier_z;
	static int firsttime = TRUE;

	if (firsttime) {
		firsttime = FALSE;
		multiplier_x = 0.9;
		multiplier_y = 0.8;
		multiplier_z = 0.83;
	}

	ballVel[1] = ballVel[1] * multiplier_x;
	ballVel[2] = ballVel[2] * multiplier_y;
	ballVel[3] = ballVel[3] * multiplier_z;

}

/*
 * Calculate the net hitting time and the landing time
 * Assuming that there was an interaction (hitting) between racket and ball
 *
 */
static void calc_times(const double *x, double *netTime, double *landTime) {

	int i;
	static Vector  xfdot;
	static Vector  normal;
	static Vector  pos;
	static Vector  racketVel;
	static int     firsttime = TRUE;
	static SL_Cstate ballPredict;
	static Vector ballVel;
	static double table_z;
	static double g;
	static double net_y;
	double distBall2TableZ;

	/* initialization of static variables */
	if (firsttime) {
		firsttime = FALSE;
		table_z = floor_level - table_height + ball_radius;
		net_y = dist_to_table - table_length/2;
		g = fabs(gravity);
		racketVel = my_vector(1,NCART);
		xfdot = my_vector(1,2*NCART);
		normal = my_vector(1,NCART);
		pos = my_vector(1,NCART);
		ballVel = my_vector(1,NCART);
		bzero((char *)&ballPredict, sizeof(ballPredict));
	}
	calc_racket_state(x,pos, xfdot, normal);
	interp_ball(&ballPredict, x[2*NDOF]);
	for (i = 1; i <= NCART; i++) {
		racketVel[i] = xfdot[i];
		ballVel[i] = ballPredict.xd[i];
	}
	racket_contact_model(racketVel, normal, ballVel);
	modify_ball_outgoing_vel(ballVel);
	distBall2TableZ = ballPredict.x[Z] - table_z;

	if (sqr(ballVel[3]) > -2*g*distBall2TableZ) {
		*landTime = fmax(ballVel[3] + sqrt(sqr(ballVel[3]) + 2*g*distBall2TableZ),
		                ballVel[3] - sqrt(sqr(ballVel[3]) + 2*g*distBall2TableZ)) / g;
	}
	else {
		// landTime is not real!
		*landTime = 1.0;
	}

	*netTime = (net_y - ballPredict.x[Y])/ballVel[2];
}

/*
 * Finalize the solution and update target SL structure and hitTime value
 */
static static void finalize_soln(const double* x, optim * params, double time_elapsed) {

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
 */
static void set_penalty_matrices(weights * pen) {

	pen = (weights *)malloc(sizeof(weights));
	double* R1 = (double*)calloc(NDOF,sizeof(double));
	double* R2 = (double*)calloc(NDOF,sizeof(double));
	double* Rwait = (double*)calloc(NDOF,sizeof(double));
	double* Rhit = (double*)calloc(NDOF,sizeof(double));
	double* Rland = (double*)calloc(2,sizeof(double));
	double Rnet;

	const_vec(NDOF,10.0,R1);
	const_vec(NDOF,1.0,R2);
	const_vec(NDOF,1.0,Rwait);
	const_vec(NCART,10.0,Rhit);
	const_vec(1,1.0,Rnet);
	const_vec(2,1.0,Rland);

	pen->R_hit = Rhit;
	pen->R_land = Rland;
	pen->R_return = R2;
	pen->R_strike = R1;
	pen->R_wait = Rwait;
	pen->R_net = Rnet;
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
