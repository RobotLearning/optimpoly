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

/*
 * Multi-threading entry point for the NLOPT optimization
 */
double nlopt_optim_lazy_run(coptim *coparams,
	                     racketdes *racketdata,
		                 optim *params) {

	static int count = 0;
	printf("==========================================\n");
	printf("Running NLOPT\n");
	//initTime = get_time();
	nlopt_optim_fixed_run(coparams,racketdata,params);
	double maxviol = nlopt_optim_lazy(coparams,params);
	printf("Optim count: %d\n", (++count));
	printf("==========================================\n");
	return maxviol;
}

/*
 * NLOPT optimization routine for table tennis LAZY PLAYER (LP)
 *
 * If constraints are violated, it will not modify the lookup values (input)
 */
double nlopt_optim_lazy(coptim *coparams, optim *params) {

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
	nlopt_set_lower_bounds(opt, coparams->lb);
	nlopt_set_upper_bounds(opt, coparams->ub);
	nlopt_set_min_objective(opt, costfunc, coparams);
	nlopt_add_equality_mconstraint(opt, EQ_HIT_CONSTR_DIM, hit_eq_constr, NULL, tol_eq_hit);
	nlopt_add_inequality_mconstraint(opt, INEQ_LAND_CONSTR_DIM,
			land_ineq_constr, coparams, tol_ineq_land);
	nlopt_add_inequality_mconstraint(opt, INEQ_JOINT_CONSTR_DIM,
			joint_lim_ineq_constr, coparams, tol_ineq_joint);
	nlopt_set_xtol_rel(opt, 1e-2);
	set_penalty_matrices();

	double init_time = get_time();
	double past_time = 0.0;
	double minf; // the minimum objective value, upon return //
	int res; // error code
	double max_violation;

	if ((res = nlopt_optimize(opt, x, &minf)) < 0) {
	    printf("NLOPT failed with exit code %d!\n", res);
	    past_time = (get_time() - init_time)/1e3;
	    max_violation = test_optim(x,coparams,racketdata,TRUE);
	}
	else {
		past_time = (get_time() - init_time)/1e3;
		printf("NLOPT success with exit code %d!\n", res);
		printf("NLOPT took %f ms\n", past_time);
	    printf("Found minimum at f = %0.10g\n", minf);
	    max_violation = test_optim(x,coparams,racketdata,TRUE);
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
double costfunc(unsigned n, double *x, double *grad, void *my_func_params) {

	int i;
	static double J1;
	static double J2;
	static double Jhit, Jland, Jwait;
	static double a1[NDOF];
	static double a2[NDOF];
	static double b1[NDOF];
	static double b2[NDOF];
	double T = x[2*NDOF];
	double T2 = 0.0;
	double netTime = 0;
	double landTime = 0;

	if (grad) {
		//TODO:
	}

	// calculate the polynomial coeffs which are used in the cost calculation
	calc_strike_poly_coeff_lp(a1,a2,x);

	// calculate the landing time
	calc_times(x, &netTime, &landTime);

	// calculate the follow up coeffs which are added to the cost
	calc_return_poly_coeff_lp(b1,b2,x,landTime);

	J1 = T * (3*T*T*inner_inv_weighted_prod_lp(a1,a1,R1,NDOF) +
			3*T*inner_inv_weighted_prod_lp(a1,a2,R1,NDOF) + inner_inv_weighted_prod_lp(a2,a2,R1,NDOF));

	T2 = landTime;
	J2 = T2 * (3*T2*T2*inner_inv_weighted_prod_lp(b1,b1,R2,NDOF) +
			3*T2*inner_inv_weighted_prod_lp(b1,b2,R2,NDOF) + inner_inv_weighted_prod_lp(b2,b2,R2,NDOF));

	Jhit = punish_hit_robot(x);
	Jland = punish_land_robot(x);
	Jwait = punish_wait_robot(x,landTime);

	return J1 + Jhit + Jland + Jwait; //J2
}

/*
 * Punish the robot to encourage hitting
 *
 * We punish the square of the deviations of the ball projected to the racket to the racket center
 */
double punish_hit_robot(double *x) {

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

	return inner_weighted_prod_lp(vecBallAlongRacket,vecBallAlongRacket,Rhit,NCART);

}

/*
 * Punish the robot sufficiently as to induce proper landing behaviour
 *
 * The desired landing location chosen is the center of the table
 */
double punish_land_robot(double *x) {

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
		Rnet[0] = 1.0;
		Rland[0] = 1.0;
		Rland[1] = 1.0;
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
	zNet = ballPredict.x[_Z_] + netTime * ballVel[_Z_] + 0.5*gravity*netTime*netTime;
	xLand = ballPredict.x[_X_] + landTime * ballVel[_X_];
	yLand = ballPredict.x[_Y_] + landTime * ballVel[_Y_];

	return sqr(xLand - xDesLand)*Rland[0] + sqr(yLand - yDesLand)*Rland[1] + sqr(zNet - zDesNet)*Rnet[0];

}

/*
 * Punish the deviations of resting joints (as a function of optim params) from a chosen init_joint_state
 */
double punish_wait_robot(double *x, double landTime) {

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

	return inner_weighted_prod(qdiff,qdiff,Rwait,NDOF);

}

/*
 * This is the constraint that makes sure we land the ball
 *
 */
void land_ineq_constr(unsigned m, double *result, unsigned n, double *x, double *grad, void *data) {

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
void hit_eq_constr(unsigned m, double *result, unsigned n, double *x, double *grad, void *data) {

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
void joint_lim_ineq_constr(unsigned m, double *result, unsigned n, double *x, double *grad, void *my_func_params) {

	int i;
	static double a1[NDOF];
	static double a2[NDOF];
	static double b1[NDOF];
	static double b2[NDOF];
	static double joint_strike_max_cand[NDOF];
	static double joint_strike_min_cand[NDOF];
	static double joint_return_max_cand[NDOF];
	static double joint_return_min_cand[NDOF];
	double netTime, landTime = 0;
	calc_times(x, &netTime, &landTime);

	// calculate the polynomial coeffs which are used for checking joint limits
	calc_strike_poly_coeff_lp(a1,a2,x);
	calc_return_poly_coeff_lp(b1,b2,x,landTime);
	// calculate the candidate extrema both for strike and return
	calc_strike_ext_cand_lp(a1,a2,x[2*NDOF],joint_strike_max_cand,joint_strike_min_cand);
	calc_return_ext_cand_lp(b1,b2,x,landTime,joint_return_max_cand, joint_return_min_cand);

	/* deviations from joint min and max */
	for (i = 0; i < NDOF; i++) {
		result[i] = joint_strike_max_cand[i] - joint_range[i+1][MAX_THETA];
		result[i+NDOF] = joint_range[i+1][MIN_THETA] - joint_strike_min_cand[i];
		result[i+2*NDOF] = joint_return_max_cand[i] - joint_range[i+1][MAX_THETA];
		result[i+3*NDOF] = joint_range[i+1][MIN_THETA] - joint_return_min_cand[i];
		//printf("%f %f %f %f\n", result[i],result[i+DOF],result[i+2*DOF],result[i+3*DOF]);
	}

}

/*
 * Calculate the polynomial coefficients from the optimized variables qf,qfdot,T
 * p(t) = a1*t^3 + a2*t^2 + a3*t + a4
 */
void calc_strike_poly_coeff(double *a1, double *a2, double *x) {

	// variables to be optimized
	int i;
	double T = x[2*NDOF];

	for (i = 0; i < NDOF; i++) {
		a1[i] = (2/pow(T,3))*(current_joint_state[i+1].th-x[i]) + (1/(T*T))*(current_joint_state[i+1].thd + x[i+NDOF]);
		a2[i] = (3/(T*T))*(x[i]-current_joint_state[i+1].th) - (1/T)*(x[i+NDOF] + 2*current_joint_state[i+1].thd);
	}
}

/*
 * Calculate the returning polynomial coefficients from the optimized variables qf,qfdot
 * and time2return variable defined in the header
 * p(t) = a1*t^3 + a2*t^2 + a3*t + a4
 *
 * Assuming that R_wait is diagonal here
 */
void calc_return_poly_coeff(double *a1, double *a2, double *x, double tau) {

	double weights[NDOF];
	int i;
	for (i = 0; i < NDOF; i++) {
		weights[i] = Rwait[i] / (6 + pow(tau,3)*Rwait[i]);
		a1[i] = weights[i] * (x[i+NDOF]*tau + 2*x[i] - 2*init_joint_state[i+1].th);
		a2[i] = -x[i+NDOF]/(2*tau) - (3/2)*tau*a1[i];
	}

	// variables to be optimized
	/*static double q0dot[NDOF];
	static double T = 1.0;
	int i;

	for (i = 0; i < NDOF; i++) {
		a1[i] = (2/pow(T,3))*(x[i]-init_joint_state[i+1].th) + (1/(T*T))*(init_joint_state[i+1].thd + x[i+NDOF]);
		a2[i] = (3/(T*T))*(init_joint_state[i+1].th-x[i]) - (1/T)*(2*x[i+NDOF] + init_joint_state[i+1].thd);
	}*/
}

/*
 * Calculate the extrema candidates for each joint (2*7 candidates in total)
 * For the striking polynomial
 * Clamp to [0,T]
 *
 */
void calc_strike_ext_cand(double *a1, double *a2, double T, double *joint_max_cand, double *joint_min_cand) {

	static double q0dot[NDOF];
	int i;
	static double cand1, cand2;

	for (i = 0; i < NDOF; i++) {
		// find the extrema time candidates
		cand1 = fmin(T,fmax(0,(-a2[i] + sqrt(a2[i]*a2[i] - 3*a1[i]*current_joint_state[i+1].thd))/(3*a1[i])));
		cand2 =  fmin(T,fmax(0,(-a2[i] - sqrt(a2[i]*a2[i] - 3*a1[i]*current_joint_state[i+1].thd))/(3*a1[i])));
		// find the joint extrema values at those candidate points
		cand1 = a1[i]*pow(cand1,3) + a2[i]*pow(cand1,2) + current_joint_state[i+1].thd*cand1 + current_joint_state[i+1].th;
		cand2 = a1[i]*pow(cand2,3) + a2[i]*pow(cand2,2) + current_joint_state[i+1].thd*cand2 + current_joint_state[i+1].th;
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
void calc_return_ext_cand(double *a1, double *a2, const double *x, double landTime,
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
int test_optim(double *x, double *params, int verbose) {

	// give info on constraint violation
	double *grad = FALSE;
	static double hit_violation[EQ_HIT_CONSTR_DIM];
	static double land_violation[INEQ_LAND_CONSTR_DIM];
	static double lim_violation[INEQ_JOINT_CONSTR_DIM]; // joint limit violations on strike and return
	joint_lim_ineq_constr_lp(INEQ_JOINT_CONSTR_DIM, lim_violation, OPTIM_DIM, x, grad, NULL);
	hit_eq_constr_lp(EQ_HIT_CONSTR_DIM, hit_violation, OPTIM_DIM, x, grad, NULL);
	land_ineq_constr_lp(INEQ_LAND_CONSTR_DIM, land_violation, OPTIM_DIM, x, grad, NULL);
	double cost = costfunc_lp(OPTIM_DIM, x, grad, params);

	if (verbose) {
		// give info on solution vector
		print_optim_vec_lp(x);
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

	int i;

	for (i = 0; i < EQ_HIT_CONSTR_DIM; i++) {
		if (hit_violation[i] > 10*SLACK) {
			printf("Hitting equality constraints violated!\n");
			printf("Not updating trajectory!\n");
			return FALSE;
		}
	}
	for (i = 0; i < INEQ_LAND_CONSTR_DIM; i++) {
		if (land_violation[i] > 10*SLACK) {
			printf("Landing inequality constraints violated!\n");
			printf("Not updating trajectory!\n");
			return FALSE;
		}
	}
	for (i = 0; i < INEQ_JOINT_CONSTR_DIM; i++) {
		if (lim_violation[i] > SLACK) {
			printf("Joint limit violated by %.2f on joint %d\n", lim_violation[i], i % NDOF + 1);
			printf("Not updating trajectory!\n");
			return FALSE;
		}
	}

	return TRUE;
}

/*
 * First order hold to interpolate linearly at time T between ball prediction matrix ballMat entries
 *
 */
void interp_ball(SL_Cstate *ballPredict, double T) {

	int i;
	int N = (int) (T/TSTEP);
	double Tdiff = T - N*TSTEP;
	static int iter;
	int Nmax = (int) TPRED/TSTEP;

	if (isnan(T)) {
		printf("Warning: T value is nan! "
				"Setting ball to first index.\n");
		for (i = 1; i <= NCART; i++) {
			ballPredict->x[i] = ballPath[0].x[i];
			ballPredict->xd[i] = ballPath[0].xd[i];
		}
		return;
	}
	if (T > TPRED) {
		printf("Warning: Extrapolation! T value is greater than TPRED = %f. "
				"Setting ball to last index \n", TPRED);
		for (i = 1; i <= NCART; i++) {
			ballPredict->x[i] = ballPath[Nmax].x[i];
			ballPredict->xd[i] = ballPath[Nmax].xd[i];
		}
		return;
	}

	for (i = 1; i <= NCART; i++) {
		if (N < Nmax) {
			ballPredict->x[i] = ballPath[N].x[i] + (Tdiff/TSTEP) * (ballPath[N+1].x[i] - ballPath[N].x[i]);
			ballPredict->xd[i] = ballPath[N].xd[i] + (Tdiff/TSTEP) * (ballPath[N+1].xd[i] - ballPath[N].xd[i]);
		}
		else {
			ballPredict->x[i] = ballPath[N].x[i];
			ballPredict->xd[i] = ballPath[N].xd[i];
		}
	}

}

/*
 * Modify ball outgoing velocity by multiplying with a constant vector (less than one)
 *
 * Necessary since the landing and net hitting locations are calculated using projectile motion
 *
 */
void modify_ball_outgoing_vel(Vector ballVel) {

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
void calc_times(double *x, double *netTime, double *landTime) {

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
		g = intern_gravity;
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
	distBall2TableZ = ballPredict.x[_Z_] - table_z;

	if (sqr(ballVel[3]) > -2*g*distBall2TableZ) {
		*landTime = fmax(ballVel[3] + sqrt(sqr(ballVel[3]) + 2*g*distBall2TableZ),
		                ballVel[3] - sqrt(sqr(ballVel[3]) + 2*g*distBall2TableZ)) / g;
	}
	else {
		// landTime is not real!
		*landTime = 1.0;
	}

	*netTime = (net_y - ballPredict.x[_Y_])/ballVel[2];
}

/*
 * Set the initial solution for NLOPT 2*dof + 1 dimensional problem
 * to initial posture with zero velocity
 *
 * Guess the final time to be 0.5;
 *
 */
void init_soln_to_rest_posture(double *x) {

	// initialize first dof entries to q0
	int i;
	for (i = 0; i < NDOF; i++) {
		x[i] = init_joint_state[i+1].th;
		x[i+NDOF] = 0.0;
	}
	x[2*NDOF] = 0.5;
}

/*
 * target and hitTime must be initialized in the lookup table
 * can be used for MPC
 */
void init_soln_to_last(double *x, SL_DJstate target[], double *hitTime) {

	// initialize first dof entries to q0
	int i;
	for (i = 0; i < NDOF; i++) {
		x[i] = target[i+1].th;
		x[i+NDOF] = target[i+1].thd;
	}
	x[2*NDOF] = *hitTime;
}

/*
 * Finalize the solution and update target SL structure and hitTime value
 */
void finalize_soln_cp(double *x, SL_DJstate target[], double *hitTime) {

	int i;
	for (i = 0; i < NDOF; i++) {
		target[i+1].th = x[i];
		target[i+1].thd = x[i+NDOF];
	}
	*hitTime = x[2*NDOF];
	//print_optim_vec(x);
}

/*
 * Finalize the solution and update target SL structure and hitTime value
 */
void finalize_soln_lp(double *x, SL_DJstate target[], double *hitTime) {

	int i;
	for (i = 0; i < NDOF; i++) {
		target[i+1].th = x[i];
		target[i+1].thd = x[i+NDOF];
	}
	double optim_time = (get_time() - initTime) / 1e3;
	printf("NLOPT took %f ms\n", optim_time);
	*hitTime = x[2*NDOF] - optim_time/1e3;
	//print_optim_vec(x);
}

/*
 * Set the elements of the double vector to 1
 */
void ones_vector(double *vec, int size) {

	int i;
	for (i = 0; i < size; i++) {
		vec[i] = 1.0;
	}

}

/*
 * Initialize the penalty matrices used to punish the robot
 */
void set_penalty_matrices() {

	int i;
	ones_vector(R1,NDOF);
	for (i = 0; i < NDOF; i++) {
		R1[i] = 10 * R1[i];
	}
	ones_vector(R2,NDOF);
	ones_vector(Rwait,NDOF);
	ones_vector(Rhit,NCART);
	for (i = 0; i < NCART; i++) {
		Rhit[i] = 10 * Rhit[i];
	}
	Rnet[0] = 1.0;
	Rland[0] = 1.0;
	Rland[1] = 1.0;
}

/*
 * Returns the inner product between two vectors of size given in last argument
 */
double inner_prod_lp(double *a1, double *a2, int size) {

	int i;
	double val = 0.0;
	for (i = 0; i < size; i++) {
		val += a1[i]*a2[i];
	}
	return val;
}

/*
 * Returns the weighted inner product between two vectors of size given in last argument
 */
double inner_weighted_prod_lp(double *a1, double *a2, double *w, int size) {

	int i;
	double val = 0.0;
	for (i = 0; i < size; i++) {
		val += a1[i]*w[i]*a2[i];
	}
	return val;
}

/*
 * Returns the inverse weighted inner product between two vectors of size given in last argument
 */
double inner_inv_weighted_prod_lp(double *a1, double *a2, double *w, int size) {

	int i;
	double val = 0.0;
	for (i = 0; i < size; i++) {
		val += a1[i]*a2[i]/w[i];
	}
	return val;
}

/*
 * Returns a1 - a2 vector into a1, assuming both have dof = 7 length
 */
void vec_minus_lp(double *a1, double *a2) {

	int i;
	for (i = 0; i < NDOF; i++) {
		a1[i] = a1[i] - a2[i];
	}
}

/*
 * Set array entries all to the same value
 */
void const_vec_lp(const double val, const int n, double* array) {

	int i;
	for (i = 0; i < n; i++) {
		array[i] = val;
	}
}



