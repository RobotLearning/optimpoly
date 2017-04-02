/*
 * optim.c
 *
 * NLOPT polynomial optimization functions are stored here.
 *
 * Table tennis tasks make use of them to generate trajectories.
 *
 *  Created on: Aug 9, 2016
 *      Author: okan
 */

#include "stdio.h"
#include "string.h"
#include "sys/time.h"

#include "../../../player/old/ccode/kinematics.h"
#include "../../../player/old/ccode/optims.h"
#include "../../../player/old/ccode/SL0.h"
#include "../../../player/old/ccode/utils.h"

/*
 * Global variables used in lp_task
 * TODO: replace with semaphore/mutex?
 */
SL_Cstate *ballPath; // predicted ball pos and vel values for T_pred time seconds
Racket *racketDes; // racket strategy
SL_DJstate current_joint_state[DOF+1]; // for correction with MPC
int busy = FALSE; // to check between threads whether optimization is running

/*
 * Static global vars confined to optim code
 */
static double initTimeOpt; // time the optimization start at
static int count = 0; // number of optimizations that ran so far

/*
 * Multi-thread synchronization function
 *
 */
int mthread_sync_optim(nlopt_thread_data data, SL_DJstate target[], double *hitTime) {

	int i;
	if (check_thread_termination(data, target, *hitTime)) { // this means thread updated the optim variables
		*hitTime = data.hitTime;
		for (i = 1; i <= DOF; i++) {
			target[i].th = data.target[i].th;
			target[i].thd = data.target[i].thd;
		}
		trj_time = DT;
		return TRUE;
	}
	return FALSE;
}

/*
 * Check thread's successful termination.
 * We know that the thread finalizes the solution [by calling finalize_soln] before being detached.
 * The last operation is on the hitting time.
 * Hence we only need to check whether the hitting time was updated or not.
 *
 */
int check_thread_termination(nlopt_thread_data data, SL_DJstate target[], double hitTime) {

	int i;
	if (data.hitTime != hitTime) {
		//printf("hittimes: %f - %f\n", data.hitTime, hitTime);
		return TRUE;
	}
	return FALSE;
}

/*
 * Multi-threading entry point for the NLOPT optimization
 */
void* mthread_nlopt_optim(void *arg) {

	busy = TRUE;
	nlopt_thread_data *data = (nlopt_thread_data *) arg;
	printf("==========================================\n");
	printf("Running NLOPT\n");
	initTimeOpt = get_time();
	double hitTimeThreadSafe = data->hitTime;
	nlopt_optim_fixed_run(data->target,&hitTimeThreadSafe);
	printf("Optim count: %d\n", (++count));
	printf("==========================================\n");
	data->hitTime = hitTimeThreadSafe;
	busy = FALSE;
	pthread_exit(NULL);

	return NULL;
}

/*
 * NLOPT optimization routine for table tennis centred player (CP)
 *
 * If constraints are violated, it will not modify the lookup values (input)
 */
void nlopt_optim_fixed_run(SL_DJstate target[], double* hitTime) {

	static double tol_eq[EQ_CONSTR_DIM] = {1e-4};
	static double tol_ineq[INEQ_CONSTR_DIM] = {1e-4};
	static double lb[OPTIM_DIM]; /* lower bounds */
	static double ub[OPTIM_DIM]; /* upper bounds */
	static double x[OPTIM_DIM]; /* initial guess */
	static int firsttime = TRUE;

	static nlopt_opt opt;

	if (firsttime) {
		firsttime = FALSE;
		set_bounds(lb,ub,10*SLACK);
		opt = nlopt_create(NLOPT_LN_COBYLA, OPTIM_DIM); /* LN = does not require gradients */
		//nlopt_set_xtol_rel(opt, 1e-2);
		//opt = nlopt_create(NLOPT_AUGLAG, OPTIM_DIM); /* algorithm and dimensionality */
		//nlopt_set_local_optimizer(opt, opt);
		nlopt_set_lower_bounds(opt, lb);
		nlopt_set_upper_bounds(opt, ub);
		//nlopt_set_min_objective(opt, const_costfunc, NULL);
		nlopt_set_min_objective(opt, costfunc, NULL);
		nlopt_add_inequality_mconstraint(opt, INEQ_CONSTR_DIM, joint_lim_ineq_constr, NULL, tol_ineq);
		nlopt_add_equality_mconstraint(opt, EQ_CONSTR_DIM, kinematics_eq_constr, NULL, tol_eq);
		nlopt_set_xtol_rel(opt, 1e-2);
	}

	//int maxeval = 20000;
	//nlopt_set_maxeval(opt, maxeval);
	//double maxtime = 10e-3;
	//nlopt_set_maxtime(opt, maxtime);

	init_soln_to_rest_posture(x); //parameters are the initial joint positions q0
	//init_soln_to_last(x, target, hitTime);
	double minf; /* the minimum objective value, upon return */
	if (nlopt_optimize(opt, x, &minf) < 0) {
	    printf("NLOPT failed!\n");
	    //freeze();
	}
	else {
	    //printf("Found minimum at f = %0.10g\n", minf);
	    if (test_optim(x,NULL,TRUE))
	    	finalize_soln(x, target, hitTime);
	}
	//nlopt_destroy(opt);
}

/*
 * Constant cost function.
 * Used to project lookup table solutions to equality constraints.
 *
 * Note: seems to be much better than costfunc when the optim is restricted to 1ms.
 */
double const_costfunc(unsigned n, const double *x, double *grad, void *my_func_params) {

	int i;
	if (grad) {
		for (i = 0; i < DOF; i++) {
			grad[i] = 0.0;
			grad[i+DOF] = 0.0;
		}
		grad[2*DOF] = 0.0;
	}

	return 1;
}

/*
 * Calculates the cost function for table tennis centred player
 * to find spline (3rd order strike+return) polynomials
 */
double costfunc(unsigned n, const double *x, double *grad, void *my_func_params) {

	int i;
	double a1[DOF];
	double a2[DOF];
	static double q0vec[DOF]; //double representation of q0 th structure
	static double q0dot[DOF]; // all zeros
	static double qfdot[DOF]; // all zeros
	static double qf[DOF]; // opt value
	double T = x[2*DOF];

	if (grad) {

		for (i = 0; i < DOF; i++) {
			qf[i] = x[i];
			qfdot[i] = x[i+DOF];
			q0dot[i] = current_joint_state[i].thd;
			q0vec[i] = current_joint_state[i].th;
			grad[i] = (6/pow(T,3))*(qf[i] - q0vec[i]) - (3/(T*T))*(q0dot[i] + qfdot[i]);
			grad[i+DOF] = (-3/(T*T))*(qf[i] - q0vec[i]) + (1/T)*(2*qfdot[i] + q0dot[i]);
		}
		//time derivative of cost J
		vec_minus(qf,q0vec);
		//assuming q0dot = 0 here;
		grad[2*DOF] = (-9/pow(T,4))*inner_prod(qf,qf) +
					  (6/pow(T,3))*inner_prod(qf,qfdot) +
					  (3/(T*T))*inner_prod(q0dot,qfdot) -
					  (1/(T*T))*inner_prod(qfdot,qfdot);
	}



	// calculate the polynomial coeffs which are used in the cost calculation
	calc_strike_poly_coeff(a1,a2,x);

	return T * (3*T*T*inner_prod(a1,a1) +
			3*T*inner_prod(a1,a2) + inner_prod(a2,a2));
}

/*
 * Returns a1 - a2 vector into a1, assuming both have dof = 7 length
 */
void vec_minus(double *a1, const double *a2) {

	int i;
	for (i = 0; i < DOF; i++) {
		a1[i] = a1[i] - a2[i];
	}
}

/*
 * This is the constraint that makes sure we hit the ball
 *
 */
void kinematics_eq_constr(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data) {

	static double ballPred[CART];
	static double racketDesVel[CART];
	static double racketDesNormal[CART];
	static SL_DJstate joint_local_state[DOF+1];
	static Matrix  racketTransform;
	static Matrix  Jacobi;
	static Vector  qfdot;
	static Vector  xfdot;
	static Vector  normal;
	static int     firsttime = TRUE;

	double T = x[2*DOF];
	int N = T/TSTEP;
	int i;

	/* initialization of static variables */
	if (firsttime) {
		firsttime = FALSE;
		Jacobi               = my_matrix(1,2*CART,1,DOF);
		racketTransform      = my_matrix(1,4,1,4);
		qfdot                = my_vector(1,DOF);
		xfdot                = my_vector(1,2*CART);
		normal               = my_vector(1,CART);
		bzero((char *)&(joint_local_state[1]), DOF * sizeof(joint_local_state[1]));

		// homogeneous transform instead of using quaternions
		racketTransform[1][1] = 1;
		racketTransform[2][3] = 1;
		racketTransform[3][2] = -1;
		racketTransform[4][4] = 1;
	}

	if (grad) {
		// compute gradient of kinematics = jacobian
		//TODO:
		grad[0] = 0.0; grad[1] = 0.0; grad[2] = 0.0;
	}

	// interpolate at time T to get the desired racket parameters
	first_order_hold(ballPred,racketDesVel,racketDesNormal,&T);

	// extract state information from array to joint_des_state structure
	for (i = 1; i <= DOF; i++) {
		joint_local_state[i].th = x[i-1];
		joint_local_state[i].thd = qfdot[i] = x[i-1+DOF];
	}

	/* compute the desired link positions */
	kinematics(joint_local_state, &base_state, &base_orient, endeff,
			           joint_cog_mpos_des, joint_axis_pos_des, joint_origin_pos_des,
			           link_pos_des, Alink_des);

	/* compute the racket normal */
	mat_mult(Alink_des[6],racketTransform,Alink_des[6]);
	for (i = 1; i <= CART; i++) {
		normal[i] = Alink_des[6][i][3];
	}

	/* compute the jacobian */
	jacobian(link_pos_des, joint_origin_pos_des, joint_axis_pos_des, Jacobi);
	mat_vec_mult(Jacobi, qfdot, xfdot);

	/* deviations from the desired racket frame */
	for (i = 1; i <= CART; i++) {
		//printf("xfdot[%d] = %.4f, racketDesVel[%d] = %.4f\n",i,xfdot[i],i,racketDesVel[i-1]);
		//printf("normal[%d] = %.4f, racketDesNormal[%d] = %.4f\n",i,normal[i],i,racketDesNormal[i-1]);
		result[i-1] = link_pos_des[6][i] - ballPred[i-1];
		result[i-1 + CART] = xfdot[i] - racketDesVel[i-1];
		result[i-1 + 2*CART] = normal[i] - racketDesNormal[i-1];
	}

}

/*
 * Set upper and lower bounds on the optimization
 */
void set_bounds(double *lb, double *ub, double slack) {

	// lower bounds and upper bounds for qf are the joint limits
	int i;
	for (i = 1; i <= DOF; i++) {
		ub[i-1] = joint_range[i][MAX_THETA] - slack;
		lb[i-1] = joint_range[i][MIN_THETA] + slack;
		ub[i-1+DOF] = MAX_VEL;
		lb[i-1+DOF] = -MAX_VEL;
	}
	// constraints on final time
	ub[2*DOF] = TPRED;
	lb[2*DOF] = 0.00;
}

/*
 * This is the inequality constraint that makes sure we never exceed the
 * joint limits during the striking and returning motion
 *
 */
void joint_lim_ineq_constr(unsigned m, double *result, unsigned n, const double *x, double *grad, void *my_func_params) {

	int i;
	static double a1[DOF];
	static double a2[DOF];
	static double a1ret[DOF]; // coefficients for the returning polynomials
	static double a2ret[DOF];
	static double joint_strike_max_cand[DOF];
	static double joint_strike_min_cand[DOF];
	static double joint_return_max_cand[DOF];
	static double joint_return_min_cand[DOF];

	// calculate the polynomial coeffs which are used for checking joint limits
	calc_strike_poly_coeff(a1,a2,x);
	calc_return_poly_coeff(a1ret,a2ret,x);
	// calculate the candidate extrema both for strike and return
	calc_strike_ext_cand(a1,a2,x[2*DOF],joint_strike_max_cand,joint_strike_min_cand);
	calc_return_ext_cand(a1ret,a2ret,x,joint_return_max_cand,joint_return_min_cand);

	/* deviations from joint min and max */
	for (i = 0; i < DOF; i++) {
		result[i] = joint_strike_max_cand[i] - joint_range[i+1][MAX_THETA];
		result[i+DOF] = joint_range[i+1][MIN_THETA] - joint_strike_min_cand[i];
		result[i+2*DOF] = joint_return_max_cand[i] - joint_range[i+1][MAX_THETA];
		result[i+3*DOF] = joint_range[i+1][MIN_THETA] - joint_return_min_cand[i];
		//printf("%f %f %f %f\n", result[i],result[i+DOF],result[i+2*DOF],result[i+3*DOF]);
	}

}

/*
 * Returns the inner product between two vectors of size DOF
 */
double inner_prod(const double *a1, const double *a2) {

	int i;
	double val = 0.0;
	for (i = 0; i < DOF; i++) {
		val += a1[i]*a2[i];
	}
	return val;
}

/*
 * Calculate the polynomial coefficients from the optimized variables qf,qfdot,T
 * p(t) = a1*t^3 + a2*t^2 + a3*t + a4
 */
void calc_strike_poly_coeff(double *a1, double *a2, const double *x) {

	// variables to be optimized
	int i;
	double T = x[2*DOF];

	for (i = 0; i < DOF; i++) {
		a1[i] = (2/pow(T,3))*(current_joint_state[i+1].th-x[i]) + (1/(T*T))*(current_joint_state[i+1].thd + x[i+DOF]);
		a2[i] = (3/(T*T))*(x[i]-current_joint_state[i+1].th) - (1/T)*(x[i+DOF] + 2*current_joint_state[i+1].thd);
	}
	return;
}

/*
 * Calculate the returning polynomial coefficients from the optimized variables qf,qfdot
 * and time2return variable defined in the header
 * p(t) = a1*t^3 + a2*t^2 + a3*t + a4
 */
void calc_return_poly_coeff(double *a1, double *a2, const double *x) {

	// variables to be optimized
	static double q0dot[DOF];
	static double T = TIME2RETURN;
	int i;

	for (i = 0; i < DOF; i++) {
		a1[i] = (2/pow(T,3))*(x[i]-init_joint_state[i+1].th) + (1/(T*T))*(init_joint_state[i+1].thd + x[i+DOF]);
		a2[i] = (3/(T*T))*(init_joint_state[i+1].th-x[i]) - (1/T)*(2*x[i+DOF] + init_joint_state[i+1].thd);
	}

	return;
}

/*
 * Calculate the extrema candidates for each joint (2*7 candidates in total)
 * For the striking polynomial
 * Clamp to [0,T]
 *
 */
void calc_strike_ext_cand(double *a1, double *a2, double T, double *joint_max_cand, double *joint_min_cand) {

	static double q0dot[DOF];
	int i;
	static double cand1, cand2;

	for (i = 0; i < DOF; i++) {
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
void calc_return_ext_cand(double *a1, double *a2, const double *x, double *joint_max_cand, double *joint_min_cand) {

	static double Tret = TIME2RETURN;
	int i;
	static double cand1, cand2;

	double T = x[2*DOF];

	for (i = 0; i < DOF; i++) {
		// find the extrema time candidates
		cand1 = fmin(Tret, fmax(0,(-a2[i] + sqrt(a2[i]*a2[i] - 3*a1[i]*x[i+DOF]))/(3*a1[i])));
		cand2 = fmin(Tret, fmax(0,(-a2[i] - sqrt(a2[i]*a2[i] - 3*a1[i]*x[i+DOF]))/(3*a1[i])));
		// find the joint extrema values at those candidate points
		cand1 = a1[i]*pow(cand1,3) + a2[i]*pow(cand1,2) + x[i+DOF]*cand1 + x[i];
		cand2 = a1[i]*pow(cand2,3) + a2[i]*pow(cand2,2) + x[i+DOF]*cand2 + x[i];
		joint_max_cand[i] = fmax(cand1,cand2);
		joint_min_cand[i] = fmin(cand1,cand2);
	}
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
	for (i = 0; i < DOF; i++) {
		x[i] = init_joint_state[i+1].th;
		x[i+DOF] = 0.0;
	}
	x[2*DOF] = 0.5;
}

/*
 * target and hitTime must be initialized in the lookup table
 * can be used for MPC
 */
void init_soln_to_last(double *x, SL_DJstate target[], double *hitTime) {

	// initialize first dof entries to q0
	int i;
	for (i = 0; i < DOF; i++) {
		x[i] = target[i+1].th;
		x[i+DOF] = target[i+1].thd;
	}
	x[2*DOF] = *hitTime;
}

/*
 * Finalize the solution and update target SL structure and hitTime value
 */
void finalize_soln(double *x, SL_DJstate target[], double *hitTime) {

	int i;
	for (i = 0; i < DOF; i++) {
		target[i+1].th = x[i];
		target[i+1].thd = x[i+DOF];
	}
	double optim_time = (get_time() - initTimeOpt)/1e6;
	*hitTime = x[2*DOF] - optim_time;
	printf("NLOPT took %f ms\n", (get_time() - initTimeOpt)/1e3);
	//print_optim_vec(x);
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
	static double kin_violation[EQ_CONSTR_DIM];
	static double lim_violation[INEQ_CONSTR_DIM]; // joint limit violations on strike and return
	joint_lim_ineq_constr(INEQ_CONSTR_DIM, lim_violation, OPTIM_DIM, x, grad, NULL);
	kinematics_eq_constr(EQ_CONSTR_DIM, kin_violation, OPTIM_DIM, x, grad, NULL);
	double cost = costfunc(OPTIM_DIM, x, grad, params);

	if (verbose) {
		// give info on solution vector
		print_optim_vec(x);
	    printf("Found minimum at f = %0.10g\n", cost);
		printf("Position constraint violation: [%.2f %.2f %.2f]\n",kin_violation[0],kin_violation[1],kin_violation[2]);
		printf("Velocity constraint violation: [%.2f %.2f %.2f]\n",kin_violation[3],kin_violation[4],kin_violation[5]);
		printf("Normal constraint violation: [%.2f %.2f %.2f]\n",kin_violation[6],kin_violation[7],kin_violation[8]);
	}

	int i;
	for (i = 0; i < EQ_CONSTR_DIM; i++) {
		if (kin_violation[i] > 5*SLACK) {
			printf("Equality constraints violated!\n");
			printf("Not updating trajectory!\n");
			return FALSE;
		}
	}

	for (i = 0; i < INEQ_CONSTR_DIM; i++) {
		if (lim_violation[i] > SLACK) {
			printf("Joint limit violated by %.2f on joint %d\n", lim_violation[i], i % DOF + 1);
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
void first_order_hold(double *ballPred, double *racketVel,
		double *racketNormal, double* T) {

	int i;
	int N = (int) ((*T)/TSTEP);
	double Tdiff = (*T) - N*TSTEP;
	static int iter;
	int Nmax = (int) TPRED/TSTEP;

	if (isnan(*T)) {
		printf("Warning: T value is nan! "
				"Setting ball and racket to first index and T = 0.5.\n");
		for (i = 1; i <= CART; i++) {
			ballPred[i-1] = ballPath[0].x[i];
			racketVel[i-1] = racketDes[0].xd[i];
			racketNormal[i-1] = racketDes[0].n[i];
		}
		*T = 0.5;
		return;
	}
	if ((*T) > TPRED) {
		printf("Warning: Extrapolation! T value is greater than TPRED = %f. "
				"Setting ball and racket to last index and T = 0.5\n", TPRED);
		for (i = 1; i <= CART; i++) {
			ballPred[i-1] = ballPath[Nmax].x[i];
			racketVel[i-1] = racketDes[Nmax].xd[i];
			racketNormal[i-1] = racketDes[Nmax].n[i];
		}
		*T = 0.5;
		return;
	}

	for (i = 1; i <= CART; i++) {
		if (N < Nmax) {
			ballPred[i-1] = ballPath[N].x[i] + (Tdiff/TSTEP) * (ballPath[N+1].x[i] - ballPath[N].x[i]);
			racketVel[i-1] = racketDes[N].xd[i] + (Tdiff/TSTEP) * (racketDes[N+1].xd[i] - racketDes[N].xd[i]);
			racketNormal[i-1] = racketDes[N].n[i] + (Tdiff/TSTEP) * (racketDes[N+1].n[i] - racketDes[N].n[i]);
		}
		else {
			ballPred[i-1] = ballPath[N].x[i];
			racketVel[i-1] = racketDes[N].xd[i];
			racketNormal[i-1] = racketDes[N].n[i];
		}
	}
}


/*
 * Prints the 2*DOF + 1 dimensional solution in user-friendly format
 */
void print_optim_vec(double *x) {

	int i;
	printf("qf = [");
	for (i = 0; i < DOF; i++) {
		printf("%.2f  ", x[i]);
	}
	printf("]\n");
	printf("qfdot = [");
	for (i = 0; i < DOF; i++) {
		printf("%.2f  ", x[i+DOF]);
	}
	printf("]\n");
	printf("T = %.2f\n", x[2*DOF]);
}

/*
 * Return time of day as micro seconds
 */
long get_time() {
	struct timeval tv;
	if (gettimeofday(&tv, (struct timezone *)0) == 0)
		return (tv.tv_sec*1000*1000 + tv.tv_usec);  //us

	return 0.;
}
