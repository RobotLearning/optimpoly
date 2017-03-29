/*
 ============================================================================
 Name        : optimpoly.c
 Author      : Okan
 Version     :
 Date        : 30/05/2016
 Description : Nonlinear optimization in C using the NLOPT library
 ============================================================================
 */

#include "constants.h"
#include "utils.h"
#include "kinematics.h"
#include "stdlib.h"
#include "math.h"
#include "optim.h"

// firsttime checking
static bool firsttime[3]; // TODO: remove this global var!

// termination
static double test_optim(const double *x, coptim *params, racketdes *racketdata, bool info);
static void finalize_soln(const double* x, optim * params, double time_elapsed);
static bool check_optim_result(const int res);

// optimization related methods
static double costfunc(unsigned n, const double *x, double *grad, void *my_func_data);
static void kinematics_eq_constr(unsigned m, double *result, unsigned n,
		                  const double *x, double *grad, void *f_data);
static void joint_limits_ineq_constr(unsigned m, double *result,
		                      unsigned n, const double *x, double *grad, void *data);

static void calc_strike_poly_coeff(const double *q0, const double *q0dot, const double *x,
		                    double *a1, double *a2);
static void calc_return_poly_coeff(const double *q0, const double *q0dot,
		                    const double *x, const double time2return,
		                    double *a1, double *a2);
static void calc_strike_extrema_cand(const double *a1, const double *a2, const double T,
		                      const double *q0, const double *q0dot,
							  double *joint_max_cand, double *joint_min_cand);
static void calc_return_extrema_cand(const double *a1, const double *a2,
		                      const double *x, const double time2return,
							  double *joint_max_cand, double *joint_min_cand);
static void init_last_soln(const optim * params, double x[OPTIM_DIM]);
static void init_rest_soln(const coptim * params, double x[OPTIM_DIM]);
static void first_order_hold(const racketdes* racketdata, const double T, double racket_pos[NCART],
		               double racket_vel[NCART], double racket_n[NCART]);
static void print_input_structs(coptim *coparams,
	      	  	  	  	  	   racketdes *racket,
							   optim * params);

/*
 * FIXED PLAYER
 *
 * NLOPT optimization routine for table tennis trajectory generation
 *
 * Returns the maximum of violations
 * Maximum of :
 * 1. kinematics equality constraint violations
 * 2. joint limit violations throughout trajectory
 *
 */
double nlopt_optim_fixed_run(coptim *coparams,
					      racketdes *racketdata,
						  optim *params) {

	firsttime[0] = firsttime[1] = firsttime[2] = true;

	//print_input_structs(coparams, racketdata, params);
	params->update = false;
	params->running = true;
	nlopt_opt opt;
	double x[OPTIM_DIM];
	double tol_eq[EQ_CONSTR_DIM];
	double tol_ineq[INEQ_CONSTR_DIM];
	const_vec(EQ_CONSTR_DIM,1e-2,tol_eq);
	const_vec(INEQ_CONSTR_DIM,1e-3,tol_ineq);
	init_rest_soln(coparams,x); //parameters are the initial joint positions q0
	//init_last_soln(params,x);
	// set tolerances equal to second argument //

	// LN = does not require gradients //
	opt = nlopt_create(NLOPT_LN_COBYLA, OPTIM_DIM);
	nlopt_set_xtol_rel(opt, 1e-2);
	nlopt_set_lower_bounds(opt, coparams->lb);
	nlopt_set_upper_bounds(opt, coparams->ub);
	nlopt_set_min_objective(opt, costfunc, coparams);
	nlopt_add_inequality_mconstraint(opt, INEQ_CONSTR_DIM, joint_limits_ineq_constr,
			                         coparams, tol_ineq);
	nlopt_add_equality_mconstraint(opt, EQ_CONSTR_DIM, kinematics_eq_constr,
			                         racketdata, tol_eq);

	double init_time = get_time();
	double past_time = 0.0;
	double minf; // the minimum objective value, upon return //
	int res; // error code
	double max_violation;

	if ((res = nlopt_optimize(opt, x, &minf)) < 0) {
		if (params->verbose)
			printf("NLOPT failed with exit code %d!\n", res);
	    past_time = (get_time() - init_time)/1e3;
	    max_violation = test_optim(x,coparams,racketdata,params->verbose);
	}
	else {
		past_time = (get_time() - init_time)/1e3;
		if (params->verbose) {
			printf("NLOPT success with exit code %d!\n", res);
			printf("NLOPT took %f ms\n", past_time);
			printf("Found minimum at f = %0.10g\n", minf);
		}
	    max_violation = test_optim(x,coparams,racketdata,params->verbose);
	    if (max_violation < 1e-2)
	    	finalize_soln(x,params,past_time);
	}
	params->running = false;
	if (params->verbose)
		check_optim_result(res);
	//nlopt_destroy(opt);
	return max_violation;
}

static void print_input_structs(coptim *coparams,
	      	  	  	  	  	   racketdes *racketdata,
							   optim * params) {

	for (int i = 0; i < NDOF; i++) {
		printf("q0[%d] = %f\n", i, coparams->q0[i]);
		printf("q0dot[%d] = %f\n", i, coparams->q0dot[i]);
		printf("lb[%d] = %f\n", i, coparams->lb[i]);
		printf("ub[%d] = %f\n", i, coparams->ub[i]);
	}
	printf("Tret = %f\n", coparams->time2return);
	for (int i = 0; i < NDOF; i++) {
		printf("qf[%d] = %f\n", i, params->qf[i]);
		printf("qfdot[%d] = %f\n", i, params->qfdot[i]);
	}
	printf("Thit = %f\n", params->T);

	print_mat_size("pos = ", racketdata->pos, NCART, 5);
	print_mat_size("vel = ", racketdata->vel, NCART, 5);
	print_mat_size("normal = ", racketdata->normal, NCART, 5);

}

/*
 * Give info about the optimization after termination
 *
 */
static bool check_optim_result(const int res) {

	bool flag = false;
	switch (res) {
	case NLOPT_SUCCESS:
		printf("Success!\n");
		flag = true;
		break;
	case NLOPT_STOPVAL_REACHED:
		printf("Optimization stopped because stopval (above) was reached.\n");
		flag = true;
		break;
	case NLOPT_FTOL_REACHED:
		printf("Optimization stopped because ftol_rel "
				"or ftol_abs (above) was reached.\n");
		flag = true;
		break;
	case NLOPT_XTOL_REACHED:
		flag = true;
		printf("Optimization stopped because xtol_rel or xtol_abs (above) was reached.\n");
		break;
	case NLOPT_MAXEVAL_REACHED:
		flag = true;
		printf("Optimization stopped because maxeval (above) was reached.\n");
		break;
	case NLOPT_MAXTIME_REACHED:
		flag = true;
		printf("Optimization stopped because maxtime (above) was reached.\n");
		break;
	case NLOPT_FAILURE:
		printf("Epic fail!\n");
		break;
	case NLOPT_INVALID_ARGS:
		printf("Invalid arguments (e.g. lower bounds are bigger than "
				"upper bounds, an unknown algorithm was specified, etcetera).\n");
		break;
	case NLOPT_OUT_OF_MEMORY:
		printf("Ran out of memory!\n");
		break;
	case NLOPT_ROUNDOFF_LIMITED:
		printf("Halted because roundoff errors limited progress."
			"(In this case, the optimization still typically returns a useful result.\n");
		break;
	case NLOPT_FORCED_STOP:
		printf("Halted because of a forced termination: "
				"the user called nlopt_force_stop(opt)"
				"on the optimization’s nlopt_opt object "
				"opt from the user’s objective function or constraints.\n");
		break;

	}
	return flag;
}

/*
 * Debug by testing the constraint violation of the solution vector
 *
 */
static double test_optim(const double *x, coptim * coparams, racketdes * racketdata, bool info) {

	// give info on constraint violation
	double *grad = 0;
	static double kin_violation[EQ_CONSTR_DIM];
	static double lim_violation[INEQ_CONSTR_DIM]; // joint limit violations on strike and return
	kinematics_eq_constr(EQ_CONSTR_DIM, kin_violation,
			             OPTIM_DIM, x, grad, racketdata);
	joint_limits_ineq_constr(INEQ_CONSTR_DIM, lim_violation,
			                 OPTIM_DIM, x, grad, coparams);
	double cost = costfunc(OPTIM_DIM, x, grad, coparams);

	if (info) {
		// give info on solution vector
		print_optim_vec(x);
		printf("f = %.2f\n",cost);
		printf("Position constraint violation: [%.2f %.2f %.2f]\n",kin_violation[0],kin_violation[1],kin_violation[2]);
		printf("Velocity constraint violation: [%.2f %.2f %.2f]\n",kin_violation[3],kin_violation[4],kin_violation[5]);
		printf("Normal constraint violation: [%.2f %.2f %.2f]\n",kin_violation[6],kin_violation[7],kin_violation[8]);
		for (int i = 0; i < INEQ_CONSTR_DIM; i++) {
			if (lim_violation[i] > 0.0)
				printf("Joint limit violated by %.2f on joint %d\n", lim_violation[i], i % NDOF + 1);
		}
	}

	return fmax(max_abs_array(kin_violation,EQ_CONSTR_DIM),
			    max_array(lim_violation,INEQ_CONSTR_DIM));
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
	params->update = true;
}

/*
 * Calculates the cost function for table tennis trajectory generation optimization
 * to find spline (3rd order strike+return) polynomials
 */
static double costfunc(unsigned n, const double *x, double *grad, void *my_func_params) {

	double a1[NDOF];
	double a2[NDOF];
	static double *q0dot; // initial joint velocity
	static double *q0; // initial joint pos
	double T = x[2*NDOF];

	if (firsttime[0]) {
		firsttime[0] = false;
		coptim* optim_data = (coptim*)my_func_params;
		q0 = optim_data->q0;
		q0dot = optim_data->q0dot;
	}

	// calculate the polynomial coeffs which are used in the cost calculation
	calc_strike_poly_coeff(q0,q0dot,x,a1,a2);

	return T * (3*T*T*inner_prod(NDOF,a1,a1) +
			3*T*inner_prod(NDOF,a1,a2) + inner_prod(NDOF,a2,a2));
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

	if (firsttime[1]) {
		firsttime[1] = false;
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
 * This is the constraint that makes sure we hit the ball
 */
static void kinematics_eq_constr(unsigned m, double *result, unsigned n,
		                  const double *x, double *grad, void *my_function_data) {

	static double racket_des_pos[NCART];
	static double racket_des_vel[NCART];
	static double racket_des_normal[NCART];
	static double pos[NCART];
	static double qfdot[NDOF];
	static double vel[NCART];
	static double normal[NCART];
	static double qf[NDOF];
	static racketdes* racket_data;
	double T = x[2*NDOF];

	/* initialization of static variables */
	if (firsttime[2]) {
		firsttime[2] = false;
		racket_data = (racketdes*) my_function_data;
	}

	// interpolate at time T to get the desired racket parameters
	first_order_hold(racket_data,T,racket_des_pos,racket_des_vel,racket_des_normal);

	// extract state information from optimization variables
	for (int i = 0; i < NDOF; i++) {
		qf[i] = x[i];
		qfdot[i] = x[i+NDOF];
	}

	// compute the actual racket pos,vel and normal
	calc_racket_state(qf,qfdot,pos,vel,normal);

	// deviations from the desired racket frame
	for (int i = 0; i < NCART; i++) {
		result[i] = pos[i] - racket_des_pos[i];
		result[i + NCART] = vel[i] - racket_des_vel[i];
		result[i + 2*NCART] = normal[i] - racket_des_normal[i];
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
static void first_order_hold(const racketdes* racketdata, const double T, double racket_pos[NCART],
		               double racket_vel[NCART], double racket_n[NCART]) {

	double deltat = racketdata->dt;
	if (isnan(T)) {
		printf("Warning: T value is nan!\n");

		for(int i = 0; i < NCART; i++) {
			racket_pos[i] = racketdata->pos[i][0];
			racket_vel[i] = racketdata->vel[i][0];
			racket_n[i] = racketdata->normal[i][0];
		}
	}
	else {
		int N = (int) (T/deltat);
		double Tdiff = T - N*deltat;
		int Nmax = racketdata->Nmax;

		for (int i = 0; i < NCART; i++) {
			if (N < Nmax) {
				racket_pos[i] = racketdata->pos[i][N] +
						(Tdiff/deltat) * (racketdata->pos[i][N+1] - racketdata->pos[i][N]);
				racket_vel[i] = racketdata->vel[i][N] +
						(Tdiff/deltat) * (racketdata->vel[i][N+1] - racketdata->vel[i][N]);
				racket_n[i] = racketdata->normal[i][N] +
						(Tdiff/deltat) * (racketdata->normal[i][N+1] - racketdata->normal[i][N]);
			}
			else {
				racket_pos[i] = racketdata->pos[i][N];
				racket_vel[i] = racketdata->vel[i][N];
				racket_n[i] = racketdata->normal[i][N];
			}
		}
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
static void calc_return_poly_coeff(const double *q0, const double *q0dot,
		                    const double *x, const double T,
		                    double *a1, double *a2) {

	for (int i = 0; i < NDOF; i++) {
		a1[i] = (2/pow(T,3))*(x[i]-q0[i]) + (1/(T*T))*(q0dot[i] + x[i+NDOF]);
		a2[i] = (3/(T*T))*(q0[i]-x[i]) - (1/T)*(2*x[i+NDOF] + q0dot[i]);
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
static void calc_return_extrema_cand(const double *a1, const double *a2,
		                      const double *x, const double Tret,
		                      double *joint_max_cand, double *joint_min_cand) {

	static double cand1, cand2;

	for (int i = 0; i < NDOF; i++) {
		cand1 = fmin(Tret, fmax(0,(-a2[i] + sqrt(a2[i]*a2[i] - 3*a1[i]*x[i+NDOF]))/(3*a1[i])));
		cand2 =  fmin(Tret, fmax(0,(-a2[i] - sqrt(a2[i]*a2[i] - 3*a1[i]*x[i+NDOF]))/(3*a1[i])));
		cand1 = a1[i]*pow(cand1,3) + a2[i]*pow(cand1,2) + x[i+NDOF]*cand1 + x[i];
		cand2 = a1[i]*pow(cand2,3) + a2[i]*pow(cand2,2) + x[i+NDOF]*cand2 + x[i];
		joint_max_cand[i] = fmax(cand1,cand2);
		joint_min_cand[i] = fmin(cand1,cand2);
	}
}

/*
 * Initialize solution for NLOPT using last solution of
 * 2*dof + 1 dimensional problem
 *
 * The closer to the optimum it is the faster alg should converge
 */
static void init_last_soln(const optim * params, double x[OPTIM_DIM]) {

	// initialize first dof entries to q0
	for (int i = 0; i < NDOF; i++) {
		x[i] = params->qf[i];
		x[i+NDOF] = params->qfdot[i];
	}
	x[2*NDOF] = params->T;
}

/*
 * Initialize solution for NLOPT using rest posture
 * 2*dof + 1 dimensional problem
 *
 * Sets T to 0.5
 *
 * The closer to the optimum it is the faster alg should converge
 */
static void init_rest_soln(const coptim * params, double x[OPTIM_DIM]) {

	// initialize first dof entries to q0
	for (int i = 0; i < NDOF; i++) {
		x[i] = params->qrest[i];
		x[i+NDOF] = 0.0;
	}
	x[2*NDOF] = 0.5;
}
