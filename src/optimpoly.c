/*
 ============================================================================
 Name        : optimpoly.c
 Author      : Okan
 Version     :
 Date        : 30/05/2016
 Description : Nonlinear optimization in C using the NLOPT library
 ============================================================================
 */

#include "optimpoly.h"
#include "table_tennis.h"
#include "extra.h"
#include "kinematics.h"
//#include "SL_kinematics_body.h"

int main(void) {

	//nlopt_example_run();
	//pthread_example_run();
	//return TRUE;

	Vector ballLand = my_vector(1,CART);
	double landTime;
	double q0[DOF];
	double b0[CART];
	double v0[CART];
	double x[OPTIM_DIM]; /* initial guess for optim */

	Matrix lookupTable = my_matrix(1, LOOKUP_TABLE_SIZE, 1, LOOKUP_COLUMN_SIZE);
	load_lookup_table(lookupTable);
	load_joint_limits();
	set_des_land_param(ballLand,&landTime);

	/* initialize ball and racket */
	// predict for T_pred seconds
	init_ball_state(b0,v0);
	init_joint_state(q0);
	predict_ball_state(b0,v0);
	calc_racket_strategy(ballLand,landTime);

	/* run NLOPT opt algorithm here */
	lookup(lookupTable,b0,v0,x); // initialize solution
	nlopt_optim_poly_run(x,q0);

	// test the lookup value to see if constraint is not violated
	printf("================== TEST ==================\n");
	printf("Lookup values:\n");
	lookup(lookupTable,b0,v0,x);
	test_optim(x,q0);

	return TRUE;
}

/*
 *
 * Loading lookup table to initialize the optimization
 * And to test constraints
 */
void load_lookup_table(Matrix lookupTable) {

	printf("Loading lookup table...\n");

	load_vec_into_mat(lookupTable, LOOKUP_TABLE_SIZE, LOOKUP_COLUMN_SIZE, LOOKUP_TABLE_NAME);
	//print_mat("Lookup: ", lookupTable);
	/*int i;
	for(i = 1; i <= LOOKUP_COLUMN_SIZE; i++)
		printf("%.2f\t",lookupTable[1][i]);
	printf("\n");*/
}


/*
 *
 * Set the ball values to a reasonable value
 * SO FAR setting it into the first lookup table entry from May 2016
 *
 */
void init_ball_state(double *b0, double *v0) {

	// initialize the ball
	b0[0] = 0.1972;
	b0[1] = -2.4895;
	b0[2] = -0.5040;
	v0[0] = -1.7689;
	v0[1] = 4.7246;
	v0[2] = -1.0867;
}

/*
 * Set the initial posture of the robot
 */
void init_joint_state(double *q0) {

	// initialize the variables
	//q0 = [1.0; -0.2; -0.1; 1.8; -1.57; 0.1; 0.3];
	q0[0] = 1.0;
	q0[1] = -0.2;
	q0[2] = -0.1;
	q0[3] = 1.8;
	q0[4] = -1.57;
	q0[5] = 0.1;
	q0[6] = 0.3;
}

/*
 * NLOPT optimization routine for table tennis traj gen
 *
 * TODO: adding inequality constraints makes problem harder. Is it necessary always?
 */
void nlopt_optim_poly_run(double *x, double *params) {

	static double tol[EQ_CONSTR_DIM];
	static double lb[OPTIM_DIM]; /* lower bounds */
	static double ub[OPTIM_DIM]; /* upper bounds */

	set_bounds(lb,ub);
	const_vec(EQ_CONSTR_DIM,1e-2,tol); /* set tolerances equal to second argument */

	nlopt_opt opt;
	opt = nlopt_create(NLOPT_LN_COBYLA, OPTIM_DIM); /* LN = does not require gradients */
	//nlopt_set_xtol_rel(opt, 1e-2);
	//opt = nlopt_create(NLOPT_AUGLAG, OPTIM_DIM); /* algorithm and dimensionality */
	//nlopt_set_local_optimizer(opt, opt);
	nlopt_set_lower_bounds(opt, lb);
	nlopt_set_upper_bounds(opt, ub);
	nlopt_set_min_objective(opt, costfunc, params);
	nlopt_add_inequality_mconstraint(opt, INEQ_CONSTR_DIM, joint_limits_ineq_constr, params, tol);
	nlopt_add_equality_mconstraint(opt, EQ_CONSTR_DIM, kinematics_eq_constr, params, tol);
	nlopt_set_xtol_rel(opt, 1e-2);

	//int maxeval = 20000;
	//nlopt_set_maxeval(opt, maxeval);

	//double maxtime = 0.001;
	//nlopt_set_maxtime(opt, maxtime);

	//init_soln_to_rest_posture(x,params); //parameters are the initial joint positions q0
	double initTime = get_time();
	double minf; /* the minimum objective value, upon return */

	if (nlopt_optimize(opt, x, &minf) < 0) {
	    printf("NLOPT failed!\n");
	}
	else {
		//nlopt_example_run();
		printf("NLOPT took %f ms\n", (get_time() - initTime)/1e3);
	    printf("Found minimum at f = %0.10g\n", minf);
	    test_optim(x,params);
	}
	nlopt_destroy(opt);
}

/*
 * Debug by testing the constraint violation of the solution vector
 *
 */
void test_optim(double *x, double *params) {
	// give info on solution vector
	print_optim_vec(x);
	// give info on constraint violation
	double *grad = FALSE;
	static double kin_violation[EQ_CONSTR_DIM];
	static double lim_violation[INEQ_CONSTR_DIM]; // joint limit violations on strike and return
	kinematics_eq_constr(EQ_CONSTR_DIM, kin_violation, OPTIM_DIM, x, grad, NULL);
	joint_limits_ineq_constr(INEQ_CONSTR_DIM, lim_violation, OPTIM_DIM, x, grad, params);
	double cost = costfunc(OPTIM_DIM, x, grad, params);
	printf("f = %.2f\n",cost);
	printf("Position constraint violation: [%.2f %.2f %.2f]\n",kin_violation[0],kin_violation[1],kin_violation[2]);
	printf("Velocity constraint violation: [%.2f %.2f %.2f]\n",kin_violation[3],kin_violation[4],kin_violation[5]);
	printf("Normal constraint violation: [%.2f %.2f %.2f]\n",kin_violation[6],kin_violation[7],kin_violation[8]);

	int i;
	for (i = 0; i < INEQ_CONSTR_DIM; i++) {
		if (lim_violation[i] > 0.0)
			printf("Joint limit violated by %.2f on joint %d\n", lim_violation[i], i % DOF + 1);
	}
}

/*
 * Set upper and lower bounds on the optimization
 */
void set_bounds(double *lb, double *ub) {

	// lower bounds and upper bounds for qf are the joint limits
	int i;
	for (i = 1; i <= DOF; i++) {
		ub[i-1] = joint_range[i][MAX_THETA];
		lb[i-1] = joint_range[i][MIN_THETA];
		ub[i-1+DOF] = MAX_VEL;
		lb[i-1+DOF] = -MAX_VEL;
	}
	// constraints on final time
	ub[2*DOF] = 1.0;
	lb[2*DOF] = 0.0;
}

/*
 * Load the joint limits from config/ file into joint_range array
 *
 */
void load_joint_limits() {

	char *fname = "SensorOffset.cf";
	read_sensor_offsets(fname);

}

/*
 * Estimate an initial solution for NLOPT
 * 2*dof + 1 dimensional problem
 *
 * The closer to the optimum it is the faster alg should converge
 */
void init_soln_to_rest_posture(double *x, double *q0) {

	// initialize first dof entries to q0
	int i;
	for (i = 0; i < DOF; i++) {
		x[i] = q0[i];
		x[i+DOF] = 0.0;
	}
	x[2*DOF] = 0.58;
}

/*
 * Calculates the cost function for table tennis trajectory generation optimization
 * to find spline (3rd order strike+return) polynomials
 */
double costfunc(unsigned n, const double *x, double *grad, void *my_func_params) {

	int i;
	double a1[DOF];
	double a2[DOF];
	static double q0dot[DOF]; // all zeros
	static double qfdot[DOF]; // all zeros
	static double qf[DOF]; // opt value
	double T = x[2*DOF];

	double *q0 = (double *) my_func_params;
	// instead of using global variables feeding parameters directly

	if (grad) {

		for (i = 0; i < DOF; i++) {
			qf[i] = x[i];
			qfdot[i] = x[i+DOF];
			grad[i] = (6/pow(T,3))*(qf[i]-q0[i]) - (3/(T*T))*(q0dot[i]+qfdot[i]);
			grad[i+DOF] = (-3/(T*T))*(qf[i]-q0[i]) + (1/T)*(2*qfdot[i]+q0dot[i]);
		}
		//time derivative of cost J
		vec_minus(qf,q0);
		//assuming q0dot = 0 here;
		grad[2*DOF] = (-9/pow(T,4))*inner_prod(qf,qf) +
					  (6/pow(T,3))*inner_prod(qf,qfdot) +
					  (3/(T*T))*inner_prod(q0dot,qfdot) -
					  (1/(T*T))*inner_prod(qfdot,qfdot);
	}



	// calculate the polynomial coeffs which are used in the cost calculation
	calc_strike_poly_coeff(a1,a2,q0,x);

	return T * (3*T*T*inner_prod(a1,a1) +
			3*T*inner_prod(a1,a2) + inner_prod(a2,a2));
}

/*
 * This is the inequality constraint that makes sure we never exceed the
 * joint limits during the striking and returning motion
 *
 */
void joint_limits_ineq_constr(unsigned m, double *result, unsigned n, const double *x, double *grad, void *my_func_params) {

	int i;
	static double a1[DOF];
	static double a2[DOF];
	static double a1ret[DOF]; // coefficients for the returning polynomials
	static double a2ret[DOF];
	static double joint_strike_max_cand[DOF];
	static double joint_strike_min_cand[DOF];
	static double joint_return_max_cand[DOF];
	static double joint_return_min_cand[DOF];

	double *q0 = (double *) my_func_params;

	// calculate the polynomial coeffs which are used for checking joint limits
	calc_strike_poly_coeff(a1,a2,q0,x);
	calc_return_poly_coeff(a1ret,a2ret,q0,x);
	// calculate the candidate extrema both for strike and return
	calc_strike_extrema_cand(a1,a2,x[2*DOF],q0,joint_strike_max_cand,joint_strike_min_cand);
	calc_return_extrema_cand(a1ret,a2ret,x,joint_return_max_cand,joint_return_min_cand);

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
 * This is the constraint that makes sure we hit the ball
 */
void kinematics_eq_constr(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data) {

	static double ballPred[CART];
	static double racketDesVel[CART];
	static double racketDesNormal[CART];
	static Matrix  link_pos_des;
	static Matrix  joint_cog_mpos_des;
	static Matrix  joint_origin_pos_des;
	static Matrix  joint_axis_pos_des;
	static Matrix  Alink_des[N_LINKS+1];
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

		link_pos_des         = my_matrix(1,N_LINKS,1,3);
		joint_cog_mpos_des   = my_matrix(1,N_DOFS,1,3);
		joint_origin_pos_des = my_matrix(1,N_DOFS,1,3);
		joint_axis_pos_des   = my_matrix(1,N_DOFS,1,3);
		Jacobi               = my_matrix(1,2*CART,1,N_DOFS);
		racketTransform      = my_matrix(1,4,1,4);
		qfdot                = my_vector(1,DOF);
		xfdot                = my_vector(1,2*CART);
		normal               = my_vector(1,CART);

		for (i = 0; i <= N_LINKS; ++i)
			Alink_des[i] = my_matrix(1,4,1,4);

		// initialize the base variables
		//taken from ParameterPool.cf
		bzero((void *) &base_state, sizeof(base_state));
		bzero((void *) &base_orient, sizeof(base_orient));
		base_orient.q[_Q1_] = 1.0;

		// the the default endeffector parameters
		setDefaultEndeffector();

		// homogeneous transform instead of using quaternions
		racketTransform[1][1] = 1;
		racketTransform[2][3] = 1;
		racketTransform[3][2] = -1;
		racketTransform[4][4] = 1;
	}

	if (grad) {
		// compute gradient of kinematics = jacobian
		//TODO:
		grad[0] = 0.0;
		grad[1] = 0.0;
		grad[2] = 0.0;
	}

	// interpolate at time T to get the desired racket parameters
	first_order_hold(ballPred,racketDesVel,racketDesNormal,T);

	// extract state information from array to joint_des_state structure
	for (i = 1; i <= DOF; i++) {
		joint_des_state[i].th = x[i-1];
		joint_des_state[i].thd = qfdot[i] = x[i-1+DOF];
	}

	/* compute the desired link positions */
	kinematics(joint_des_state, &base_state, &base_orient, endeff,
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
 * Calculate the polynomial coefficients from the optimized variables qf,qfdot,T
 * p(t) = a1*t^3 + a2*t^2 + a3*t + a4
 */
void calc_strike_poly_coeff(double *a1, double *a2, const double *q0, const double *x) {

	// variables to be optimized
	static double q0dot[DOF];
	int i;
	double T = x[2*DOF];

	for (i = 0; i < DOF; i++) {
		a1[i] = (2/pow(T,3))*(q0[i]-x[i]) + (1/(T*T))*(q0dot[i] + x[i+DOF]);
		a2[i] = (3/(T*T))*(x[i]-q0[i]) - (1/T)*(x[i+DOF] + 2*q0dot[i]);
	}

	return;
}

/*
 * Calculate the returning polynomial coefficients from the optimized variables qf,qfdot
 * and time2return variable defined in the header
 * p(t) = a1*t^3 + a2*t^2 + a3*t + a4
 */
void calc_return_poly_coeff(double *a1, double *a2, const double *q0, const double *x) {

	// variables to be optimized
	static double q0dot[DOF];
	static double T = TIME2RETURN;
	int i;

	for (i = 0; i < DOF; i++) {
		a1[i] = (2/pow(T,3))*(x[i]-q0[i]) + (1/(T*T))*(q0dot[i] + x[i+DOF]);
		a2[i] = (3/(T*T))*(q0[i]-x[i]) - (1/T)*(2*x[i+DOF] + q0dot[i]);
	}

	return;
}

/*
 * Calculate the extrema candidates for each joint (2*7 candidates in total)
 * For the striking polynomial
 * Clamp to [0,T]
 * TODO: discard if complex
 *
 */
void calc_strike_extrema_cand(double *a1, double *a2, double T, double *q0,
		                      double *joint_max_cand, double *joint_min_cand) {

	static double q0dot[DOF];
	int i;
	static double cand1, cand2;

	for (i = 0; i < DOF; i++) {
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
 * Clamp to [T,T + TIME2RETURN]
 * TODO: discard if complex
 *
 */
void calc_return_extrema_cand(double *a1, double *a2, const double *x, double *joint_max_cand, double *joint_min_cand) {

	static double Tret = TIME2RETURN;
	int i;
	static double cand1, cand2;

	double T = x[2*DOF];

	for (i = 0; i < DOF; i++) {
		cand1 = fmin(T + Tret, fmax(T,(-a2[i] + sqrt(a2[i]*a2[i] - 3*a1[i]*x[i+DOF]))/(3*a1[i])));
		cand2 =  fmin(T + Tret, fmax(0,(-a2[i] - sqrt(a2[i]*a2[i] - 3*a1[i]*x[i+DOF]))/(3*a1[i])));
		cand1 = a1[i]*pow(cand1,3) + a2[i]*pow(cand1,2) + x[i+DOF]*cand1 + x[i];
		cand2 = a1[i]*pow(cand2,3) + a2[i]*pow(cand2,2) + x[i+DOF]*cand2 + x[i];
		joint_max_cand[i] = fmax(cand1,cand2);
		joint_min_cand[i] = fmin(cand1,cand2);
	}
}
