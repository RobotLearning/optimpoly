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
#include "SL_kinematics_body.h"

int main(void) {

	Vector ballLand = my_vector(1,CART);
	static double landTime;

	// initialize ball
	// predict for T_pred seconds
	set_des_land_param(ballLand,&landTime);
	init_ball_state();
	init_joint_state();
	predict_ball_state();
	calc_racket_strategy(ballLand,landTime);

	//double initTime = get_time();
	//nlopt_example_run();
	optim_poly_nlopt_run();
	//printf("NLOPT took %f ms\n", (get_time() - initTime)/1e3);

	// test the lookup value to see if constraint is not violated
	printf("================== TEST ==================\n");
	printf("Lookup values:\n");
	static double x[OPTIM_DIM];
	lookup(x);
	test_constraint(x);

	return TRUE;
}


/*
 *
 * Set the ball values to a reasonable value
 * SO FAR setting it into the first lookup table entry from March 2016
 *
 * Using the ballPred structure from table_tennis_common.h
 */
void init_ball_state() {

	bzero((char *)&(ballPred), sizeof(ballPred));
	// initialize the ball

	ballPred.x[1] = 0.1972;
	ballPred.x[2] = -2.4895;
	ballPred.x[3] = -0.5040;
	ballPred.xd[1] = -1.7689;
	ballPred.xd[2] = 4.7246;
	ballPred.xd[3] = -1.0867;
}

/*
 * Set the initial posture of the robot
 */
void init_joint_state() {

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
 *
 * Test the solution found by using a lookup table
 * TODO: expand to actually use the saved lookup table
 */
void lookup(double *x) {

	// load an x value

	// qf values
//	x[0] = 0.5580;
//	x[1] = 0.2266;
//	x[2] = 0.0179;
//	x[3] = 1.6754;
//	x[4] = -1.3887;
//	x[5] = -0.8331;
//	x[6] = 0.3118;
//	// qfdot values
//	x[7] = -1.8601;
//	x[8] = 2.1229;
//	x[9] = -0.4704;
//	x[10] = 0.2180;
//	x[11] = 0.2867;
//	x[12] = -1.7585;
//	x[13] = 0.0281;
//	// T values
//	x[14] = 0.6300;

	x[0] = 0.4477;
	x[1] = 0.4065;
    x[2] = -0.0825;
    x[3] = 1.6478;
	x[4] = -1.2868;
	x[5] = -1.0515;
    x[6] = 0.1807;
	// qfdot values
	x[7] = -1.6092;
	x[8] = 2.9100;
	x[9] = -0.2404;
	x[10] = -0.5597;
	x[11] = 0.9379;
	x[12] = -2.4477;
	x[13] = -0.3082;
	// T values
	x[14] = 0.5800;

}

/*
 * NLOPT optimization routine for table tennis traj gen
 */
void optim_poly_nlopt_run() {

	static double tol[CONSTR_DIM];
	static double lb[OPTIM_DIM]; /* lower bounds */
	static double ub[OPTIM_DIM]; /* upper bounds */
	static double x[OPTIM_DIM]; /* initial guess */

	load_joint_limits();
	set_bounds(lb,ub);
	const_vec(CONSTR_DIM,1e-2,tol);

	nlopt_opt opt;
	opt = nlopt_create(NLOPT_LN_COBYLA, OPTIM_DIM); /* LN = does not require gradients */
	//nlopt_set_xtol_rel(opt, 1e-2);
	//opt = nlopt_create(NLOPT_AUGLAG, OPTIM_DIM); /* algorithm and dimensionality */
	//nlopt_set_local_optimizer(opt, opt);
	nlopt_set_lower_bounds(opt, lb);
	nlopt_set_upper_bounds(opt, ub);
	nlopt_set_min_objective(opt, costfunc, NULL);
	nlopt_add_equality_mconstraint(opt, CONSTR_DIM, kinematics_constr, NULL, tol);
	nlopt_set_xtol_rel(opt, 1e-2);

	//int maxeval = 20000;
	//nlopt_set_maxeval(opt, maxeval);

	//double maxtime = 5.0;
	//nlopt_set_maxtime(opt, maxtime);

	init_soln(x);
	double initTime = get_time();
	double minf; /* the minimum objective value, upon return */

	if (nlopt_optimize(opt, x, &minf) < 0) {
	    printf("NLOPT failed!\n");
	}
	else {
		//nlopt_example_run();
		printf("NLOPT took %f ms\n", (get_time() - initTime)/1e3);
	    printf("Found minimum at f = %0.10g\n", minf);
	    test_constraint(x);
	}
	nlopt_destroy(opt);
}

/*
 * Debug by testing the constraint violation of the solution vector
 *
 */
void test_constraint(double *x) {
	// give info on solution vector
	print_optim_vec(x);
	// give info on constraint violation
	double *grad = FALSE;
	static double violation[CONSTR_DIM];
	kinematics_constr(CONSTR_DIM, violation, OPTIM_DIM, x, grad, NULL);
	double cost = costfunc(OPTIM_DIM, x, grad, NULL);
	printf("f = %.2f\n",cost);
	printf("Position constraint violation: [%.2f %.2f %.2f]\n",violation[0],violation[1],violation[2]);
	printf("Velocity constraint violation: [%.2f %.2f %.2f]\n",violation[3],violation[4],violation[5]);
	printf("Normal constraint violation: [%.2f %.2f %.2f]\n",violation[6],violation[7],violation[8]);
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
 * Copied from SL_common. Dont want to call the function from SL because
 * we have to include a lot of extra SL files
 *
 */
int read_sensor_offsets(char *fname) {

  int j,i,rc;
  char   string[100];
  FILE  *in;

  /* get the max, min, and offsets of the position sensors */

  sprintf(string,"%s/robolab/barrett/%s%s",getenv("HOME"),CONFIG,fname);
  in = fopen(string,"r");
  if (in == NULL) {
    printf("ERROR: Cannot open file >%s<!\n",string);
    return FALSE;
  }

  /* find all joint variables and read them into the appropriate array */

  for (i=1; i<= n_dofs; ++i) {
    if (!find_keyword(in, &(joint_names[i][0]))) {
      printf("ERROR: Cannot find offset for >%s<!\n",joint_names[i]);
      fclose(in);
      return FALSE;
    }
    rc=fscanf(in,"%lf %lf %lf %lf %lf %lf",
	&joint_range[i][MIN_THETA], &joint_range[i][MAX_THETA],
	   &(joint_default_state[i].th),
	   &(joint_opt_state[i].th),
	   &(joint_opt_state[i].w),
	   &joint_range[i][THETA_OFFSET]);
    joint_default_state[i].thd = 0;
    joint_default_state[i].uff = 0;
  }

  fclose(in);

  return TRUE;

}

/*
 *
 * Copied from SL_common. Dont want to call the function from SL because
 * we have to include a lot of extra SL files
 *

 computes one column for the geometric jacobian of a revolute joint
 from the given input vectors

 *******************************************************************************
 Function Parameters: [in]=input,[out]=output

 \param[in]     p    : position of endeffector
 \param[in]     pi   : position of joint origin
 \param[in]     zi   : unit vector of joint axis
 \param[out]    c    : column vector of Jacobian

 ******************************************************************************/
void revoluteGJacColumn(Vector p, Vector pi, Vector zi, Vector c) {
  int i,j;

  c[1] = zi[2] * (p[3]-pi[3]) - zi[3] * (p[2]-pi[2]);
  c[2] = zi[3] * (p[1]-pi[1]) - zi[1] * (p[3]-pi[3]);
  c[3] = zi[1] * (p[2]-pi[2]) - zi[2] * (p[1]-pi[1]);
  c[4] = zi[1];
  c[5] = zi[2];
  c[6] = zi[3];

}

/*
 * Estimate an initial solution for NLOPT
 * 2*dof + 1 dimensional problem
 *
 * The closer to the optimum it is the faster alg should converge
 * TODO: load values from a lookup table
 */
void init_soln(double * x) {

//	x[0] = 0.45;
//	x[1] = 0.41;
//	x[2] = -0.08;
//	x[3] = 1.65;
//	x[4] = -1.29;
//	x[5] = -1.05;
//	x[6] = 0.18;
//
//	x[7] = -1.61;
//	x[8] = 2.91;
//	x[9] = -0.24;
//	x[10] = -0.56;
//	x[11] = 0.94;
//	x[12] = -2.45;
//	x[13] = -0.31;

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
double costfunc(unsigned n, const double *x, double *grad, void *my_func_data) {

	int i;
	double a1[DOF];
	double a2[DOF];
	static double q0dot[DOF]; // all zeros
	static double qfdot[DOF]; // all zeros
	static double qf[DOF]; // opt value
	double T = x[2*DOF];

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
	calc_poly_coeff(a1,a2,q0,x);

	return T * (3*T*T*inner_prod(a1,a1) +
			3*T*inner_prod(a1,a2) + inner_prod(a2,a2));
}

/*
 * This is the constraint that makes sure we hit the ball
 */
void kinematics_constr(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data) {

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
	linkInformationDes(joint_des_state, &base_state, &base_orient, endeff,
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
 * Copied from SL_user_common.c for convenience
 * TODO: is this correct?
 *
 */
void setDefaultEndeffector(void) {


	endeff[RIGHT_HAND].m       = 0.0;
	endeff[RIGHT_HAND].mcm[_X_]= 0.0;
	endeff[RIGHT_HAND].mcm[_Y_]= 0.0;
	endeff[RIGHT_HAND].mcm[_Z_]= 0.0;
	endeff[RIGHT_HAND].x[_X_]  = 0.0;
	endeff[RIGHT_HAND].x[_Y_]  = 0.0;
	endeff[RIGHT_HAND].x[_Z_]  = 0.08531+0.06;
	endeff[RIGHT_HAND].a[_A_]  = 0.0;
	endeff[RIGHT_HAND].a[_B_]  = 0.0;
	endeff[RIGHT_HAND].a[_G_]  = 0.0;

	// attach the racket
	endeff[RIGHT_HAND].x[_Z_] = .3;

}

/*
 * Calculate the polynomial coefficients from the optimized variables qf,qfdot,T
 * p(t) = a1*t^3 + a2*t^2 + a3*t + a4
 */
void calc_poly_coeff(double *a1, double *a2, const double *q0, const double *x) {

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
