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
#include "SL_kinematics_body.h"

// initialization needs to be done for this mapping - used in SL_kinematics_body.h
int  link2endeffmap[] = {0,PALM};

int main(void) {

	Vector ballLand = my_vector(1,CART);
	static double landTime;

	// initialize ball
	// predict for T_pred seconds
	set_land_parameters(ballLand,&landTime);
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
 * Set desired landing position to the centre of the table and time to a
 * reasonable value
 */
void set_land_parameters(Vector ballLand, double *landTime) {

	*landTime = 0.8;

	ballLand[_X_] = 0.0;
	ballLand[_Y_] = dist_to_table - 3*table_length/4; // centre of opponents court
	ballLand[_Z_] = floor_level - table_height;// + ball_radius;
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
 * Predict the ball for Tpred seconds
 * Using ballflight + bounce model
 */
void predict_ball_state() {

	int N = TPRED/TSTEP;
	ballMat = my_matrix(1, N, 1, 2*CART);
	int i,j;

	for (j = 1; j <= CART; j++) {
		ballMat[0][j] = ballPred.x[j];
		ballMat[0][j+CART] = ballPred.xd[j];
	}

	// predict Tpred seconds into the future
	int bounce = FALSE;

	for (i = 1; i <= N; i++) {
		integrateBallState(ballPred,&ballPred,TSTEP,&bounce);
		for (j = 1; j <= CART; j++) {
			ballMat[i][j] = ballPred.x[j];
			ballMat[i][j+CART] = ballPred.xd[j];
		}
	}

	//print_mat("Ball pred matrix: ", ballMat);
	/*for (i = 1; i <= N; i++) {
		for (j = 1; j <= 2*CART; j++)
			printf("%.4f  ", ballMat[i][j]);
		printf("\n");
	}*/

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

	guesstimate_soln(x);
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
 * Integrate the ball state dt seconds later. Used for prediction
 * Using symplectic Euler.
 *
 */
void integrateBallState(SL_Cstate ballState, SL_Cstate *ballPred, double deltat, int *bounce) {

	int i;
	double velBall;
	static double slack = 0.0001;

	// does it touch the floor?
	if (ballState.x[_Z_] >= floor_level) {

		/*****************************************/
		// Symplectic Euler for No-Contact-Situation (Flight model)

		velBall = sqrt(sqr(ballState.xd[_X_]) + sqr(ballState.xd[_Y_]) + sqr(ballState.xd[_Z_]));

		if (fabs(velBall) > slack) {

			ballPred->xdd[_X_] = -ballState.xd[_X_] * Cdrag * velBall;
			ballPred->xdd[_Y_] = -ballState.xd[_Y_] * Cdrag * velBall;
			ballPred->xdd[_Z_] = -intern_gravity - ballState.xd[_Z_] * Cdrag * velBall;
		}
		else {
			ballPred->xdd[_X_] = 0.;
			ballPred->xdd[_Y_] = 0.;
			ballPred->xdd[_Z_] = -intern_gravity;
		}

		ballPred->xd[_X_] = ballState.xd[_X_] + ballPred->xdd[_X_] * deltat;
		ballPred->xd[_Y_] = ballState.xd[_Y_] + ballPred->xdd[_Y_] * deltat;
		ballPred->xd[_Z_] = ballState.xd[_Z_] + ballPred->xdd[_Z_] * deltat;
		ballPred->x[_X_] =  ballState.x[_X_] + ballPred->xd[_X_] * deltat;
		ballPred->x[_Y_] =  ballState.x[_Y_] + ballPred->xd[_Y_] * deltat;
		ballPred->x[_Z_] =  ballState.x[_Z_] + ballPred->xd[_Z_] * deltat;

		if (checkForBallTableContact(*ballPred)) {

			//printf("Expecting a bounce!\n");
			ballPred->x[_Z_]  = floor_level - table_height  + ball_radius;
			ballPred->xd[_Z_] = -CRT * ballPred->xd[_Z_];
			ballPred->xd[_Y_] = CFTY * ballPred->xd[_Y_];
			ballPred->xd[_X_] = CFTX * ballPred->xd[_X_];
			*bounce = TRUE;
		}
	}
	else { // it touches the floor
		for (i = 1; i <= N_CART; i++) {
			ballPred->xd[i] = 0.;
			ballPred->x[i]  = ballState.x[i];
		}
	}

}

/*
 * Condition to determine if ball hits the table
 * Useful for prediction including a rebound model
 * Useful also in KF/EKF filtering.
 */
int checkForBallTableContact(SL_Cstate state) {

	int sign;
	if (dist_to_table > 0)
		sign = 1;
	else
		sign = -1;

	double center_table = dist_to_table + sign * 0.5 * table_length;
	double dist_to_center_y = fabs(center_table - state.x[_Y_]);
	double dist_to_center_x = fabs(table_center - state.x[_X_]);
	double dist_to_table_plane = state.x[_Z_] - (floor_level - table_height + ball_radius);

	if (dist_to_center_y < table_length/2.
		&& dist_to_center_x <= table_width/2.
		&& dist_to_table_plane <= 0 && state.xd[_Z_] <= 0)
		return TRUE;
	else
		return FALSE;
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
 * Returns constant vector of val value from 1 to n
 */
void const_vec(const int n, const double val, double * vec) {

	int i;
	for (i = 0; i < n; i++) {
		vec[i] = val;
	}
}

/*
 * Estimate an initial solution for NLOPT
 * 2*dof + 1 dimensional problem
 *
 * The closer to the optimum it is the faster alg should converge
 * TODO: load values from a lookup table
 */
void guesstimate_soln(double * x) {

	x[0] = 0.45;
	x[1] = 0.41;
	x[2] = -0.08;
	x[3] = 1.65;
	x[4] = -1.29;
	x[5] = -1.05;
	x[6] = 0.18;

	x[7] = -1.61;
	x[8] = 2.91;
	x[9] = -0.24;
	x[10] = -0.56;
	x[11] = 0.94;
	x[12] = -2.45;
	x[13] = -0.31;

	// initialize first dof entries to q0
//	int i;
//	for (i = 0; i < DOF; i++) {
//		x[i] = q0[i];
//		x[i+DOF] = 0.0;
//	}
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
 * Calculate desired racket normal using the mirror law
 */
void calc_racket_normal(Vector bin, Vector bout, Vector normal) {

	vec_sub(bout, bin, normal);
	// normalize
	vec_mult_scalar(normal, 1./sqrt(vec_mult_inner(normal,normal)), normal);
}

/*
 *
 *  Computes the desired outgoing velocity of the ball after contact
 *		            to hit the goal on the opponents court
 *  Input:
 *
 *		SL_CState	hitPoint:   hitting point at contact
 *		Vector		landPoint: 	landing point on the opponents court
 *
 *	Output:
 *
 *		Vector		velOut:     desired velocity of the ball
 *
 */
void calc_ball_vel_out(SL_Cstate hitPoint, Vector landPoint, double time2land, Vector velOut) {

	double x, y, z, time2Land;
	double ynet, znet, zoc;
	double alpha;
	int sign, check, i;

	static int firsttime = TRUE;
	static double zTable;

	if (firsttime) {
		firsttime = FALSE;
		zTable = floor_level - table_height + ball_radius;
	}

	velOut[_X_] = (landPoint[1] - hitPoint.x[_X_]) / time2land;
	velOut[_Y_] = (landPoint[2] - hitPoint.x[_Y_]) / time2land;
	velOut[_Z_] = (zTable - hitPoint.x[_Z_] + 0.5 * intern_gravity * sqr(time2land)) / time2land;

	//TODO: consider the air drag case
	// hack for now
	velOut[_X_] = 1.1 * velOut[_X_];
	velOut[_Y_] = 1.1 * velOut[_Y_];
	velOut[_Z_] = 1.2 * velOut[_Z_];

}

/*
 * First order hold to interpolate linearly at time T between ball prediction matrix ballMat entries
 *
 * TODO: NO Checking for extrapolation! Only checking for NaN value.
 *
 */
void first_order_hold(double *ballPred, double *racketVel, double *racketNormal, double T) {

	int i;

	if (isnan(T)) {
		printf("Warning: T value is nan! Setting ballPred to initial value\n");
		for (i = 1; i <= CART; i++) {
			ballPred[i-1] = ballMat[0][i];
		}
		return;
	}

	int N = (int) (T/TSTEP);
	double Tdiff = T - N*TSTEP;
	static int iter;
	int Nmax = (int) TPRED/TSTEP;

	//printf("T = %f\t", T);
	//printf("Iter no = %d\n", iter++);

	for (i = 1; i <= CART; i++) {
		if (N < Nmax) {
			ballPred[i-1] = ballMat[N][i] + (Tdiff/TSTEP) * (ballMat[N+1][i] - ballMat[N][i]);
			racketVel[i-1] = racketMat[N][i+CART] + (Tdiff/TSTEP) * (racketMat[N+1][i+CART] - racketMat[N][i+CART]);
			racketNormal[i-1] = racketMat[N][i+2*CART] + (Tdiff/TSTEP) * (racketMat[N+1][i+2*CART] - racketMat[N][i+2*CART]);
		}
		else {
			ballPred[i-1] = ballMat[N][i];
			racketVel[i-1] = racketMat[N][i+CART];
			racketNormal[i-1] = racketMat[N][i+2*CART];
		}
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
 * Returns a1 - a2 vector into a1, assuming both have dof = 7 length
 */
void vec_minus(double *a1, const double *a2) {

	int i;
	for (i = 0; i < DOF; i++) {
		a1[i] = a1[i] - a2[i];
	}
}

/*
 * Returns a1 + a2 vector into a1, assuming both have dof = 7 length
 */
void vec_plus(double *a1, const double *a2) {

	int i;
	for (i = 0; i < DOF; i++) {
		a1[i] = a1[i] + a2[i];
	}
}

/*
 * Returns the inner product between two vectors of size N_DOFS
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

/*
 * Return time of day as micro seconds
 */
long get_time() {
	struct timeval tv;
	if (gettimeofday(&tv, (struct timezone *)0) == 0)
		return (tv.tv_sec*1000*1000 + tv.tv_usec);  //us

	return 0.;
}
