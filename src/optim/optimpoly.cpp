/**
 * @file optimpoly.cpp
 * @brief Nonlinear optimization in C using the NLOPT library
 * @author Okan
 * @date 30/05/2016
 *
 */

#include <armadillo>
#include "constants.h"
#include "utils.h"
#include "stdlib.h"
#include "math.h"
#include "kinematics.h"
#include "optim.h"
#include <adolc/adouble.h>            // use of active doubles
#include <adolc/drivers/drivers.h>    // use of "Easy to Use" drivers
#include <adolc/taping.h>             // use of taping

// termination
static bool check_optim_result(const int res);

// optimization related methods
static double costfunc(unsigned n, const double *x, double *grad, void *my_func_data);
static void kinematics_eq_constr(unsigned m, double *result, unsigned n,
		                  const double *x, double *grad, void *f_data);
static void first_order_hold(const optim_des* racketdata, const double T, double racket_pos[NCART],
		               double racket_vel[NCART], double racket_n[NCART]);

// autodiff versions
static adouble costfunc(unsigned n, const adouble *x, double *grad, void *my_func_params);
static void kinematics_eq_constr(unsigned m, adouble *result, unsigned n,
		                  const adouble *x, double *grad, void *my_function_data);

/**
 * @brief Update the initial state of optimization to PLAYER's current joint states.
 * @param qact
 */
void Optim::update_init_state(const joint & qact) {
	for (int i = 0; i < NDOF; i++) {
		q0[i] = qact.q(i);
		q0dot[i] = qact.qd(i);
	}
}

/**
 * @brief Tells the player optimization thread is still BUSY.
 *
 * If the (detached) thread is still running then table tennis player does not
 * update/launch new trajectories.
 * @return
 */
bool Optim::check_running() {
	return running;
}

/**
 * @brief If the optimization was successful notify the Player class
 *
 * If the optimization was successful, update is turned ON and the table tennis
 * player can launch/update the polynomial trajectories.
 * @return
 */
bool Optim::check_update() {
	return update;
}

/**
 * @brief If the robot starts moving the optimization is notified via this function
 *
 * If the robot is moving, this means last optimization was feasible, hence
 * we can initialize the new optimization from the last optimized parameter values.
 * @param flag_move
 */
void Optim::set_moving(bool flag_move) {
	moving = flag_move;
}

/**
 * @brief Detach the optimization thread.
 *
 * If the optimization is performed on SL or REAL_ROBOT then the optimization
 * thread should be detached.
 * @param flag_detach
 */
void Optim::set_detach(bool flag_detach) {
	detach = flag_detach;
}

/**
 * @brief Print verbose optimization output (detailed optimization results are printed)
 * @param flag_verbose
 */
void Optim::set_verbose(bool flag_verbose) {
	verbose = flag_verbose;
}

/**
 * @brief If optimization succeeded, update polynomial parameters p
 *
 * If the optimizers finished running and terminated successfully,
 * then generates the striking and returning polynomial parameters
 * from qf, qfdot and T, and given the actual joint states qact, qactdot
 *
 */
bool Optim::get_params(const joint & qact, spline_params & p) {

	bool flag = false;
	if (update && !running) {
		vec7 qf_, qfdot_, qrest_;
		for (int i = 0; i < NDOF; i++) {
			qf_(i) = qf[i];
			qfdot_(i) = qfdot[i];
			qrest_(i) = qrest[i];
		}
		vec7 qnow = qact.q;
		vec7 qdnow = qact.qd;
		p.a.col(0) = 2.0 * (qnow - qf_) / pow(T,3) + (qfdot_ + qdnow) / pow(T,2);
		p.a.col(1) = 3.0 * (qf_ - qnow) / pow(T,2) - (qfdot_ + 2.0*qdnow) / T;
		p.a.col(2) = qdnow;
		p.a.col(3) = qnow;
		//cout << "A = \n" << p.a << endl;
		p.b.col(0) = 2.0 * (qf_ - qrest_) / pow(time2return,3) + (qfdot_) / pow(time2return,2);
		p.b.col(1) = 3.0 * (qrest_ - qf_) / pow(time2return,2) - (2.0*qfdot_) / time2return;
		p.b.col(2) = qfdot_;
		p.b.col(3) = qf_;
		p.time2hit = T;
		//cout << "B = \n" << p.b << endl;
		flag = true;
		update = false;
	}
	return flag;
}

void Optim::set_des_params(optim_des *params_) {
	delete param_des;
	param_des = params_;
}

void Optim::run() {

	std::thread t = std::thread(&Optim::optim,this);
	if (detach) {
		t.detach();
	}
	else {
		t.join();
	}
}


void Optim::optim() {

	update = false;
	running = true;
	double x[OPTIM_DIM];

	if (moving) {
		init_last_soln(x);
	}
	else {
		init_rest_soln(x);
	}

	double init_time = get_time();
	double past_time = 0.0;
	double minf; // the minimum objective value, upon return //
	int res; // error code

	if ((res = nlopt_optimize(opt, x, &minf)) < 0) {
		if (verbose)
			printf("NLOPT failed with exit code %d!\n", res);
	    past_time = (get_time() - init_time)/1e3;
	}
	else {
		past_time = (get_time() - init_time)/1e3;
		if (verbose) {
			printf("NLOPT success with exit code %d!\n", res);
			printf("NLOPT took %f ms\n", past_time);
			printf("Found minimum at f = %0.10g\n", minf);
		}
	    if (test_soln(x) < 1e-2)
	    	finalize_soln(x,past_time);
	}
	if (verbose)
		check_optim_result(res);
	running = false;
}

FocusedOptim::FocusedOptim(double qrest_[NDOF], double lb_[NDOF], double ub_[NDOF], bool grad) {

	GRAD_BASED_OPT = grad;
	double tol_eq[EQ_CONSTR_DIM];
	double tol_ineq[INEQ_CONSTR_DIM];
	const_vec(EQ_CONSTR_DIM,1e-2,tol_eq);
	const_vec(INEQ_CONSTR_DIM,1e-3,tol_ineq);
	// set tolerances equal to second argument

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
	nlopt_set_xtol_rel(opt, 1e-2);
	nlopt_set_lower_bounds(opt, lb_);
	nlopt_set_upper_bounds(opt, ub_);
	nlopt_set_min_objective(opt, costfunc, this);
	//nlopt_add_inequality_mconstraint(opt, INEQ_CONSTR_DIM, joint_limits_ineq_constr, this, tol_ineq);
	nlopt_add_equality_mconstraint(opt, EQ_CONSTR_DIM, kinematics_eq_constr, this, tol_eq);

	for (int i = 0; i < NDOF; i++) {
		qrest[i] = qrest_[i];
		ub[i] = ub_[i];
		lb[i] = lb_[i];
	}
}

/*
 * Release auto-diff tape (jacobian double-pointer)
 */
FocusedOptim::~FocusedOptim() {

	//nlopt_destroy(opt);
	if (GRAD_BASED_OPT) {
		for (int i = 0; i < EQ_CONSTR_DIM; i++)
		    delete[] jac[i];
		delete[] jac;
	}
}

void FocusedOptim::init_last_soln(double x[]) const {

	// initialize first dof entries to q0
	for (int i = 0; i < NDOF; i++) {
		x[i] = qf[i];
		x[i+NDOF] = qfdot[i];
	}
	x[2*NDOF] = T;
	//cout << "Initialization from T = " << T << endl;

}

void FocusedOptim::init_rest_soln(double x[]) const {

	// initialize first dof entries to q0
	for (int i = 0; i < NDOF; i++) {
		x[i] = qrest[i];
		x[i+NDOF] = 0.0;
	}
	x[2*NDOF] = 0.5;
}

void FocusedOptim::finalize_soln(const double x[], double time_elapsed) {

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
}

double FocusedOptim::test_soln(const double x[]) const {

	// give info on constraint violation
	double *grad = 0;
	static double kin_violation[EQ_CONSTR_DIM];
	static double lim_violation[INEQ_CONSTR_DIM]; // joint limit violations on strike and return
	kinematics_eq_constr(EQ_CONSTR_DIM, kin_violation,
			             OPTIM_DIM, x, grad, (void*)this);
	joint_limits_ineq_constr(INEQ_CONSTR_DIM, lim_violation,
			                 OPTIM_DIM, x, grad, (void*)this);
	double cost = costfunc(OPTIM_DIM, x, grad, (void*)this);

	if (verbose) {
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
 * Initialize ADOLC operations
 * Objective has easily computed gradients so no need to auto-diff
 */
void FocusedOptim::generate_tape() {

	int n = OPTIM_DIM;
	int m = EQ_CONSTR_DIM;

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

/*
 * Calculates the cost function for table tennis trajectory generation optimization
 * to find spline (3rd order strike+return) polynomials
 */
static double costfunc(unsigned n, const double *x, double *grad, void *my_func_params) {

	double a1[NDOF];
	double a2[NDOF];
	double T = x[2*NDOF];

	FocusedOptim *opt = (FocusedOptim*) my_func_params;
	double *q0 = opt->q0;
	double *q0dot = opt->q0dot;

	if (grad) {
		gradient(1,n,x,grad);
	}

	// calculate the polynomial coeffs which are used in the cost calculation
	calc_strike_poly_coeff(q0,q0dot,x,a1,a2);

	return T * (3*T*T*inner_prod(NDOF,a1,a1) +
			3*T*inner_prod(NDOF,a1,a2) + inner_prod(NDOF,a2,a2));
}

/*
 * Calculates the cost function for table tennis trajectory generation optimization
 * to find spline (3rd order strike+return) polynomials
 */
static adouble costfunc(unsigned n, const adouble *x, double *grad, void *my_func_params) {

	adouble a1[NDOF];
	adouble a2[NDOF];
	adouble T = x[2*NDOF];

	FocusedOptim *opt = (FocusedOptim*) my_func_params;
	double *q0 = opt->q0;
	double *q0dot = opt->q0dot;

	// calculate the polynomial coeffs which are used in the cost calculation
	calc_strike_poly_coeff(q0,q0dot,x,a1,a2);

	return T * (3*T*T*inner_prod(NDOF,a1,a1) +
			3*T*inner_prod(NDOF,a1,a2) + inner_prod(NDOF,a2,a2));
}

/*
 * This is the constraint that makes sure we hit the ball
 */
static void kinematics_eq_constr(unsigned m, double *result, unsigned n,
		                  const double *x, double *grad, void *my_function_data) {

	static double racket_des_pos[NCART];
	static double racket_des_vel[NCART];
	static double racket_des_normal[NCART];
	static double racket_des_pos_diff[NCART];
	static double racket_des_vel_diff[NCART];
	static double racket_des_normal_diff[NCART];
	static double racket_des_pos_diff_pre[NCART];
	static double racket_des_vel_diff_pre[NCART];
	static double racket_des_normal_diff_pre[NCART];
	static double pos[NCART];
	static double qfdot[NDOF];
	static double vel[NCART];
	static double normal[NCART];
	static double qf[NDOF];
	static const double g = -9.8;
	double T = x[2*NDOF];

	FocusedOptim *opt = (FocusedOptim*) my_function_data;
	optim_des* racket_data(opt->param_des);

	// interpolate at time T to get the desired racket parameters
	first_order_hold(racket_data,T,racket_des_pos,racket_des_vel,racket_des_normal);

	if (grad) {
        jacobian(2,m,n,x,opt->jac);
        int idx = 0;
        for(int i = 0; i < m; i++) {
        	for(int j = 0; j < n; j++) {
        		grad[idx++] = opt->jac[i][j];
        		//printf("%.2f\t",vhp->jac[i][j]);
        	}
        	//printf("\n");
        }
        // change the T dependencies!
        // last column all entries
        double h = 1e-6;
        first_order_hold(racket_data,T+h,racket_des_pos_diff,racket_des_vel_diff,racket_des_normal_diff);
        first_order_hold(racket_data,T-h,racket_des_pos_diff_pre,racket_des_vel_diff_pre,racket_des_normal_diff_pre);
        grad[(0+X)*n + n-1] = -((racket_des_pos_diff[X] - racket_des_pos_diff_pre[X]) / (2*h));
        grad[(0+Y)*n + n-1] = -((racket_des_pos_diff[Y] - racket_des_pos_diff_pre[Y]) / (2*h));
        grad[(0+Z)*n + n-1] = -((racket_des_pos_diff[Z] - racket_des_pos_diff_pre[Z]) / (2*h));
        grad[(3+X)*n + n-1] = -((racket_des_vel_diff[X] - racket_des_vel_diff_pre[X]) / (2*h));
        grad[(3+Y)*n + n-1] = -((racket_des_vel_diff[Y] - racket_des_vel_diff_pre[Y]) / (2*h));
        grad[(3+Z)*n + n-1] = -((racket_des_vel_diff[Z] - racket_des_vel_diff_pre[Z]) / (2*h));
        grad[(6+X)*n + n-1] = -((racket_des_normal_diff[X] - racket_des_normal_diff_pre[X]) / (2*h));
        grad[(6+Y)*n + n-1] = -((racket_des_normal_diff[Y] - racket_des_normal_diff_pre[Y]) / (2*h));
        grad[(6+Z)*n + n-1] = -((racket_des_normal_diff[Z] - racket_des_normal_diff_pre[Z]) / (2*h));
        idx = 0;
        for(int i = 0; i < m; i++) {
        	for(int j = 0; j < n; j++) {
        		printf("%.2f\t",grad[idx++]);
        	}
        	printf("\n");
        }
        printf("\n");
	}

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
 * This is the constraint that makes sure we hit the ball
 */
static void kinematics_eq_constr(unsigned m, adouble *result, unsigned n,
		                  const adouble *x, double *grad, void *my_function_data) {

	static double racket_des_pos[NCART];
	static double racket_des_vel[NCART];
	static double racket_des_normal[NCART];
	static adouble pos[NCART];
	static adouble qfdot[NDOF];
	static adouble vel[NCART];
	static adouble normal[NCART];
	static adouble qf[NDOF];
	adouble T = x[2*NDOF];

	FocusedOptim *opt = (FocusedOptim*) my_function_data;
	optim_des* racket_data(opt->param_des);

	// interpolate at time T to get the desired racket parameters
	//first_order_hold(racket_data,T,racket_des_pos,racket_des_vel,racket_des_normal);

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
static void first_order_hold(const optim_des * data, const double T, double racket_pos[NCART],
		               double racket_vel[NCART], double racket_n[NCART]) {

	double deltat = data->dt;
	if (std::isnan(T)) {
		printf("Warning: T value is nan!\n");

		for(int i = 0; i < NCART; i++) {
			racket_pos[i] = data->racket_pos(i,0);
			racket_vel[i] = data->racket_vel(i,0);
			racket_n[i] = data->racket_normal(i,0);
		}
	}
	else {
		int N = (int) (T/deltat);
		double Tdiff = T - N*deltat;
		int Nmax = data->Nmax;

		for (int i = 0; i < NCART; i++) {
			if (N < Nmax - 1) {
				racket_pos[i] = data->racket_pos(i,N) +
						(Tdiff/deltat) * (data->racket_pos(i,N+1) - data->racket_pos(i,N));
				racket_vel[i] = data->racket_vel(i,N) +
						(Tdiff/deltat) * (data->racket_vel(i,N+1) - data->racket_vel(i,N));
				racket_n[i] = data->racket_normal(i,N) +
						(Tdiff/deltat) * (data->racket_normal(i,N+1) - data->racket_normal(i,N));
			}
			else {
				racket_pos[i] = data->racket_pos(i,Nmax-1);
				racket_vel[i] = data->racket_vel(i,Nmax-1);
				racket_n[i] = data->racket_normal(i,Nmax-1);
			}
		}
	}
}

/*
 * This is the inequality constraint that makes sure we never exceed the
 * joint limits during the striking and returning motion
 *
 */
void joint_limits_ineq_constr(unsigned m, double *result,
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

	FocusedOptim *opt = (FocusedOptim*) my_func_params;
	double *q0 = opt->q0;
	double *q0dot = opt->q0dot;
	double *qrest = opt->qrest;
	double *ub = opt->ub;
	double *lb = opt->lb;
	double Tret = opt->time2return;

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
void calc_strike_poly_coeff(const double *q0, const double *q0dot, const double *x,
		                    double *a1, double *a2) {

	double T = x[2*NDOF];

	for (int i = 0; i < NDOF; i++) {
		a1[i] = (2/pow(T,3))*(q0[i]-x[i]) + (1/(T*T))*(q0dot[i] + x[i+NDOF]);
		a2[i] = (3/(T*T))*(x[i]-q0[i]) - (1/T)*(x[i+NDOF] + 2*q0dot[i]);
	}

	return;
}

/*
 * Calculate the polynomial coefficients from the optimized variables qf,qfdot,T
 * p(t) = a1*t^3 + a2*t^2 + a3*t + a4
 */
void calc_strike_poly_coeff(const double *q0, const double *q0dot, const adouble *x,
		                    adouble *a1, adouble *a2) {

	adouble T = x[2*NDOF];

	for (int i = 0; i < NDOF; i++) {
		a1[i] = (2/(T*T*T))*(q0[i]-x[i]) + (1/(T*T))*(q0dot[i] + x[i+NDOF]);
		a2[i] = (3/(T*T))*(x[i]-q0[i]) - (1/T)*(x[i+NDOF] + 2*q0dot[i]);
	}

	return;
}

/*
 * Calculate the returning polynomial coefficients from the optimized variables qf,qfdot
 * and time to return constant T
 * p(t) = a1*t^3 + a2*t^2 + a3*t + a4
 */
void calc_return_poly_coeff(const double *q0, const double *q0dot,
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
void calc_strike_extrema_cand(const double *a1, const double *a2, const double T,
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
void calc_return_extrema_cand(const double *a1, const double *a2,
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
