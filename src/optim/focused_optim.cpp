/**
 * @file focused_optim.cpp
 * @brief Nonlinear optimization in C/C++ using the NLOPT library
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
#include "tabletennis.h"
#include "lookup.h"

namespace optim {

/*
 * Calculates the cost function for table tennis trajectory generation optimization
 * to find spline (3rd order strike+return) polynomials
 */
static double costfunc(unsigned n,
                        const double *x,
                        double *grad,
                        void *my_func_data);

/*
 * This is the constraint that makes sure we hit the ball
 */
static void kinematics_eq_constr(unsigned m,
                                double *result,
                                unsigned n,
                                const double *x,
                                double *grad,
                                void *f_data);

/*
 * First order hold to interpolate linearly at time T
 * between racket pos,vel,normal entries
 *
 * IF T is nan, racket variables are assigned to zero-element of
 * relevant racket entries
 *
 */
static void first_order_hold(const optim_des* racketdata,
                            const double T,
                            double racket_pos[NCART],
                            double racket_vel[NCART],
                            double racket_n[NCART]);

FocusedOptim::FocusedOptim(const vec7 & qrest_,
                            double lb_[2*NDOF+1],
                            double ub_[2*NDOF+1]) {

	double tol_eq[EQ_CONSTR_DIM];
	double tol_ineq[INEQ_CONSTR_DIM];
	const_vec(EQ_CONSTR_DIM,1e-2,tol_eq);
	const_vec(INEQ_CONSTR_DIM,1e-3,tol_ineq);
	// set tolerances equal to second argument

	opt = nlopt_create(NLOPT_LN_COBYLA, OPTIM_DIM);
	nlopt_set_xtol_rel(opt, 1e-2);
	nlopt_set_lower_bounds(opt, lb_);
	nlopt_set_upper_bounds(opt, ub_);
	nlopt_set_min_objective(opt, costfunc, this);
	nlopt_add_inequality_mconstraint(opt, INEQ_CONSTR_DIM, joint_limits_ineq_constr, this, tol_ineq);
	nlopt_add_equality_mconstraint(opt, EQ_CONSTR_DIM, kinematics_eq_constr, this, tol_eq);

	for (int i = 0; i < NDOF; i++) {
		qrest[i] = qrest_(i);
	}
	for (int i = 0; i < OPTIM_DIM; i++) {
		ub[i] = ub_[i];
		lb[i] = lb_[i];
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
	static double max_acc_violation; // at hitting time
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
	max_acc_violation = calc_max_acc_violation(x,q0,q0dot);

	return fmax(fmax(max_abs_array(kin_violation,EQ_CONSTR_DIM),
			    max_array(lim_violation,INEQ_CONSTR_DIM)),
			    max_acc_violation);
}

static double costfunc(unsigned n,
                        const double *x,
                        double *grad,
                        void *my_func_params) {

	double a1[NDOF];
	double a2[NDOF];
	double T = x[2*NDOF];

	if (grad) {
		static double h = 1e-6;
		static double val_plus, val_minus;
		static double xx[2*NDOF+1];
		for (unsigned i = 0; i < n; i++)
			xx[i] = x[i];
		for (unsigned i = 0; i < n; i++) {
			xx[i] += h;
			val_plus = costfunc(n, xx, NULL, my_func_params);
			xx[i] -= 2*h;
			val_minus = costfunc(n, xx, NULL, my_func_params);
			grad[i] = (val_plus - val_minus) / (2*h);
			xx[i] += h;
		}
	}

	FocusedOptim *opt = (FocusedOptim*) my_func_params;
	double *q0 = opt->q0;
	double *q0dot = opt->q0dot;

	// calculate the polynomial coeffs which are used in the cost calculation
	calc_strike_poly_coeff(q0,q0dot,x,a1,a2);

	return T * (3*T*T*inner_prod(NDOF,a1,a1) +
			3*T*inner_prod(NDOF,a1,a2) + inner_prod(NDOF,a2,a2));
}

static void kinematics_eq_constr(unsigned m,
                                double *result,
                                unsigned n,
                                const double *x,
                                double *grad,
                                void *my_function_data) {

	static double racket_des_pos[NCART];
	static double racket_des_vel[NCART];
	static double racket_des_normal[NCART];
	static double pos[NCART];
	static double qfdot[NDOF];
	static double vel[NCART];
	static double normal[NCART];
	static double qf[NDOF];
	double T = x[2*NDOF];

	FocusedOptim *opt = (FocusedOptim*) my_function_data;
	optim_des* racket_data = opt->param_des;

	if (grad) {
		static double h = 1e-6;
		static double res_plus[EQ_CONSTR_DIM], res_minus[EQ_CONSTR_DIM];
		static double xx[2*NDOF+1];
		for (unsigned i = 0; i < n; i++)
			xx[i] = x[i];
		for (unsigned i = 0; i < n; i++) {
			xx[i] += h;
			kinematics_eq_constr(m, res_plus, n, xx, NULL, my_function_data);
			xx[i] -= 2*h;
			kinematics_eq_constr(m, res_minus, n, xx, NULL, my_function_data);
			xx[i] += h;
			for (unsigned j = 0; j < m; j++)
				grad[j*n + i] = (res_plus[j] - res_minus[j]) / (2*h);
		}
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

static void first_order_hold(const optim_des* data,
                            const double T,
                            double racket_pos[NCART],
                            double racket_vel[NCART],
                            double racket_n[NCART]) {

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

void joint_limits_ineq_constr(unsigned m,
                            double *result,
                            unsigned n,
                            const double *x,
                            double *grad,
                            void *my_func_params) {

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

	if (grad) {
		static double h = 1e-6;
		static double res_plus[INEQ_CONSTR_DIM], res_minus[INEQ_CONSTR_DIM];
		static double xx[2*NDOF+1];
		for (unsigned i = 0; i < n; i++)
			xx[i] = x[i];
		for (unsigned i = 0; i < n; i++) {
			xx[i] += h;
			joint_limits_ineq_constr(m, res_plus, n, xx, NULL, my_func_params);
			xx[i] -= 2*h;
			joint_limits_ineq_constr(m, res_minus, n, xx, NULL, my_func_params);
			xx[i] += h;
			for (unsigned j = 0; j < m; j++)
				grad[j*n + i] = (res_plus[j] - res_minus[j]) / (2*h);
		}
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

void calc_strike_poly_coeff(const double *q0,
                            const double *q0dot,
                            const double *x,
		                    double *a1,
		                    double *a2) {

	double T = x[2*NDOF];

	for (int i = 0; i < NDOF; i++) {
		a1[i] = (2/pow(T,3))*(q0[i]-x[i]) + (1/(T*T))*(q0dot[i] + x[i+NDOF]);
		a2[i] = (3/(T*T))*(x[i]-q0[i]) - (1/T)*(x[i+NDOF] + 2*q0dot[i]);
	}

	return;
}

void calc_return_poly_coeff(const double *q0,
                            const double *q0dot,
		                    const double *x,
		                    const double T,
		                    double *a1,
		                    double *a2) {

	for (int i = 0; i < NDOF; i++) {
		a1[i] = (2/pow(T,3))*(x[i]-q0[i]) + (1/(T*T))*(q0dot[i] + x[i+NDOF]);
		a2[i] = (3/(T*T))*(q0[i]-x[i]) - (1/T)*(2*x[i+NDOF] + q0dot[i]);
	}

}

void calc_strike_extrema_cand(const double *a1,
                                const double *a2,
                                const double T,
		                        const double *q0,
		                        const double *q0dot,
		                        double *joint_max_cand,
		                        double *joint_min_cand) {

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

void calc_return_extrema_cand(const double *a1,
                                const double *a2,
		                        const double *x,
		                        const double Tret,
		                        double *joint_max_cand,
		                        double *joint_min_cand) {

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

double calc_max_acc_violation(const double x[2*NDOF+1],
		                        const double q0[NDOF],
		                        const double q0dot[NDOF]) {

	double T = x[2*NDOF];
	double a1[NDOF], a2[NDOF];
	double acc_abs_max = 0.0;
	double acc_max_cand;

	calc_strike_poly_coeff(q0,q0dot,x,a1,a2); // get a1,a2 out

	for (int i = 0; i < NDOF; i++) {
		acc_max_cand = fmax(fabs(6*a1[i]*T + a2[i]),fabs(a2[i]));
		//printf("qdd_max[%d] = %f\n", i, acc_max_cand);
		if (acc_max_cand > MAX_ACC && acc_max_cand > acc_abs_max) {
			acc_abs_max = acc_max_cand;
		}
	}
	return acc_abs_max;
}

void set_bounds(double *lb, double *ub, double SLACK, double Tmax) {

    using namespace std;
    string env = getenv("HOME");
    string filename = env + "/projects/table-tennis/config/Limits.cfg";
    mat joint_limits;
    joint_limits.load(filename);
    vec7 lb_ = joint_limits.col(0);
    vec7 ub_ = joint_limits.col(1);
    //read_joint_limits(lb,ub);
    // lower bounds and upper bounds for qf are the joint limits
    for (int i = 0; i < NDOF; i++) {
        ub[i] = ub_(i) - SLACK;
        lb[i] = lb_(i) + SLACK;
        ub[i+NDOF] = MAX_VEL;
        lb[i+NDOF] = -MAX_VEL;
    }
    // constraints on final time
    ub[2*NDOF] = Tmax;
    lb[2*NDOF] = 0.01;
}

}
