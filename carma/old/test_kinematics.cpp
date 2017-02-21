/*
 * kinematics.cpp
 *
 * Unit tests for Barrett WAM kinematics.
 *
 *  Created on: Feb 14, 2017
 *      Author: okoc
 */

#include <boost/test/unit_test.hpp>
#include <armadillo>
#include <string>
#include "kinematics.hpp"

using namespace std;
using namespace arma;

// all needed to compare with SL
inline double Sin(double x) { return sin(x); };
inline double Cos(double x) { return cos(x); };
inline double Power(double x, double y)	{return pow(x,y); };

void print_mat(const string & comm, double a[4+1][4+1]);
double** my_matrix(int nrl, int nrh, int ncl, int nch);
void sl_kin(double * state, double ** Xlink, double ** Xorigin, double ** Xaxis,
		    double *** Ahmat);

/*
 * Comparing with MATLAB the racket states calculated given q, qd
 *
 */
BOOST_AUTO_TEST_CASE( test_racket_state_calculation ) {

	cout << "Comparing racket state calculations with MATLAB..." << endl;
	vec3 x;
	vec3 xdot;
	vec3 n;

	vec7 q = ones<vec>(7);
	vec7 qdot = zeros<vec>(7);

	calc_racket_state(q,qdot,x,xdot,n);

	//cout << x << xdot << n << endl;

	vec3 x_MATLAB = {-0.8309, 0.0861, -0.6672};
	vec3 xdot_MATLAB = zeros<vec>(3);
	vec3 n_MATLAB = {0.3857, 0.0174, -0.9225};

	BOOST_TEST(approx_equal(x,x_MATLAB,"absdiff", 0.002));
	BOOST_TEST(approx_equal(xdot,xdot_MATLAB,"absdiff", 0.002));
	BOOST_TEST(approx_equal(n,n_MATLAB,"absdiff", 0.002));

}

/*
 * Testing the times it takes for kinematics computations
 * CARMA vs. SL
 */
BOOST_AUTO_TEST_CASE( compare_kin_times ) {

	wall_clock timer;
	mat::fixed<3,7> origin = zeros<mat>(3,7);
	mat::fixed<3,7> axis = zeros<mat>(3,7);
	mat::fixed<3,6> link = zeros<mat>(3,6);
	cube::fixed<4,4,7> amats = zeros<cube>(4,4,7);
	vec7 q = randu<vec>(7);

	kinematics(q,link,origin,axis,amats);
	timer.tic();
	kinematics(q,link,origin,axis,amats);
	double t_carma = timer.toc();

	typedef double** Matrix;
	#define N_LINKS 6
	#define N_DOFS 7
	#define _X_ 1
	#define _Y_ 2
	#define _Z_ 3
	#define N_ENDEFFS 1

	static double** Xorigin;
	static Matrix Xaxis;
	static Matrix Xlink;
	static Matrix Ahmat[N_LINKS+1];
	static double state[N_DOFS+1];
	static bool firsttime = true;

	if (firsttime) {
		firsttime = false;
		Xlink = my_matrix(0,N_LINKS,1,3);
		Xorigin = my_matrix(0,N_DOFS,1,3);
		Xaxis = my_matrix(0,N_DOFS,1,3);
		for (int i = 0; i <= N_LINKS; ++i) {
			Ahmat[i] = my_matrix(1,4,1,4);
		}
		for (int i = 1; i <= N_DOFS; i++) {
			state[i] = q(i-1);
		}
	}
	sl_kin(state,Xlink,Xorigin,Xaxis,Ahmat);
	timer.tic();
	sl_kin(state,Xlink,Xorigin,Xaxis,Ahmat);
	double t_sl = timer.toc();

	cout << "Kinematics comparison..." << endl;
	cout << "t_sl = " << t_sl << endl;
	cout << "t_carma = " << t_carma << endl;

}

/*
 * Testing whether kinematics works well
 */
BOOST_AUTO_TEST_CASE( test_forward_kinematics_matlab ) {

	cout << "Comparing kinematics with MATLAB..." << endl;

	//int seed = 1;
	//arma_rng::set_seed(seed);
	vec7 q = ones<vec>(7);

	mat::fixed<3,7> origin = zeros<mat>(3,7);
	mat::fixed<3,7> axis = zeros<mat>(3,7);
	mat::fixed<3,6> link = zeros<mat>(3,6);
	cube::fixed<4,4,7> amats = zeros<cube>(4,4,7);

	kinematics(q,link,origin,axis,amats);

	mat::fixed<6,3> link_MATLAB = {{0,0,-0.3460},
			{-0.3576,-0.2296,-0.6189},
			{-0.4210, -0.2253, -0.6227},
			{-0.4745, -0.2461, -0.6501},
			{-0.7223, -0.1907, -0.6270},
			{-0.8309, 0.0861, -0.6672}};

	mat::fixed<7,3> origin_MATLAB =
	{{0, 0, -0.3460},
			{0, 0, -0.3460},
			{-0.3576, -0.2296, -0.6189},
			{-0.4210, -0.2253, -0.6227},
			{-0.4745, -0.2461, -0.6501},
			{-0.7223, -0.1907, -0.6270},
			{-0.7223, -0.1907, -0.6270}};

	mat::fixed<7,3> axis_MATLAB =
	{{0, 0, -1.0000},
			{-0.5403, 0.8415, 0},
			{-0.7081, -0.4546, -0.5403},
			{0.0906, 0.7003, -0.7081},
			{-0.9719, 0.2175, 0.0906},
			{-0.1340, -0.1938, -0.9719},
			{-0.3620, 0.9225, -0.1340}};

	//cout << link << endl;
	BOOST_TEST(approx_equal(link,link_MATLAB.t(),"absdiff", 0.002));
	BOOST_TEST(approx_equal(origin,origin_MATLAB.t(),"absdiff", 0.002));
	BOOST_TEST(approx_equal(axis,axis_MATLAB.t(),"absdiff", 0.002));

}

/*
 * Compare with SL's kinematics functions
 *
 */
BOOST_AUTO_TEST_CASE( test_forward_kinematics_sl ) {

	cout << "Comparing kinematics with SL..." << endl;

typedef double** Matrix;
#define N_LINKS 6
#define N_DOFS 7
#define _X_ 1
#define _Y_ 2
#define _Z_ 3
#define N_ENDEFFS 1

	static double** Xorigin;
	static Matrix Xaxis;
	static Matrix Xlink;
	static Matrix Ahmat[N_LINKS+1];
	static double state[N_DOFS+1];
	static bool firsttime = true;

	if (firsttime) {
		firsttime = false;
		Xlink = my_matrix(0,N_LINKS,1,3);
		Xorigin = my_matrix(0,N_DOFS,1,3);
		Xaxis = my_matrix(0,N_DOFS,1,3);
		for (int i = 0; i <= N_LINKS; ++i) {
			Ahmat[i] = my_matrix(1,4,1,4);
		}
		for (int i = 1; i <= N_DOFS; i++) {
			state[i] = 1.0;
		}
	}
	sl_kin(state,Xlink,Xorigin,Xaxis,Ahmat);

	mat::fixed<N_LINKS,3> Xlink_sl;
	mat::fixed<N_DOFS,3> Xorigin_sl;
	mat::fixed<N_DOFS,3> Xaxis_sl;
	cube::fixed<4,4,N_LINKS+1> Amat_sl;

	for (int i = 0; i < N_LINKS; i++) {
		Xlink_sl(i,0) = Xlink[i+1][1];
		Xlink_sl(i,1) = Xlink[i+1][2];
		Xlink_sl(i,2) = Xlink[i+1][3];
	}
	for (int i = 0; i < N_DOFS; i++) {
		Xorigin_sl(i,0) = Xorigin[i+1][1];
		Xorigin_sl(i,1) = Xorigin[i+1][2];
		Xorigin_sl(i,2) = Xorigin[i+1][3];
		Xaxis_sl(i,0) = Xaxis[i+1][1];
		Xaxis_sl(i,1) = Xaxis[i+1][2];
		Xaxis_sl(i,2) = Xaxis[i+1][3];
		for (int j = 0; j < 4; j++) {
			for (int k = 0; k < 4; k++) {
				Amat_sl(j,k,i) = Ahmat[i][j+1][k+1];
			}
		}
	}
	//cout << Xlink_sl.t() << endl;

	vec7 q = ones<vec>(7);

	mat::fixed<3,7> origin = zeros<mat>(3,7);
	mat::fixed<3,7> axis = zeros<mat>(3,7);
	mat::fixed<3,6> link = zeros<mat>(3,6);
	cube::fixed<4,4,7> amats = zeros<cube>(4,4,7);

	kinematics(q,link,origin,axis,amats);

	//cout << link << endl;

	BOOST_TEST(approx_equal(link,Xlink_sl.t(),"absdiff", 0.002));
	BOOST_TEST(approx_equal(origin,Xorigin_sl.t(),"absdiff", 0.002));
	BOOST_TEST(approx_equal(axis,Xaxis_sl.t(),"absdiff", 0.002));
	BOOST_TEST(approx_equal(amats,Amat_sl,"absdiff", 0.002));
}

/*
 * Prints homogeneous transformation matrices
 */
void print_mat(const string & comm, double a[4+1][4+1]) {

	int i,j;

	cout << "Matrix >" << comm << endl;

	for (i=1; i<=4; ++i) {
		printf("          ");
		for (j=1; j<=4; ++j) {
			printf("% 8.4f ",a[i][j]);
		}
		printf("\n");
	}
	printf("\n");


}

void sl_kin(double state[N_DOFS+1], double ** Xlink, double ** Xorigin, double ** Xaxis,
		           double *** Ahmat) {

	typedef struct {
		double x[4];
		double a[4];
	} SL_endeff;

	static double basec[3+1] = {0.0};
	static double baseo[4+1] = {0.0};

	baseo[2] = 1.0;

	/*
	 * Copied from SL_user_common.c for convenience
	 *
	 */
	SL_endeff * eff = (SL_endeff*)calloc(2,sizeof(SL_endeff));

	// attach the racket
	eff[1].x[_Z_] = .3;


	static double  sstate1th;
	static double  cstate1th;
	static double  sstate2th;
	static double  cstate2th;
	static double  sstate3th;
	static double  cstate3th;
	static double  sstate4th;
	static double  cstate4th;
	static double  sstate5th;
	static double  cstate5th;
	static double  sstate6th;
	static double  cstate6th;
	static double  sstate7th;
	static double  cstate7th;

	static double  rseff1a1;
	static double  rceff1a1;
	static double  rseff1a2;
	static double  rceff1a2;
	static double  rseff1a3;
	static double  rceff1a3;

	static double  Hi00[4+1][4+1];
	static double  Hi01[4+1][4+1];
	static double  Hi12[4+1][4+1];
	static double  Hi23[4+1][4+1];
	static double  Hi34[4+1][4+1];
	static double  Hi45[4+1][4+1];
	static double  Hi56[4+1][4+1];
	static double  Hi67[4+1][4+1];
	static double  Hi78[4+1][4+1];

	static double  Ai01[4+1][4+1];
	static double  Ai02[4+1][4+1];
	static double  Ai03[4+1][4+1];
	static double  Ai04[4+1][4+1];
	static double  Ai05[4+1][4+1];
	static double  Ai06[4+1][4+1];
	static double  Ai07[4+1][4+1];
	static double  Ai08[4+1][4+1];

	/* Need [n_joints+1]x[3+1] matrices: Xorigin,Xmcog,Xaxis, and Xlink[nLinks+1][3+1] */

	/* sine and cosine precomputation */
	sstate1th=Sin(state[1]);
	cstate1th=Cos(state[1]);

	sstate2th=Sin(state[2]);
	cstate2th=Cos(state[2]);

	sstate3th=Sin(state[3]);
	cstate3th=Cos(state[3]);

	sstate4th=Sin(state[4]);
	cstate4th=Cos(state[4]);

	sstate5th=Sin(state[5]);
	cstate5th=Cos(state[5]);

	sstate6th=Sin(state[6]);
	cstate6th=Cos(state[6]);

	sstate7th=Sin(state[7]);
	cstate7th=Cos(state[7]);

	/* rotation matrix sine and cosine precomputation */


	rseff1a1=Sin(eff[1].a[1]);
	rceff1a1=Cos(eff[1].a[1]);

	rseff1a2=Sin(eff[1].a[2]);
	rceff1a2=Cos(eff[1].a[2]);

	rseff1a3=Sin(eff[1].a[3]);
	rceff1a3=Cos(eff[1].a[3]);



	/* inverse homogeneous rotation matrices */
	Hi00[1][1]=-1 + 2*Power(baseo[1],2) + 2*Power(baseo[2],2);
	Hi00[1][2]=2*(baseo[2]*baseo[3] - baseo[1]*baseo[4]);
	Hi00[1][3]=2*(baseo[1]*baseo[3] + baseo[2]*baseo[4]);
	Hi00[1][4]=basec[1];

	Hi00[2][1]=2*(baseo[2]*baseo[3] + baseo[1]*baseo[4]);
	Hi00[2][2]=-1 + 2*Power(baseo[1],2) + 2*Power(baseo[3],2);
	Hi00[2][3]=2*(-(baseo[1]*baseo[2]) + baseo[3]*baseo[4]);
	Hi00[2][4]=basec[2];

	Hi00[3][1]=2*(-(baseo[1]*baseo[3]) + baseo[2]*baseo[4]);
	Hi00[3][2]=2*(baseo[1]*baseo[2] + baseo[3]*baseo[4]);
	Hi00[3][3]=-1 + 2*Power(baseo[1],2) + 2*Power(baseo[4],2);
	Hi00[3][4]=basec[3];


	Hi01[1][1]=cstate1th;
	Hi01[1][2]=-sstate1th;

	Hi01[2][1]=sstate1th;
	Hi01[2][2]=cstate1th;

	Hi01[3][4]=ZSFE;


	Hi12[2][1]=sstate2th;
	Hi12[2][2]=cstate2th;

	Hi12[3][1]=cstate2th;
	Hi12[3][2]=-sstate2th;


	Hi23[1][4]=ZHR;

	Hi23[2][1]=sstate3th;
	Hi23[2][2]=cstate3th;

	Hi23[3][1]=-cstate3th;
	Hi23[3][2]=sstate3th;


	Hi34[2][1]=sstate4th;
	Hi34[2][2]=cstate4th;
	Hi34[2][4]=YEB;

	Hi34[3][1]=cstate4th;
	Hi34[3][2]=-sstate4th;
	Hi34[3][4]=ZEB;


	Hi45[1][4]=ZWR;

	Hi45[2][1]=sstate5th;
	Hi45[2][2]=cstate5th;
	Hi45[2][4]=YWR;

	Hi45[3][1]=-cstate5th;
	Hi45[3][2]=sstate5th;


	Hi56[2][1]=sstate6th;
	Hi56[2][2]=cstate6th;

	Hi56[3][1]=cstate6th;
	Hi56[3][2]=-sstate6th;
	Hi56[3][4]=ZWFE;


	Hi67[2][1]=sstate7th;
	Hi67[2][2]=cstate7th;

	Hi67[3][1]=-cstate7th;
	Hi67[3][2]=sstate7th;


	Hi78[1][1]=rceff1a2*rceff1a3;
	Hi78[1][2]=-(rceff1a2*rseff1a3);
	Hi78[1][3]=rseff1a2;
	Hi78[1][4]=eff[1].x[1];

	Hi78[2][1]=rceff1a3*rseff1a1*rseff1a2 + rceff1a1*rseff1a3;
	Hi78[2][2]=rceff1a1*rceff1a3 - rseff1a1*rseff1a2*rseff1a3;
	Hi78[2][3]=-(rceff1a2*rseff1a1);
	Hi78[2][4]=eff[1].x[2];

	Hi78[3][1]=-(rceff1a1*rceff1a3*rseff1a2) + rseff1a1*rseff1a3;
	Hi78[3][2]=rceff1a3*rseff1a1 + rceff1a1*rseff1a2*rseff1a3;
	Hi78[3][3]=rceff1a1*rceff1a2;
	Hi78[3][4]=eff[1].x[3];


	/*print_mat("Hi00", Hi00);
	print_mat("Hi01", Hi01);
	print_mat("Hi12", Hi12);
	print_mat("Hi23", Hi23);
	print_mat("Hi34", Hi34);
	print_mat("Hi45", Hi45);
	print_mat("Hi56", Hi56);
	print_mat("Hi67", Hi67);
	print_mat("Hi78", Hi78);*/


	/* per link inverse homogeneous rotation matrices */
	Ai01[1][1]=Hi00[1][1]*Hi01[1][1] + Hi00[1][2]*Hi01[2][1];
	Ai01[1][2]=Hi00[1][1]*Hi01[1][2] + Hi00[1][2]*Hi01[2][2];
	Ai01[1][3]=Hi00[1][3];
	Ai01[1][4]=Hi00[1][4] + Hi00[1][3]*Hi01[3][4];

	Ai01[2][1]=Hi00[2][1]*Hi01[1][1] + Hi00[2][2]*Hi01[2][1];
	Ai01[2][2]=Hi00[2][1]*Hi01[1][2] + Hi00[2][2]*Hi01[2][2];
	Ai01[2][3]=Hi00[2][3];
	Ai01[2][4]=Hi00[2][4] + Hi00[2][3]*Hi01[3][4];

	Ai01[3][1]=Hi00[3][1]*Hi01[1][1] + Hi00[3][2]*Hi01[2][1];
	Ai01[3][2]=Hi00[3][1]*Hi01[1][2] + Hi00[3][2]*Hi01[2][2];
	Ai01[3][3]=Hi00[3][3];
	Ai01[3][4]=Hi00[3][4] + Hi00[3][3]*Hi01[3][4];


	Ai02[1][1]=Ai01[1][2]*Hi12[2][1] + Ai01[1][3]*Hi12[3][1];
	Ai02[1][2]=Ai01[1][2]*Hi12[2][2] + Ai01[1][3]*Hi12[3][2];
	Ai02[1][3]=-Ai01[1][1];
	Ai02[1][4]=Ai01[1][4];

	Ai02[2][1]=Ai01[2][2]*Hi12[2][1] + Ai01[2][3]*Hi12[3][1];
	Ai02[2][2]=Ai01[2][2]*Hi12[2][2] + Ai01[2][3]*Hi12[3][2];
	Ai02[2][3]=-Ai01[2][1];
	Ai02[2][4]=Ai01[2][4];

	Ai02[3][1]=Ai01[3][2]*Hi12[2][1] + Ai01[3][3]*Hi12[3][1];
	Ai02[3][2]=Ai01[3][2]*Hi12[2][2] + Ai01[3][3]*Hi12[3][2];
	Ai02[3][3]=-Ai01[3][1];
	Ai02[3][4]=Ai01[3][4];


	Ai03[1][1]=Ai02[1][2]*Hi23[2][1] + Ai02[1][3]*Hi23[3][1];
	Ai03[1][2]=Ai02[1][2]*Hi23[2][2] + Ai02[1][3]*Hi23[3][2];
	Ai03[1][3]=Ai02[1][1];
	Ai03[1][4]=Ai02[1][4] + Ai02[1][1]*Hi23[1][4];

	Ai03[2][1]=Ai02[2][2]*Hi23[2][1] + Ai02[2][3]*Hi23[3][1];
	Ai03[2][2]=Ai02[2][2]*Hi23[2][2] + Ai02[2][3]*Hi23[3][2];
	Ai03[2][3]=Ai02[2][1];
	Ai03[2][4]=Ai02[2][4] + Ai02[2][1]*Hi23[1][4];

	Ai03[3][1]=Ai02[3][2]*Hi23[2][1] + Ai02[3][3]*Hi23[3][1];
	Ai03[3][2]=Ai02[3][2]*Hi23[2][2] + Ai02[3][3]*Hi23[3][2];
	Ai03[3][3]=Ai02[3][1];
	Ai03[3][4]=Ai02[3][4] + Ai02[3][1]*Hi23[1][4];


	Ai04[1][1]=Ai03[1][2]*Hi34[2][1] + Ai03[1][3]*Hi34[3][1];
	Ai04[1][2]=Ai03[1][2]*Hi34[2][2] + Ai03[1][3]*Hi34[3][2];
	Ai04[1][3]=-Ai03[1][1];
	Ai04[1][4]=Ai03[1][4] + Ai03[1][2]*Hi34[2][4] + Ai03[1][3]*Hi34[3][4];

	Ai04[2][1]=Ai03[2][2]*Hi34[2][1] + Ai03[2][3]*Hi34[3][1];
	Ai04[2][2]=Ai03[2][2]*Hi34[2][2] + Ai03[2][3]*Hi34[3][2];
	Ai04[2][3]=-Ai03[2][1];
	Ai04[2][4]=Ai03[2][4] + Ai03[2][2]*Hi34[2][4] + Ai03[2][3]*Hi34[3][4];

	Ai04[3][1]=Ai03[3][2]*Hi34[2][1] + Ai03[3][3]*Hi34[3][1];
	Ai04[3][2]=Ai03[3][2]*Hi34[2][2] + Ai03[3][3]*Hi34[3][2];
	Ai04[3][3]=-Ai03[3][1];
	Ai04[3][4]=Ai03[3][4] + Ai03[3][2]*Hi34[2][4] + Ai03[3][3]*Hi34[3][4];


	Ai05[1][1]=Ai04[1][2]*Hi45[2][1] + Ai04[1][3]*Hi45[3][1];
	Ai05[1][2]=Ai04[1][2]*Hi45[2][2] + Ai04[1][3]*Hi45[3][2];
	Ai05[1][3]=Ai04[1][1];
	Ai05[1][4]=Ai04[1][4] + Ai04[1][1]*Hi45[1][4] + Ai04[1][2]*Hi45[2][4];

	Ai05[2][1]=Ai04[2][2]*Hi45[2][1] + Ai04[2][3]*Hi45[3][1];
	Ai05[2][2]=Ai04[2][2]*Hi45[2][2] + Ai04[2][3]*Hi45[3][2];
	Ai05[2][3]=Ai04[2][1];
	Ai05[2][4]=Ai04[2][4] + Ai04[2][1]*Hi45[1][4] + Ai04[2][2]*Hi45[2][4];

	Ai05[3][1]=Ai04[3][2]*Hi45[2][1] + Ai04[3][3]*Hi45[3][1];
	Ai05[3][2]=Ai04[3][2]*Hi45[2][2] + Ai04[3][3]*Hi45[3][2];
	Ai05[3][3]=Ai04[3][1];
	Ai05[3][4]=Ai04[3][4] + Ai04[3][1]*Hi45[1][4] + Ai04[3][2]*Hi45[2][4];


	Ai06[1][1]=Ai05[1][2]*Hi56[2][1] + Ai05[1][3]*Hi56[3][1];
	Ai06[1][2]=Ai05[1][2]*Hi56[2][2] + Ai05[1][3]*Hi56[3][2];
	Ai06[1][3]=-Ai05[1][1];
	Ai06[1][4]=Ai05[1][4] + Ai05[1][3]*Hi56[3][4];

	Ai06[2][1]=Ai05[2][2]*Hi56[2][1] + Ai05[2][3]*Hi56[3][1];
	Ai06[2][2]=Ai05[2][2]*Hi56[2][2] + Ai05[2][3]*Hi56[3][2];
	Ai06[2][3]=-Ai05[2][1];
	Ai06[2][4]=Ai05[2][4] + Ai05[2][3]*Hi56[3][4];

	Ai06[3][1]=Ai05[3][2]*Hi56[2][1] + Ai05[3][3]*Hi56[3][1];
	Ai06[3][2]=Ai05[3][2]*Hi56[2][2] + Ai05[3][3]*Hi56[3][2];
	Ai06[3][3]=-Ai05[3][1];
	Ai06[3][4]=Ai05[3][4] + Ai05[3][3]*Hi56[3][4];


	Ai07[1][1]=Ai06[1][2]*Hi67[2][1] + Ai06[1][3]*Hi67[3][1];
	Ai07[1][2]=Ai06[1][2]*Hi67[2][2] + Ai06[1][3]*Hi67[3][2];
	Ai07[1][3]=Ai06[1][1];
	Ai07[1][4]=Ai06[1][4];

	Ai07[2][1]=Ai06[2][2]*Hi67[2][1] + Ai06[2][3]*Hi67[3][1];
	Ai07[2][2]=Ai06[2][2]*Hi67[2][2] + Ai06[2][3]*Hi67[3][2];
	Ai07[2][3]=Ai06[2][1];
	Ai07[2][4]=Ai06[2][4];

	Ai07[3][1]=Ai06[3][2]*Hi67[2][1] + Ai06[3][3]*Hi67[3][1];
	Ai07[3][2]=Ai06[3][2]*Hi67[2][2] + Ai06[3][3]*Hi67[3][2];
	Ai07[3][3]=Ai06[3][1];
	Ai07[3][4]=Ai06[3][4];


	Ai08[1][1]=Ai07[1][1]*Hi78[1][1] + Ai07[1][2]*Hi78[2][1] + Ai07[1][3]*Hi78[3][1];
	Ai08[1][2]=Ai07[1][1]*Hi78[1][2] + Ai07[1][2]*Hi78[2][2] + Ai07[1][3]*Hi78[3][2];
	Ai08[1][3]=Ai07[1][1]*Hi78[1][3] + Ai07[1][2]*Hi78[2][3] + Ai07[1][3]*Hi78[3][3];
	Ai08[1][4]=Ai07[1][4] + Ai07[1][1]*Hi78[1][4] + Ai07[1][2]*Hi78[2][4] + Ai07[1][3]*Hi78[3][4];

	Ai08[2][1]=Ai07[2][1]*Hi78[1][1] + Ai07[2][2]*Hi78[2][1] + Ai07[2][3]*Hi78[3][1];
	Ai08[2][2]=Ai07[2][1]*Hi78[1][2] + Ai07[2][2]*Hi78[2][2] + Ai07[2][3]*Hi78[3][2];
	Ai08[2][3]=Ai07[2][1]*Hi78[1][3] + Ai07[2][2]*Hi78[2][3] + Ai07[2][3]*Hi78[3][3];
	Ai08[2][4]=Ai07[2][4] + Ai07[2][1]*Hi78[1][4] + Ai07[2][2]*Hi78[2][4] + Ai07[2][3]*Hi78[3][4];

	Ai08[3][1]=Ai07[3][1]*Hi78[1][1] + Ai07[3][2]*Hi78[2][1] + Ai07[3][3]*Hi78[3][1];
	Ai08[3][2]=Ai07[3][1]*Hi78[1][2] + Ai07[3][2]*Hi78[2][2] + Ai07[3][3]*Hi78[3][2];
	Ai08[3][3]=Ai07[3][1]*Hi78[1][3] + Ai07[3][2]*Hi78[2][3] + Ai07[3][3]*Hi78[3][3];
	Ai08[3][4]=Ai07[3][4] + Ai07[3][1]*Hi78[1][4] + Ai07[3][2]*Hi78[2][4] + Ai07[3][3]*Hi78[3][4];



	/* joint ID: 0 */
	Xorigin[0][1]=Hi00[1][4];
	Xorigin[0][2]=Hi00[2][4];
	Xorigin[0][3]=Hi00[3][4];

	/* link: {basec$0$$x[[1]], basec$0$$x[[2]], basec$0$$x[[3]]} */
	Xlink[0][1]=Hi00[1][4];
	Xlink[0][2]=Hi00[2][4];
	Xlink[0][3]=Hi00[3][4];

	Ahmat[0][1][1]=Hi00[1][1];
	Ahmat[0][1][2]=Hi00[1][2];
	Ahmat[0][1][3]=Hi00[1][3];
	Ahmat[0][1][4]=Hi00[1][4];

	Ahmat[0][2][1]=Hi00[2][1];
	Ahmat[0][2][2]=Hi00[2][2];
	Ahmat[0][2][3]=Hi00[2][3];
	Ahmat[0][2][4]=Hi00[2][4];

	Ahmat[0][3][1]=Hi00[3][1];
	Ahmat[0][3][2]=Hi00[3][2];
	Ahmat[0][3][3]=Hi00[3][3];
	Ahmat[0][3][4]=Hi00[3][4];

	Ahmat[0][4][4]=1;


	/* joint ID: 1 */
	Xorigin[1][1]=Ai01[1][4];
	Xorigin[1][2]=Ai01[2][4];
	Xorigin[1][3]=Ai01[3][4];

	Xaxis[1][1]=Ai01[1][3];
	Xaxis[1][2]=Ai01[2][3];
	Xaxis[1][3]=Ai01[3][3];

	/* link: {0, 0, ZSFE} */
	Xlink[1][1]=Ai01[1][4];
	Xlink[1][2]=Ai01[2][4];
	Xlink[1][3]=Ai01[3][4];

	Ahmat[1][1][1]=Ai02[1][1];
	Ahmat[1][1][2]=Ai02[1][2];
	Ahmat[1][1][3]=Ai02[1][3];
	Ahmat[1][1][4]=Ai02[1][4];

	Ahmat[1][2][1]=Ai02[2][1];
	Ahmat[1][2][2]=Ai02[2][2];
	Ahmat[1][2][3]=Ai02[2][3];
	Ahmat[1][2][4]=Ai02[2][4];

	Ahmat[1][3][1]=Ai02[3][1];
	Ahmat[1][3][2]=Ai02[3][2];
	Ahmat[1][3][3]=Ai02[3][3];
	Ahmat[1][3][4]=Ai02[3][4];

	Ahmat[1][4][4]=1;


	/* joint ID: 2 */
	Xorigin[2][1]=Ai02[1][4];
	Xorigin[2][2]=Ai02[2][4];
	Xorigin[2][3]=Ai02[3][4];

	Xaxis[2][1]=Ai02[1][3];
	Xaxis[2][2]=Ai02[2][3];
	Xaxis[2][3]=Ai02[3][3];

	/* joint ID: 3 */
	Xorigin[3][1]=Ai03[1][4];
	Xorigin[3][2]=Ai03[2][4];
	Xorigin[3][3]=Ai03[3][4];

	Xaxis[3][1]=Ai03[1][3];
	Xaxis[3][2]=Ai03[2][3];
	Xaxis[3][3]=Ai03[3][3];

	/* link: {ZHR, 0, 0} */
	Xlink[2][1]=Ai03[1][4];
	Xlink[2][2]=Ai03[2][4];
	Xlink[2][3]=Ai03[3][4];

	Ahmat[2][1][1]=Ai03[1][1];
	Ahmat[2][1][2]=Ai03[1][2];
	Ahmat[2][1][3]=Ai03[1][3];
	Ahmat[2][1][4]=Ai03[1][4];

	Ahmat[2][2][1]=Ai03[2][1];
	Ahmat[2][2][2]=Ai03[2][2];
	Ahmat[2][2][3]=Ai03[2][3];
	Ahmat[2][2][4]=Ai03[2][4];

	Ahmat[2][3][1]=Ai03[3][1];
	Ahmat[2][3][2]=Ai03[3][2];
	Ahmat[2][3][3]=Ai03[3][3];
	Ahmat[2][3][4]=Ai03[3][4];

	Ahmat[2][4][4]=1;


	/* joint ID: 4 */
	Xorigin[4][1]=Ai04[1][4];
	Xorigin[4][2]=Ai04[2][4];
	Xorigin[4][3]=Ai04[3][4];

	Xaxis[4][1]=Ai04[1][3];
	Xaxis[4][2]=Ai04[2][3];
	Xaxis[4][3]=Ai04[3][3];

	/* link: {0, YEB, ZEB} */
	Xlink[3][1]=Ai04[1][4];
	Xlink[3][2]=Ai04[2][4];
	Xlink[3][3]=Ai04[3][4];

	Ahmat[3][1][1]=Ai04[1][1];
	Ahmat[3][1][2]=Ai04[1][2];
	Ahmat[3][1][3]=Ai04[1][3];
	Ahmat[3][1][4]=Ai04[1][4];

	Ahmat[3][2][1]=Ai04[2][1];
	Ahmat[3][2][2]=Ai04[2][2];
	Ahmat[3][2][3]=Ai04[2][3];
	Ahmat[3][2][4]=Ai04[2][4];

	Ahmat[3][3][1]=Ai04[3][1];
	Ahmat[3][3][2]=Ai04[3][2];
	Ahmat[3][3][3]=Ai04[3][3];
	Ahmat[3][3][4]=Ai04[3][4];

	Ahmat[3][4][4]=1;


	/* joint ID: 5 */
	Xorigin[5][1]=Ai05[1][4];
	Xorigin[5][2]=Ai05[2][4];
	Xorigin[5][3]=Ai05[3][4];

	Xaxis[5][1]=Ai05[1][3];
	Xaxis[5][2]=Ai05[2][3];
	Xaxis[5][3]=Ai05[3][3];

	/* link: {ZWR, YWR, 0} */
	Xlink[4][1]=Ai05[1][4];
	Xlink[4][2]=Ai05[2][4];
	Xlink[4][3]=Ai05[3][4];

	Ahmat[4][1][1]=Ai05[1][1];
	Ahmat[4][1][2]=Ai05[1][2];
	Ahmat[4][1][3]=Ai05[1][3];
	Ahmat[4][1][4]=Ai05[1][4];

	Ahmat[4][2][1]=Ai05[2][1];
	Ahmat[4][2][2]=Ai05[2][2];
	Ahmat[4][2][3]=Ai05[2][3];
	Ahmat[4][2][4]=Ai05[2][4];

	Ahmat[4][3][1]=Ai05[3][1];
	Ahmat[4][3][2]=Ai05[3][2];
	Ahmat[4][3][3]=Ai05[3][3];
	Ahmat[4][3][4]=Ai05[3][4];

	Ahmat[4][4][4]=1;


	/* joint ID: 6 */
	Xorigin[6][1]=Ai06[1][4];
	Xorigin[6][2]=Ai06[2][4];
	Xorigin[6][3]=Ai06[3][4];

	Xaxis[6][1]=Ai06[1][3];
	Xaxis[6][2]=Ai06[2][3];
	Xaxis[6][3]=Ai06[3][3];

	/* link: {0, 0, ZWFE} */
	Xlink[5][1]=Ai06[1][4];
	Xlink[5][2]=Ai06[2][4];
	Xlink[5][3]=Ai06[3][4];

	Ahmat[5][1][1]=Ai07[1][1];
	Ahmat[5][1][2]=Ai07[1][2];
	Ahmat[5][1][3]=Ai07[1][3];
	Ahmat[5][1][4]=Ai07[1][4];

	Ahmat[5][2][1]=Ai07[2][1];
	Ahmat[5][2][2]=Ai07[2][2];
	Ahmat[5][2][3]=Ai07[2][3];
	Ahmat[5][2][4]=Ai07[2][4];

	Ahmat[5][3][1]=Ai07[3][1];
	Ahmat[5][3][2]=Ai07[3][2];
	Ahmat[5][3][3]=Ai07[3][3];
	Ahmat[5][3][4]=Ai07[3][4];

	Ahmat[5][4][4]=1;


	/* joint ID: 7 */
	Xorigin[7][1]=Ai07[1][4];
	Xorigin[7][2]=Ai07[2][4];
	Xorigin[7][3]=Ai07[3][4];

	Xaxis[7][1]=Ai07[1][3];
	Xaxis[7][2]=Ai07[2][3];
	Xaxis[7][3]=Ai07[3][3];

	/* link: {eff$1$$x[[1]], eff$1$$x[[2]], eff$1$$x[[3]]} */
	Xlink[6][1]=Ai08[1][4];
	Xlink[6][2]=Ai08[2][4];
	Xlink[6][3]=Ai08[3][4];

	Ahmat[6][1][1]=Ai08[1][1];
	Ahmat[6][1][2]=Ai08[1][2];
	Ahmat[6][1][3]=Ai08[1][3];
	Ahmat[6][1][4]=Ai08[1][4];

	Ahmat[6][2][1]=Ai08[2][1];
	Ahmat[6][2][2]=Ai08[2][2];
	Ahmat[6][2][3]=Ai08[2][3];
	Ahmat[6][2][4]=Ai08[2][4];

	Ahmat[6][3][1]=Ai08[3][1];
	Ahmat[6][3][2]=Ai08[3][2];
	Ahmat[6][3][3]=Ai08[3][3];
	Ahmat[6][3][4]=Ai08[3][4];

	Ahmat[6][4][4]=1;

}

/*
 * Allocate memory for a simple double array structure
 * in one chunk.
 *
 * Puts some info about matrix size into (0,0) and (0,1) entries
 *
 * From numerical recipes.
 */
double** my_matrix(int nrl, int nrh, int ncl, int nch) {

#define FALSE 0
#define TRUE 1
#define NR 0
#define NC 1
#define N_MAT_INFO 3

	double  *chunk;
	int      info = FALSE;

	if (nrl==1 && ncl == 1) {
		info = TRUE;
	}

	double **m = (double **) calloc((size_t) (nrh-nrl+1+info),sizeof(double*));

	if (info) {
		m[0] = (double *) calloc((size_t) N_MAT_INFO,sizeof(double));
		m[0][NR]       = nrh-nrl+1;
		m[0][NC]       = nch-ncl+1;

	}
	else {
		m -= nrl;
	}

	chunk = (double *) calloc( (size_t) (nrh-nrl+1) * (nch-ncl+1),sizeof(double));

	for(int i = nrl ; i <= nrh; i++) {
		m[i] = (double *) &(chunk[(i-nrl)*(nch-ncl+1)]);
		m[i] -= ncl;
	}
	return m;
}
