/*
 * optim.cpp
 *
 *  Created on: Feb 15, 2017
 *      Author: okoc
 */

#include <boost/test/unit_test.hpp>
#include <armadillo>
#include <string>
#include "constants.h"
#include "tabletennis.h"
#include "player.hpp"
#include "optim.hpp"
#include "nlopt.hpp"

using namespace std;
using namespace arma;

#define LOOKUP_TABLE_SIZE 769 //2639 // 4002
#define LOOKUP_COLUMN_SIZE 2*NDOF + 1 + 2*NCART
#define LOOKUP_TABLE_NAME "LookupTable-March-2016" //"LookupTable-April-2016"

inline void load_lookup_table(mat & lookup) {

	lookup = zeros<mat>(LOOKUP_TABLE_SIZE * LOOKUP_COLUMN_SIZE, 1);
	string env = getenv("HOME");
	string filename = env + "/robolab/barrett/saveData/" +
			          LOOKUP_TABLE_NAME + ".txt";

	lookup.load(filename);
	lookup.reshape(LOOKUP_TABLE_SIZE,LOOKUP_COLUMN_SIZE);
	//cout << lookup(span(0,5),span(0,5)) << endl;
}

inline void init_right_posture(vec7 & q0) {

	q0(0) = 1.0;
	q0(1) = -0.2;
	q0(2) = -0.1;
	q0(3) = 1.8;
	q0(4) = -1.57;
	q0(5) = 0.1;
	q0(6) = 0.3;
}

/*
 * Prepare racket desired positions, velocities and normals
 * for the optimization
 */
inline racket prepare_racket(const vec6 & ballstate,
		const vec2 & ball_land_des, const double time_land_des) {

	EKF filter = init_filter();
	mat66 P; P.eye();
	filter.set_prior(ballstate,P);
	double tpred = 1.0;
	double dt = 0.002;
	int N = (int)(tpred/dt);

	mat ballpred = filter.predict_path(dt,N);

	static double z_table = floor_level - table_height + ball_radius;

	// elementwise division
	mat balls_out_vel = zeros<mat>(3,N);
	balls_out_vel.row(X) = (ball_land_des(X) - ballpred.row(X)) / time_land_des;
	balls_out_vel.row(Y) = (ball_land_des(Y) - ballpred.row(Y)) / time_land_des;
	balls_out_vel.row(Z) = (z_table - ballpred.row(Z) +
			                0.5 * gravity * pow(time_land_des,2)) / time_land_des;

	balls_out_vel.row(X) *= 1.1;
	balls_out_vel.row(Y) *= 1.1;
	balls_out_vel.row(Z) *= 1.2;

	mat balls_in_vel = ballpred.rows(DX,DZ);
	mat normal = balls_out_vel - balls_in_vel;
	normalise(normal);

	mat racketvel = zeros<mat>(3,N);
	for (int i = 0; i < ballpred.n_cols; i++) {
		racketvel.col(i) = dot((balls_out_vel.col(i) + CRR * balls_in_vel.col(i) / (1 + CRR)),
								normal.col(i)) * normal.col(i);
	}

	racket racketdes;
	racketdes.pos = ballpred.rows(X,Z);
	racketdes.normal = normal;
	racketdes.vel = racketvel;
	return racketdes;
}

/*
 * Setup NLOPT to use derivative free COBYLA method
 */
inline PolyOptim setup_nlopt(const racket & racket_params) {

	vec7 q0;
	double time2return = 1.0;
	init_right_posture(q0);
	optim strike_params = {q0, zeros<vec>(7), 1.0, false};
	opt minimizer = opt(LN_COBYLA, OPTIM_DIM);
	PolyOptim polyopt = PolyOptim(q0,racket_params,
			                      time2return,strike_params,minimizer);
	polyopt.setup();

	return polyopt;
}

BOOST_AUTO_TEST_CASE( test_nlopt ) {

	cout << "Setting up NLOPT Optimizer..." << endl;

	double time_land_des = 0.8;
	vec2 ball_land_des;
	ball_land_des(X) = 0.0;
	ball_land_des(Y) = dist_to_table - 3*table_length/4;

	mat lookup;
	load_lookup_table(lookup);
	int entry = as_scalar(randi<vec>(1,distr_param(0,LOOKUP_TABLE_SIZE-1)));
	vec lookup_state = lookup.row(entry).t();
	vec6 ball_state = lookup_state(span(X,DZ));
	mat::fixed<15,1> joint_state = lookup_state(span(6,LOOKUP_COLUMN_SIZE-1));

	racket racketdes = prepare_racket(ball_state,ball_land_des,time_land_des);
	PolyOptim polyopt = setup_nlopt(racketdes);
	//polyopt();


}


