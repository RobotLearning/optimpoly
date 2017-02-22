/*
 * lookup.cpp
 *
 * All lookup related functions are located here.
 *
 *  Created on: Feb 22, 2017
 *      Author: okoc
 */

#include <armadillo>
#include "constants.h"
#include "lookup.h"

using namespace arma;
using namespace std;

/*
 * Load the lookup table of coparameters and main parameters of interest
 */
void load_lookup_table(mat & lookup) {

	lookup = zeros<mat>(LOOKUP_TABLE_SIZE * LOOKUP_COLUMN_SIZE, 1);
	string env = getenv("HOME");
	string filename = env + "/robolab/barrett/saveData/" +
			          LOOKUP_TABLE_NAME + ".txt";

	lookup.load(filename);
	lookup.reshape(LOOKUP_TABLE_SIZE,LOOKUP_COLUMN_SIZE);
	//cout << lookup(span(0,5),span(0,5)) << endl;
}

/*
 * Lookup a random entry with optimization coparameters (b0,v0) and optimization
 * main parameters (qf,qfdot,T)
 */
void lookup_random_entry(vec & coparams, vec & params) {

	mat lookup;
	load_lookup_table(lookup);
	int entry = as_scalar(randi<vec>(1,distr_param(0,LOOKUP_TABLE_SIZE-1)));
	vec lookup_state = lookup.row(entry).t();
	coparams = lookup_state(span(X,DZ));
	params = lookup_state(span(DZ+1,LOOKUP_COLUMN_SIZE-1));

}

/*
 *
 * K-nearest-neighbours method for looking up optimization values
 *
 * Find the closest ball states in the lookup table
 * and lookup the corresponding qf,qfdot,T values
 * and average them
 *
 * INPUT:
 *
 * 	Datapoint:     ball position and velocity estimates
 * 	Val      :     optimization parameters to be loaded
 *
 */
void knn(const mat & lookupt, const vec & datapoint, vec & val, int k) {

	// find the closest entry
	static int datalength = datapoint.n_cols;
	static bool firsttime = true;
	static vec dots;
	static mat A;

	if (firsttime) {
		firsttime = false;
		A = lookupt.cols(span(X,DZ));
		for (unsigned i = 0; i < A.n_rows; i++) {
			dots(i) = dot(A.row(i), A.row(i));
		}
	}

	uvec idx = sort_index(dots - 2*A.t()*datapoint, "descend");
	vec fullvec = zeros<vec>(lookupt.n_cols);
	for (int i = 0; i < k; i++) {
		fullvec += lookupt.row(idx(i));
	}
	val = fullvec(span(datalength,lookupt.n_cols-1))/k;
}


