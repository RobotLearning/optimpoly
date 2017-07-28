/**
 * @file lookup.cpp
 *
 * @brief All lookup table related functions are located here.
 *
 *  Created on: Feb 22, 2017
 *      Author: okoc
 */

#include <armadillo>
#include "lookup.h"
#include "constants.h"

using namespace arma;
using namespace std;

/**
 * @brief Load the lookup table of coparameters and main parameters of interest.
 *
 * Lookup table is located in the saveData/ folder of SL (robolab/barrett/).
 *
 * @param lookup Matrix used to load lookup table
 */
void load_lookup_table(mat & lookup) {

	string env = getenv("HOME");
	string filename = env + "/polyoptim/" + LOOKUP_TABLE_NAME + ".txt";

	lookup.load(filename);
	int row_size = lookup.n_elem / LOOKUP_COLUMN_SIZE;
	lookup.reshape(row_size,LOOKUP_COLUMN_SIZE);
	//cout << lookup(span(0,5),span(0,5)) << endl;
}

/**
 * @brief Lookup a random entry with optimization coparameters (b0,v0) and optimization
 * main parameters (qf,qfdot,T).
 *
 * Random entry is uniformly distributed between 1 and lookup table size.
 * Coparameters are 6-dimensional and main parameters are 15-dimensional.
 *
 * @param coparams 6-dimensional vector with stored ball positions and velocities.
 * @param params 15-dimensional vector with (qf,qfdot,T) store optimization params.
 */
void lookup_random_entry(vec & coparams, vec & params) {

	mat lookup;
	load_lookup_table(lookup);
	int entry = as_scalar(randi<vec>(1,distr_param(0,lookup.n_rows-1)));
	vec lookup_state = lookup.row(entry).t();
	coparams = lookup_state(span(X,DZ));
	params = lookup_state(span(DZ+1,LOOKUP_COLUMN_SIZE-1));

}

/**
 * @brief K-nearest-neighbours method for looking up optimization values.
 *
 * Find the closest ball states in the lookup table
 * and lookup the corresponding qf,qfdot,T values
 * and average them.
 *
 * @param lookupt Training set with each row = [coparams,params]
 * @param testpoint Test time co-parameters
 * @param val Test parameters to be loaded
 * @param k Value of the K-nearest-neighbor
 *
 * TODO: hasnt been checked for correctness for k > 1!
 *
 */
void knn(const mat & lookupt, const vec & testpoint, const int k, vec & val) {

	// find the closest entry
	static int coparam_length = testpoint.n_rows;
	static bool firsttime = true;
	static vec dots = zeros<vec>(lookupt.n_rows);
	static mat A = zeros<mat>(lookupt.n_rows,coparam_length);

	if (firsttime) {
		firsttime = false;
		A = lookupt.cols(span(0,coparam_length-1));
		for (unsigned i = 0; i < A.n_rows; i++) {
			dots(i) = dot(A.row(i), A.row(i));
		}
	}

	uvec idx = sort_index(dots - 2*A*testpoint, "descend");
	vec fullvec = zeros<vec>(lookupt.n_cols);
	for (int i = 0; i < k; i++) {
		fullvec += lookupt.row(idx(i)).t();
	}
	val = fullvec(span(coparam_length,lookupt.n_cols-1))/k;
}


