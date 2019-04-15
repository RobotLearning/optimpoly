/**
 * @file lookup.h
 *
 * @brief Lookup-table based functions are declared here.
 *
 *  Created on: Feb 22, 2017
 *      Author: okoc
 */
#pragma once
#include "constants.h"
#include <string>
#include <armadillo>

const int LOOKUP_COLUMN_SIZE =
    2 * const_tt::NDOF + 1 +
    2 * const_tt::NCART; // ball state and optimization parameters (6 + 15)

static const std::string LOOKUP_TABLE_NAME = "lookup_March_2016";

using arma::mat;
using arma::vec;

namespace player {

/**
 * @brief Load the lookup table of coparameters and main parameters of interest.
 *
 * Lookup table is located in the saveData/ folder of SL (robolab/barrett/).
 *
 * @param lookup Matrix used to load lookup table
 */
void load_lookup_table(mat &lookup);

/**
 * @brief Lookup a random entry with optimization coparameters (b0,v0) and
 * optimization main parameters (qf,qfdot,T).
 *
 * Random entry is uniformly distributed between 1 and lookup table size.
 * Coparameters are 6-dimensional and main parameters are 15-dimensional.
 *
 * @param coparams 6-dimensional vector with stored ball positions and
 * velocities.
 * @param params 15-dimensional vector with (qf,qfdot,T) store optimization
 * params.
 */
void lookup_random_entry(vec &coparams, vec &params);

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
void knn(const mat &lookupt, const vec &ballstate, const int k, vec &params);

} // namespace player
