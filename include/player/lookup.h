/**
 * @file lookup.h
 *
 * @brief Lookup-table based functions are declared here.
 *
 *  Created on: Feb 22, 2017
 *      Author: okoc
 */

#ifndef PLAYER_INCLUDE_LOOKUP_H_
#define PLAYER_INCLUDE_LOOKUP_H_

#define LOOKUP_TABLE_SIZE 4002 //769
#define LOOKUP_COLUMN_SIZE 2*NDOF + 1 + 2*NCART // ball state and optimization parameters (6 + 15)

static const std::string LOOKUP_TABLE_NAME = "LookupTable-16-May-2016";

using namespace arma;

void load_lookup_table(mat & lookup);
void lookup_random_entry(vec & coparams, vec & params);
void knn(const mat & lookupt, const vec & ballstate, const int k, vec & params);

#endif /* PLAYER_INCLUDE_LOOKUP_H_ */
