/*
 * location_constraint_set_operators.h
 *
 *  Created on: Oct 2, 2013
 *      Author: notroot
 */

#ifndef LOCATION_CONSTRAINT_SET_OPERATORS_H_
#define LOCATION_CONSTRAINT_SET_OPERATORS_H_

#include "location_constraint_set.h"

namespace hybrid_automata {

/** Return only the constraints changed between two sets
 *
 * Returns the constraints in B that aren't already implied by A. */
location_constraint_set compute_changes(const location_constraint_set& A, const location_constraint_set& B);


}



#endif /* LOCATION_CONSTRAINT_SET_OPERATORS_H_ */
