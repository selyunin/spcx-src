/*
 * location_constraint_set_operators.cpp
 *
 *  Created on: Oct 2, 2013
 *      Author: notroot
 */

#include "location_constraint_set_operators.h"

namespace hybrid_automata {

location_constraint_set compute_changes(const location_constraint_set& A,
		const location_constraint_set& B) {

	location_constraint_set C;

	// Check for every constraint in B whether it's already implied by A.
	for (location_constraint_set::const_iterator it = B.begin(); it != B.end();
			++it) {
		location_constraint_set con;
		con.add_constraint(it->first, it->second);

		if (!con.contains(A)) {
			C.add_constraint(it->first, it->second);
		}
	}

	return C;
}

}

