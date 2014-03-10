/*
 * location_constraint.h
 *
 *  Created on: Dec 29, 2009
 *      Author: frehse
 */

#ifndef LOCATION_CONSTRAINT_H_
#define LOCATION_CONSTRAINT_H_

#include "core/hybrid_automata/location_id.h"

namespace hybrid_automata {

/** A location constraint <loc,pos> implicitly attributes a set of
 * locations. loc is location_id and pos is boolean.
 *
 * If the constraint is positive (pos=true), the set consists of only loc.
 * If it is negative (pos=false), the set consists of all locations except loc.
 */
class location_constraint {
public:
	explicit location_constraint(location_id id, bool s = true) :
		my_loc_id(id), my_sign(s) {
	}
	;

	const location_id& get_id() const {
		return my_loc_id;
	}
	;

	void set_id(location_id id) {
		my_loc_id = id;
	}
	;

	bool get_sign() const {
		return my_sign;
	}
	;

	void set_sign(bool s) {
		my_sign = s;
	}
	;

	/** Returns true if c2 implies *this (is stronger than *this),
	 *  i.e. whether *this is implied by c2 (is weaker than c2). */
	bool contains(const location_constraint& c2) const {
		if (my_sign) {
			// a positive constraint needs to be matched by the same loc
			return c2.get_sign() && my_loc_id == c2.get_id();
		} else {
			// *this is negative
			if (c2.get_sign() == false) {
				// it needs to match exactly
				return my_loc_id == c2.get_id();
			} else {
				// the loc of c2 needs to differ from the forbidden one
				return my_loc_id != c2.get_id();
			}
		}
	}
	;

	/** Returns true if *this and c2 contradict (can have no common
	 * locations). */
	bool is_disjoint_from(const location_constraint& c2) const {
		if ((my_sign && c2.my_sign) || (!my_sign && !c2.my_sign))
			return my_loc_id != c2.my_loc_id;
		else
			return my_loc_id == c2.my_loc_id;
	}
	;

	bool operator==(const location_constraint& c2) const {
		return my_sign == c2.my_sign && my_loc_id == c2.my_loc_id;
	}
	;

	bool operator!=(const location_constraint& c2) const {
		return !operator==(c2);
	}
	;

	/** Comparison operator, order true before false */
	bool operator<(const location_constraint& c2) const {
		if (my_sign && !c2.my_sign)
			return true;
		else if (!my_sign && c2.my_sign)
			return false;
		else
			//both have the same sign
			return my_loc_id < c2.my_loc_id;
	}
	;

private:
	location_id my_loc_id;
	bool my_sign;
	location_constraint(); // Disallow the default constructor
};

}

#endif /* LOCATION_CONSTRAINT_H_ */
