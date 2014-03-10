/*
 * explicit_transition.cpp
 *
 *  Created on: Aug 27, 2009
 *      Author: frehse
 */

#include "core/hybrid_automata/explicit_transition.h"

namespace hybrid_automata {

explicit_transition::explicit_transition(location_id source_loc, label_id l,
		jump_constraints jcons, location_id target_loc) :
	my_source_loc(source_loc), my_target_loc(target_loc), my_label(l), my_jump_constraints(jcons) {
}

explicit_transition::explicit_transition(location_id source_loc, std::string label_name,
		jump_constraints jcons, location_id target_loc) :
	my_source_loc(source_loc), my_target_loc(target_loc), my_label(
			named_label::get_or_add_label_id(label_name)), my_jump_constraints(jcons) {
}

explicit_transition::~explicit_transition() {
}

explicit_transition* explicit_transition::create(location_id source_loc, label_id l,
		jump_constraints jcons, location_id target_loc) const {
	return new explicit_transition(source_loc, l, jcons, target_loc);
}

label_id explicit_transition::get_label() const {
	return my_label;
}

const location_id& explicit_transition::get_source() const {
	return my_source_loc;
}

const location_id& explicit_transition::get_target() const {
	return my_target_loc;
}

const jump_constraints& explicit_transition::get_jump_constraints() const {
	return my_jump_constraints;
}

/*
void explicit_transition::set_label(label_id a) {
	my_label = a;
}
void explicit_transition::set_source(const location_id& l) {
	my_source_loc = l;
}

void explicit_transition::set_target(const location_id& l) {
	my_target_loc = l;
}
*/
void explicit_transition::set_jump_constraints(const jump_constraints& t) {
	my_jump_constraints = t;
}

}
