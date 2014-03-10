/*
 * positional_vdomain.hpp
 *
 *  Created on: May 12, 2010
 *      Author: frehse
 */

#ifndef POSITIONAL_VDOMAIN_HPP_
#define POSITIONAL_VDOMAIN_HPP_

#include "positional_vdomain.h"

#include <limits>

inline positional_vdomain::positional_vdomain(
		index_to_variable_id_map_ptr iimap) :
	my_iimap(iimap) {
}

inline positional_vdomain::positional_vdomain() :
	my_iimap(index_to_variable_id_map::empty_map()) {
}


inline positional_vdomain::positional_vdomain(const variable_id_set& vars) {
	my_iimap = index_to_variable_id_map::get_map_with_ids(vars);
}

inline
void positional_vdomain::add_variable(const variable& x) {
	my_iimap = my_iimap->index_to_variable_id_map::get_map_with_id_added(
			x.get_id());
}

inline
const variable_id_set& positional_vdomain::get_variable_ids() const {
	return my_iimap->get_ids();
}

inline
positional_vdomain::size_type positional_vdomain::size() const {
	return my_iimap->dimensions();
}

inline
positional_vdomain::size_type positional_vdomain::pos(
		const variable& x) const {
	return my_iimap->get_index(x.get_id());
}

inline
variable positional_vdomain::get_variable(const size_type& i) const {
	return variable(my_iimap->get_id(i));
}

inline
bool positional_vdomain::in_domain(const variable& x) const {
	return my_iimap->has_id(x.get_id());
}

inline
bool positional_vdomain::in_domain(const variable& x, size_type& pos) const {
	bool e;
	pos = my_iimap->check_for_index(x.get_id(), e);
	return e;
}

inline
bool positional_vdomain::contains_variables(const positional_vdomain& dom) const {
	return set_contains(get_variable_ids(),dom.get_variable_ids());
}

inline
const index_to_variable_id_map_ptr& positional_vdomain::get_index_to_variable_id_map() const {
	return my_iimap;
}

inline
void positional_vdomain::swap(positional_vdomain& D) {
	std::swap(my_iimap, D.my_iimap);
}

inline
void swap(positional_vdomain& D1, positional_vdomain& D2) {
	D1.swap(D2);
}

inline
bool positional_vdomain::operator==(const positional_vdomain& D) const {
	return my_iimap == D.my_iimap;
}

inline
bool positional_vdomain::operator!=(const positional_vdomain& D) const {
	return my_iimap != D.my_iimap;
}

inline
positional_vdomain::size_type positional_vdomain::invalid_pos() {
	return std::numeric_limits<positional_vdomain::size_type>::max();
}

#endif /* POSITIONAL_VDOMAIN_HPP_ */
