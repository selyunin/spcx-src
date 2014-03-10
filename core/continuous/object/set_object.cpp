/*
 * set_object.cpp
 *
 *  Created on: Mar 25, 2010
 *      Author: frehse
 */

#include "set_object.h"
#include "core/continuous/polyhedra/constr_polyhedron/constr_polyhedron_constructors.h"
#include "core/continuous/polyhedra/constr_polyhedron/constr_polyhedron_operators.h"
#include "io/common_input/predicate_parser.h"

namespace continuous {
namespace object {

set_object::set_object() : my_number_type(FLOAT) {
	my_cset = continuous_set_ptr();
}

template<typename T>
struct create_polyhedron_implementor {
	static continuous_set_ptr implement(tree::node_ptr p) {
		return construct_constr_polyhedron<T> (p);
	}
	;
};

continuous_set_ptr create_polyhedron(tree::node_ptr p,
		global_types::coefficient_type t) {
	return global_types::coefficient_type_caller<continuous_set_ptr,
			tree::node_ptr, create_polyhedron_implementor>::call(p, t);
}

set_object::set_object(std::string s, number_type t) : my_number_type(t) {
	/* Get a predicate representation */
	parse_type_chooser::set_bool(global_types::STD_BOOL);
	parse_type_chooser::set_number(global_types::coefficient_type(t));
	tree::node_ptr pred_rep = predicate_parser::parse_predicate(s);

	/* Decide what's the best way to represent the object */
	// For now: just make a constraint polyhedron
	my_cset = create_polyhedron(pred_rep, global_types::coefficient_type(t));
}

//set_object::set_object(continuous_set_ptr p) : my_cset(p) {
//	// @todo assign number type (continous_set visitor)
//}

set_object::set_object(continuous_set_ptr p, number_type t) : my_cset(p),my_number_type(t) {
}

set_object set_object::clone() const  {
	continuous_set_ptr new_cset = continuous_set::ptr(my_cset->clone());
	assert(new_cset);
	return set_object(new_cset,my_number_type);
}

set_object::set_object(const set_object& X) {
	my_cset = X.my_cset;
	my_number_type = X.my_number_type;
	assert(my_cset);
}

set_object& set_object::operator=(const set_object&X) {
	if (this != &X && my_cset != X.my_cset) {
		my_cset = X.my_cset;
		my_number_type = X.my_number_type;
		assert(my_cset);
	}
	return *this;
}

bool set_object::is_empty() const {
	return my_cset->is_empty();
}

bool set_object::is_universe() const {
	return my_cset->is_universe();
}

bool set_object::is_disjoint_from(const set_object& X) const{
	return my_cset->is_disjoint_from(X.my_cset);
}

bool set_object::contains(const set_object& X) const{
	return my_cset->contains(X.my_cset);
}

void set_object::simplify() {
	my_cset->simplify();
}

continuous_set_const_ptr set_object::get_impl() const {
	return my_cset;
}

continuous_set_ptr set_object::get_impl() {
	return my_cset;
}

set_object::number_type set_object::get_number_type() const {
	return my_number_type;
}

void set_object::set_number_type(number_type t) {
	my_number_type = t;
}

positional_vdomain set_object::get_positional_domain() const {
	/* If the set defines a positional domain, return it.
	 * Otherwise, invent one.
	 */
	if (index_to_variable_id_map_provider* p=dynamic_cast<index_to_variable_id_map_provider*>(my_cset.get())) {
		return p->domain();
	} else {
		return positional_vdomain (my_cset->get_variable_ids());
	}
}

}
}

std::ostream& operator<<(std::ostream& os, const continuous::object::set_object& p) {
	assert(p.get_impl());
	os << p.get_impl();
	return os;
}
