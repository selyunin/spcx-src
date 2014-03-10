/*
 * location_eq_node.cpp
 *
 *  Created on: Sep 4, 2009
 *      Author: frehse
 */

#include "core/hybrid_automata/location_eq_node.h"

#include "core/hybrid_automata/hybrid_automaton.h"
#include "core/hybrid_automata/automaton_cache.h"

namespace hybrid_automata {

location_eq_node::location_eq_node(const std::string& aut, const std::string& loc, bool equal) {
	hybrid_automaton_ptr a = hybrid_automaton_cache::get_automaton(aut);
	if (a) {
		my_automaton = a->get_id();
		my_location = a->get_location_id(loc);
		my_equal = equal;
	}
	else {
		throw std::runtime_error("Unknown automaton "+aut+".");
	}
}

location_eq_node::location_eq_node(automaton_id aut, location_id loc, bool equal) {
	my_automaton = aut;
	my_location = loc;
	my_equal = equal;
}

location_eq_node::~location_eq_node() {
}

const automaton_id&  location_eq_node::get_automaton_id() const {
	return my_automaton;
}

const location_id&  location_eq_node::get_location_id() const {
	return my_location;
}

const bool& location_eq_node::get_equal() const{
	return my_equal;
}

void location_eq_node::accept(tree::node_visitor& v) const {
	v.visit(this);
}

bool location_node_creator::locked = false;

}
