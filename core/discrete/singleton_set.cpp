/*
 * singleton_set.cpp
 *
 *  Created on: Sep 2, 2009
 *      Author: frehse
 */

#include "core/discrete/singleton_set.h"

#include "utility/singleton_iterator.h"

namespace discrete {

singleton_set::singleton_set() :
	my_kind(EMPTY) {
}

singleton_set::singleton_set(const object_type& obj) :
	my_kind(SINGLE), my_location(obj) {
}

singleton_set::~singleton_set() {
}

discrete_set::ptr singleton_set::get_ptr() {
	discrete_set::ptr p = boost::enable_shared_from_this<singleton_set>::shared_from_this();
	return p;
}

discrete_set::const_ptr singleton_set::get_const_ptr() const {
	discrete_set::const_ptr p = boost::enable_shared_from_this<singleton_set>::shared_from_this();
	return p;
}

singleton_set* singleton_set::create_empty() const {
	return new singleton_set();
}

singleton_set* singleton_set::clone() const {
	return new singleton_set(*this);
}

const singleton_set::const_iterator& singleton_set::begin() const {
	return get_const_begin(singleton_const_iterator<object_type> (&my_location));
}

const singleton_set::const_iterator& singleton_set::end() const {
	return get_const_end(singleton_const_iterator<object_type> (&my_location, false));
}

const singleton_set::iterator& singleton_set::begin() {
	return get_begin(singleton_iterator<object_type> (&my_location));
}

const singleton_set::iterator& singleton_set::end() {
	return get_end(singleton_iterator<object_type> (&my_location, false));
}

const singleton_set::object_type& singleton_set::get_object() const {
	return my_location;
}

bool singleton_set::is_empty() const {
	return my_kind == EMPTY;
}

bool singleton_set::is_universe() const {
	return my_kind == UNIVERSE;
}
bool singleton_set::is_disjoint_from(const singleton_set& ds) const {
	if (is_universe())
		return ds.is_empty();
	else if (is_empty())
		return true;
	else {
		if (ds.is_universe())
			return is_empty();
		else if (ds.is_empty())
			return true;
		else { // neither is empty or universe
			assert(my_kind==SINGLE);
			assert(ds.my_kind==SINGLE);
			return my_location.is_disjoint_from(ds.my_location);
		}
	}
}

bool singleton_set::contains(const object_type& loc) const {
	if (is_universe() || loc.is_empty())
		return true;
	else if (is_empty())
		return false;
	else if (loc.is_universe()) {
		return false;
	} else {
		return my_location.contains(loc);
	}
}

bool singleton_set::contains(const singleton_set& ds) const {
	if (is_universe() || ds.is_empty())
		return true;
	else if (is_empty())
		return false;
	else {
		if (ds.is_universe())
			return false;
		else { // neither is empty or universe
			assert(my_kind==SINGLE);
			assert(ds.my_kind==SINGLE);
			return my_location.contains(ds.my_location);
		}
	}
}

void singleton_set::set_empty() {
	my_kind = EMPTY;
	my_location = object_type();
}

void singleton_set::set_universe() {
	my_kind = UNIVERSE;
	my_location = object_type();
}

void singleton_set::set_single(const object_type& loc) {
	my_location = loc;
	my_kind = SINGLE;
}

void singleton_set::intersection_assign(const singleton_set& ds) {
	if (!ds.is_universe() && !is_empty()) {
		if (ds.is_empty()) {
			set_empty();
		} else if (is_universe()) {
			*this = ds;
		} else {
			// neither *this nor ds are empty or universe
			my_location.intersection_assign(ds.my_location);
			my_kind = SINGLE;
			if (my_location.is_empty())
				set_empty();
		}
	} // otherwise nothing changes
}


void singleton_set::existentially_quantify(automaton_id aut_id) {
	if (my_kind == UNIVERSE)
		;
	else if (my_kind == EMPTY)
		;
	else {
//std::cout << "loc constraint before quant over " << aut_id << " : " << my_location << std::endl;
		my_location.existentially_quantify(aut_id);
//std::cout << "loc constraint after quant over " << aut_id << " : " << my_location << std::boolalpha << ", is_empty:" <<  my_location.is_empty() << ", is_universe:" << my_location.is_universe() << std::endl;
		/* if (my_location.is_empty()) {
			set_empty();
		}
		if (my_location.is_universe()) {
			set_universe();
		} */
	}
}

void singleton_set::print(std::ostream& os) const {
	if (my_kind == UNIVERSE)
		os << "true";
	else if (my_kind == EMPTY)
		os << "false";
	else
		os << my_location;
}

void singleton_set::accept(dispatching::dispatcher<discrete_set_typelist>& d) const {
	d.dispatch(this);
}

}

