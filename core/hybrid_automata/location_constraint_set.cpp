/*
 * location_constraint_set.cpp
 *
 *  Created on: Sep 1, 2009
 *      Author: frehse
 */

#include "core/hybrid_automata/location_constraint_set.h"

#include "core/hybrid_automata/automaton_cache.h"
#include "core/hybrid_automata/hybrid_automaton.h"
#include "core/hybrid_automata/location.h"
#include "core/hybrid_automata/automaton_name_formatter.h"

namespace hybrid_automata {

location_constraint_set::location_constraint_set() {
}

location_constraint_set::location_constraint_set(const automaton_id& aut_id, const location_id& id) {
	add_constraint(aut_id, location_constraint(id));
}

location_constraint_set::~location_constraint_set() {
}

location_constraint_set::size_type location_constraint_set::size() const {
	return my_container.size();
}

location_constraint_set::const_iterator_pair location_constraint_set::get_constraints(
		const automaton_id& aut_id) const {
	return my_container.equal_range(aut_id);
}

location_constraint_set::iterator_pair location_constraint_set::get_constraints(
		const automaton_id& aut_id) {
	return my_container.equal_range(aut_id);
}

location_constraint_set::iterator location_constraint_set::begin() {
	return my_container.begin();
}

location_constraint_set::iterator location_constraint_set::end() {
	return my_container.end();
}

location_constraint_set::const_iterator location_constraint_set::begin() const {
	return my_container.begin();
}

location_constraint_set::const_iterator location_constraint_set::end() const {
	return my_container.end();
}

bool location_constraint_set::is_empty() const {
	return my_container.find(0) != my_container.end();
}

bool location_constraint_set::is_universe() const {
	return my_container.empty();
}

bool location_constraint_set::contains(const location_constraint_set& lcons) const {
	if (!lcons.is_empty()) {
		// Check for every constraint in *this whether there is a corresponding
		// constraint that is contained in it (is stronger than it).
		for (const_iterator it = begin(); it != end(); ++it) {
			const automaton_id& aut = it->first;
			const location_constraint& con = it->second;
			const_iterator_pair jtp = lcons.my_container.equal_range(aut);
			if (jtp.first == jtp.second) {
				return false; // constraint found that doesn't exist in *this
			} else {
				const_iterator jt = jtp.first;
				bool found = false;
				while (!found && jt != jtp.second) {
					if (con.contains(jt->second))
						found = true; // found a constraint that is at lea
					++jt;
				}
				if (!found)
					return false;
			}
		}
	}
	return true;
}

bool location_constraint_set::strictly_contradicts(const location_constraint_set& lcons) const {
	// Check for every constraint in *this if there is a contradicting
	// one in lcons.
	for (container_impl_type::const_iterator it = begin(); it != end(); ++it) {
		const automaton_id& aut = it->first;
		const location_constraint& con = it->second;
		const_iterator_pair jtp = lcons.my_container.equal_range(aut);
		for (container_impl_type::const_iterator jt = jtp.first; jt != jtp.second; ++jt) {
			if (jt->second.is_disjoint_from(con))
				return true; // contradicting constraint found
		}
	}
	return false; // no contradicting constraint found
}

bool location_constraint_set::is_disjoint_from(const location_constraint_set& lcons) const {
	if (!is_empty() && !lcons.is_empty()) {
		bool is_disj= strictly_contradicts(lcons) || lcons.strictly_contradicts(*this);
//std::cerr << "check " << *this << " disj " << lcons << ":" << is_disj << std::endl;
		return is_disj;
	}
	return true;
}

void location_constraint_set::set_empty() {
	if (!is_empty()) {
		my_container = container_impl_type();
		//my_container[0] = location_constraint(0);
		my_container.insert(std::make_pair(0, location_constraint(0)));
	}
}

bool location_constraint_set::map(automaton_id aut1, automaton_id aut2) {
	if (aut1==aut2)
		return false;

	location_constraint_set::iterator_pair itp = get_constraints(aut1);
	if (itp.first != itp.second) {
		// add constraints from aut1 as constraints on aut2
		for (location_constraint_set::iterator it=itp.first;it!=itp.second;++it) {
			add_constraint(aut2,it->second);
		}
		// get rid of aut1
		existentially_quantify(aut1);
		return true;
	} else {
		return false;
	}
}

void location_constraint_set::add_constraint(automaton_id aut_id, location_constraint con) {
	if (aut_id == 0)
		set_empty();
	else {
		// if there is any contradicting constraint, *this is empty
		iterator_pair itp = my_container.equal_range(aut_id);
		bool found_disj = false; // found a contradiction
		bool found_ident = false; // found identical constraint
		for (container_impl_type::const_iterator it = itp.first; !found_disj && !found_ident && it
				!= itp.second; ++it) {
			if (it->second.is_disjoint_from(con))
				found_disj = true;
			if (it->second == con)
				found_ident = true;
		}
		if (found_disj) {
			// the result is empty
			set_empty();
		} else if (!found_ident) {
			if (con.get_sign()) {
				// a positive constraint replaces all negative constraints
				// since it's stronger (there can't be another positive
				// since it would be contradicting)
				my_container.erase(itp.first, itp.second);
			}
			my_container.insert(std::make_pair(aut_id, con));
		}
	}
}
void location_constraint_set::existentially_quantify(automaton_id aut_id) {
/* wrong:
  	iterator_pair itp = my_container.equal_range(aut_id);
	for (container_impl_type::iterator it = itp.first;it!= itp.second; ++it)
		my_container.erase(it);
*/
	size_type res = my_container.erase(aut_id);
	// res is the number of elements removed;
}

void location_constraint_set::intersection_assign(const location_constraint_set& lcons) {
	if (!is_empty())
		for (container_impl_type::const_iterator it = lcons.begin(); it != lcons.end(); ++it) {
			add_constraint(it->first, it->second);
		}
}




automaton_id_set location_constraint_set::get_automata() const {
	automaton_id_set aset;
	if (!is_empty()) {
		for (container_impl_type::const_iterator it = begin(); it != end(); ++it) {
			aset.insert(it->first);
		}
	}
	return aset;
}

bool location_constraint_set::is_complete(const automaton_id_set& aset) const {
	if (is_empty())
		return false;
	for (automaton_id_set::const_iterator it = aset.begin(); it != aset.end(); ++it) {
		const_iterator_pair itp = my_container.equal_range(*it);
		container_impl_type::const_iterator jt = my_container.find(*it);
		// not true if unconstrained or negative constraint
		// recall that if there is one negative constraint then they are
		// negative, so it suffices to check the first
		if (itp.first == itp.second || itp.first->second.get_sign() == false) {
			return false;
		}
	}
	return true;
}

bool location_constraint_set::operator==(const location_constraint_set& lcs) const {
	return my_container == lcs.my_container;
}

bool location_constraint_set::operator!=(const location_constraint_set& lcs) const {
	return my_container != lcs.my_container;
}

bool location_constraint_set::operator<(const location_constraint_set& lcs) const {
	return my_container < lcs.my_container;
}

void location_constraint_set::print(std::ostream& os) const {
	if (!is_empty()) {
		if (begin() == end()) {
			os << "true";
		} else {
			for (container_impl_type::const_iterator it = begin(); it != end(); ++it) {
				if (it != begin())
					os << " & ";
				const hybrid_automaton* aut =
						hybrid_automaton_cache::get_automaton(it->first).get();
				os << "loc(";
				if (aut)
					os << automaton_name_formatter(*aut);
				else
					os << it->first;
				os << ")";
				if (it->second.get_sign())
					os << "==";
				else
					os << "!=";
				if (aut)
					os << aut->get_location(it->second.get_id())->get_name();
				else
					os << it->second.get_id();
			}
		}
	} else
		os << "false";
}

}

