/*
 * polyhedron_collection.hpp
 *
 *  Created on: Oct 18, 2010
 *      Author: frehse
 */

#ifndef POLYHEDRON_COLLECTION_HPP_
#define POLYHEDRON_COLLECTION_HPP_

#include "polyhedron_collection.h"

namespace continuous {

template<typename scalar_type>
polyhedron_collection<scalar_type>::polyhedron_collection(typename polyhedron<
		scalar_type>::ptr element) {
	my_container.push_back(element);
}

template<typename scalar_type>
polyhedron_collection<scalar_type>::~polyhedron_collection() {
}

template<typename scalar_type>
polyhedron_collection<scalar_type>* polyhedron_collection<scalar_type>::clone() const {
	polyhedron_collection<scalar_type>* new_coll = new polyhedron_collection<
			scalar_type> ();
	for (const_iterator it = begin(); it != end(); ++it) {
		typename polyhedron<scalar_type>::ptr element_clone =
				typename polyhedron<scalar_type>::ptr((*it)->clone());
		new_coll->push_back(element_clone);
	}
	return new_coll;
}

template<typename scalar_type>
polyhedron_collection<scalar_type>* polyhedron_collection<scalar_type>::create_universe() const {
	typename polyhedron<scalar_type>::ptr universe_element(
			default_element()->create_universe());
	return new polyhedron_collection<scalar_type> (universe_element);
}

template<typename scalar_type>
polyhedron_collection<scalar_type>* polyhedron_collection<scalar_type>::create_empty() const {
	typename polyhedron<scalar_type>::ptr empty_element(
			default_element()->create_empty());
	return new polyhedron_collection<scalar_type> (empty_element);
}

template<typename scalar_type>
typename polyhedron<scalar_type>::const_ptr polyhedron_collection<scalar_type>::default_element() const {
	assert(size() > 0);
	return *begin();
}

template<typename scalar_type>
typename polyhedron_collection<scalar_type>::iterator polyhedron_collection<
		scalar_type>::begin() {
	return my_container.begin();
}

template<typename scalar_type>
typename polyhedron_collection<scalar_type>::iterator polyhedron_collection<
		scalar_type>::end() {
	return my_container.end();
}

template<typename scalar_type>
typename polyhedron_collection<scalar_type>::const_iterator polyhedron_collection<
		scalar_type>::begin() const {
	return my_container.begin();
}

template<typename scalar_type>
typename polyhedron_collection<scalar_type>::const_iterator polyhedron_collection<
		scalar_type>::end() const {
	return my_container.end();
}

template<typename scalar_type>
unsigned int polyhedron_collection<scalar_type>::size() const {
	return my_container.size();
}

template<typename scalar_type>
int polyhedron_collection<scalar_type>::get_memory() const {
	int mem = 0;
	for (const_iterator it = begin(); it != end(); ++it) {
		mem += (*it)->get_memory();
	}
	return mem;
}

template<typename scalar_type>
const variable_id_set& polyhedron_collection<scalar_type>::get_variable_ids() const {
	static variable_id_set vis;
	vis = variable_id_set();
	for (const_iterator it = begin(); it != end(); ++it) {
		variable_id_set element_vars = (*it)->get_variable_ids();
		vis.insert(element_vars.begin(), element_vars.end());
	}
	return vis;
}

template<typename scalar_type>
unsigned int polyhedron_collection<scalar_type>::get_dim() const {
	return get_variable_ids().size();
}

template<typename scalar_type>
math::tribool polyhedron_collection<scalar_type>::is_empty() const {
	math::tribool result(true);
	for (const_iterator it=begin(); it != end() && math::maybe(result); ++it) {
		result = result && (*it)->is_empty();
	}
	return result;
}

template<typename scalar_type>
void polyhedron_collection<scalar_type>::embed_variables(
		const variable_id_set& id_set) {
	for (iterator it = begin(); it != end(); ++it) {
		(*it)->embed_variables(id_set);
	}
}

template<typename scalar_type>
void polyhedron_collection<scalar_type>::existentially_quantify_variables(
		const variable_id_set& id_set) {
	if (!id_set.empty())
		for (iterator it = begin(); it != end(); ++it) {
			(*it)->existentially_quantify_variables(id_set);
		}
}

template<typename scalar_type>
math::tribool polyhedron_collection<scalar_type>::element_wise_contains(
		typename polyhedron<scalar_type>::const_ptr p) const {
	math::tribool contains_res(false);
	for (const_iterator it = begin(); it != end(); ++it) {
		contains_res = contains_res || (*it)->contains(p);
	}
	return contains_res;
}

template<typename scalar_type>
math::tribool polyhedron_collection<scalar_type>::element_wise_contains(
		const polyhedron_collection<scalar_type>& polys) const {
	math::tribool contains_res(true);
	for (const_iterator it = polys.begin(); it != polys.end(); ++it) {
		contains_res = contains_res && this->element_wise_contains((*it));
	}
	return contains_res;
}

template<typename scalar_type>
void polyhedron_collection<scalar_type>::delete_if_contained_in(
		typename polyhedron<scalar_type>::const_ptr p) {
	for (iterator it = begin(); it != end();) {
		if (p->contains(*it))
			it = my_container.erase(it);
		else
			++it;
	}
}

template<typename scalar_type>
void polyhedron_collection<scalar_type>::reassign_primedness(unsigned int d,
		unsigned int p) {
	for (iterator it = begin(); it != end(); ++it) {
		(*it)->reassign_primedness(d, p);
	}
}

template<typename scalar_type>
void polyhedron_collection<scalar_type>::increase_primedness(unsigned int d) {
	for (iterator it = begin(); it != end(); ++it) {
		(*it)->increase_primedness(d);
	}
}

template<typename scalar_type>
void polyhedron_collection<scalar_type>::decrease_primedness(unsigned int d) {
	for (iterator it = begin(); it != end(); ++it) {
		(*it)->decrease_primedness(d);
	}
}

template<typename scalar_type>
void polyhedron_collection<scalar_type>::insert(
		typename polyhedron<scalar_type>::ptr p, bool redundancy_check) {
//	if (my_container.empty() || !math::maybe(p->is_empty())) {
	// GF: The emptiness check should be redundant in most applications
	{
		/* test if p is already in the elements of *this.
		 The test is done element-wise since the actual operation involves
		 the difference operator, which is excessively costly. */
		if (!redundancy_check || !math::definitely(element_wise_contains(p))) {
			/* Test if p contains any element of *this. */
			if (redundancy_check) {
				delete_if_contained_in(p);
			}
			push_back(p);
		}
	}
}

template<typename scalar_type> void polyhedron_collection<scalar_type>::insert(
		const polyhedron_collection<scalar_type>& polys) {
	for (const_iterator it = polys.begin(); it != polys.end(); ++it) {
		insert(*it);
	}
}

template<typename scalar_type>
continuous_set_predicate::ptr polyhedron_collection<scalar_type>::get_predicate() const {
	throw std::runtime_error("missing implementation get_predicate");
	return continuous_set_predicate::ptr();
}

template<typename scalar_type>
void polyhedron_collection<scalar_type>::accept(dispatching::dispatcher<
		continuous_set_typelist>& d) const {
	// The following is dangerous, since the dispatcher might be stateful.
	// Accepting it with more than one element could lead to wrong results.
	d.dispatch(this);
//	for (const_iterator it = begin(); it != end(); ++it) {
//		(*it)->accept(d);
//	}
}

template<typename scalar_type>
void polyhedron_collection<scalar_type>::print(std::ostream& os) const {
	output_format of = get_output_format();
	os << of.preamble;
	for (const_iterator it = begin(); it != end(); ++it) {
		if (it != begin())
			os << of.element_separator;
		os << *it;
	}
	os << of.epilogue;
}

template<typename scalar_type>
polyhedron_collection<scalar_type>::polyhedron_collection() {
}

template<typename scalar_type>
void polyhedron_collection<scalar_type>::push_back(typename polyhedron<
		scalar_type>::ptr p) {
	my_container.push_back(p);
}

template<typename scalar_type>
void polyhedron_collection<scalar_type>::add_constraint(
		const math::lin_constraint<scalar_type> &c, bool check_redundancy) {
	for (iterator it = begin(); it != end(); ++it) {
		(*it)->add_constraint(c,check_redundancy);
	}
}

template<typename scalar_type>
void polyhedron_collection<scalar_type>::add_constraints(
		const math::lin_constraint_system<scalar_type>& con_set, bool check_redundancy) {
	for (typename math::lin_constraint_system<scalar_type>::const_iterator i =
			con_set.begin(); i != con_set.end(); ++i) {
		add_constraint(*i,check_redundancy);
	}
	// remove empty ones
	remove_empty(begin(), end());
}

template<typename scalar_type>
void polyhedron_collection<scalar_type>::remove_redundant_constraints() {
	for (iterator it = begin(); it != end(); ++it) {
		(*it)->remove_redundant_constraints();
	}
}

template<typename scalar_type>
bool polyhedron_collection<scalar_type>::computes_support_vector() const {
	// @attention this only makes since if the elements aren't empty
	for (const_iterator it = begin(); it != end(); ++it) {
		if (!(*it)->computes_support_vector())
			return false;
	}
	return true;
}

template<typename scalar_type>
void polyhedron_collection<scalar_type>::compute_support(
		const math::vdom_vector<Rational>& l, Rational& max_value,
		math::vdom_vector<Rational>& support_vec, bool& is_empty,
		bool& is_bounded) const {

	if (begin() == end())
		throw basic_exception(
				"calling compute support for a polyhedron_collection with no elements");

	/** Get the values for the first element */
	const_iterator it = begin();
	is_empty = true;
	// repeat until a non-empty one is found or the end of the list is reached
	while (is_empty && it!= end()) {
		(*it)->compute_support(l, max_value, support_vec, is_empty, is_bounded);
		++it;
	}
	/** Combine with the values for the other elements. */
	if (!is_empty && is_bounded && it!= end()) {
		for (; it != end(); ++it) {
			Rational max_value2;
			math::vdom_vector<Rational> support_vec2;
			bool is_empty2 = false;
			bool is_bounded2 = true;
			(*it)->compute_support(l, max_value2, support_vec2, is_empty2,
					is_bounded2);
			is_bounded &= is_bounded2;
			if (is_bounded2 && !is_empty2) {
				if (math::numeric::is_LT(max_value, max_value2)) {
					max_value = max_value2;
					support_vec = support_vec2;
				}
			}
		}
	}
}

template<typename scalar_type>
void polyhedron_collection<scalar_type>::compute_support(
		const math::vdom_vector<double>& l, double& max_value,
		math::vdom_vector<double>& support_vec, bool& is_empty,
		bool& is_bounded) const {

	if (begin() == end())
		throw basic_exception(
				"calling compute support for a polyhedron_collection with no elements");

	/** Get the values for the first element */
	const_iterator it = begin();
	is_empty = true;
	// repeat until a non-empty one is found or the end of the list is reached
	while (is_empty && it!= end()) {
		(*it)->compute_support(l, max_value, support_vec, is_empty, is_bounded);
		++it;
	}
	/** Combine with the values for the other elements. */
	if (!is_empty && is_bounded && it!= end()) {
		for (; it != end(); ++it) {
			double max_value2;
			math::vdom_vector<double> support_vec2;
			bool is_empty2 = false;
			bool is_bounded2 = true;
			(*it)->compute_support(l, max_value2, support_vec2, is_empty2,
					is_bounded2);
			is_bounded &= is_bounded2;
			if (is_bounded2 && !is_empty2) {
				if (math::numeric::is_LT(max_value, max_value2)) {
					max_value = max_value2;
					support_vec = support_vec2;
				}
			}
		}
	}
}

template<typename scalar_type>
void polyhedron_collection<scalar_type>::remove_empty(const iterator& ibeg,
		const iterator& iend) {
	iterator it = ibeg;
	while (it != iend && size() > 1) {
		if ((*it)->is_empty())
			it = my_container.erase(it);
		else
			++it;
	}
}

template<typename scalar_type>
void polyhedron_collection<scalar_type>::remove_empty() {
	iterator it = begin();
	while (it != end() && size() > 1) {
		if ((*it)->is_empty())
			it = my_container.erase(it);
		else
			++it;
	}
}

template<typename scalar_type>
void polyhedron_collection<scalar_type>::swap(polyhedron_collection<scalar_type>& polys) {
	if (this != &polys) {
		my_container.swap(polys.my_container);
	}
}

}

#endif /* POLYHEDRON_COLLECTION_HPP_ */
