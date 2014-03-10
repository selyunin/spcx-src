/*
 * sf_unary.hpp
 *
 *  Created on: Apr 21, 2010
 *      Author: frehse
 */

#ifndef SF_UNARY_HPP_
#define SF_UNARY_HPP_

#include "sf_unary.h"

#include "math/vdom/vdom_matrix_operators.h"
#include "math/vdom/vdom_vector_operators.h"
#include "math/vdom/affine_map_utility.h"
#include "math/vdom/positional_vdomain.h"

namespace continuous {
namespace support_function {

template<typename scalar_type> sf_unary<scalar_type>::sf_unary(
		const support_function_provider::const_ptr& s) :
	my_set(s) {
}
;

template<typename scalar_type> sf_unary<scalar_type>::sf_unary(
		const support_function_provider::const_ptr& s, const affine_map& M) :
	sf_set<scalar_type> (M), my_set(s) {
	// if s is a sf_unary, we can merge the maps
	if (const_ptr b = boost::dynamic_pointer_cast<const sf_unary<scalar_type> >(s)) {
		my_set = b->my_set;
		if (b->get_map()) {
			set_map(math::concatenate(M, *b->get_map()));
		}
	}
}
;

template<typename scalar_type> sf_unary<scalar_type>::~sf_unary() {
}
;

template<typename scalar_type> sf_unary<scalar_type>* sf_unary<scalar_type>::clone() const {
	support_function_provider::const_ptr new_root(my_set->clone());
	if (this->get_map())
		return new sf_unary<scalar_type> (new_root, *this->get_map());
	else
		return new sf_unary<scalar_type> (new_root);
}

template<typename scalar_type> support_function_provider::const_ptr sf_unary<scalar_type>::get_unmapped_set() const {
	return my_set;
}

template<typename scalar_type> int sf_unary<scalar_type>::get_memory() const {
	throw std::runtime_error("sf_unary : missing implementation get_memory");
	return 0;
}

template<typename scalar_type> continuous_set_predicate::ptr sf_unary<
		scalar_type>::get_predicate() const {
	throw std::runtime_error("sf_unary : missing implementation get_predicate");
	return continuous_set_predicate::ptr();
}

template<typename scalar_type>
math::tribool sf_unary<scalar_type>::is_universe() const {
	return my_set->is_universe();
}

template<typename scalar_type>
math::tribool sf_unary<scalar_type>::is_empty() const {
	return my_set->is_empty();
}

template<typename scalar_type> void sf_unary<scalar_type>::print(
		std::ostream& os) const {
	os << "support_function( " << my_set;
	if (this->get_map()) {
		os << ", mapped by " << this->get_map();
	}
	os << " )";
}

template<typename scalar_type> const variable_id_set& sf_unary<scalar_type>::get_variable_ids() const {
	if (!this->get_map() || this->get_map()->is_empty())
		return my_set->get_variable_ids();
	else
		return this->get_map()->codomain().get_variable_ids();
}
template<typename scalar_type> void sf_unary<scalar_type>::reassign_primedness(
		unsigned int, unsigned int) {
	throw std::runtime_error(
			"sf_unary : missing implementation reassign_primedness");
}
template<typename scalar_type> void sf_unary<scalar_type>::increase_primedness(
		unsigned int) {
	throw std::runtime_error(
			"sf_unary : missing implementation increase_primedness");
}
template<typename scalar_type> void sf_unary<scalar_type>::decrease_primedness(
		unsigned int) {
	throw std::runtime_error(
			"sf_unary : missing implementation decrease_primedness");
}

template<typename scalar_type> bool sf_unary<scalar_type>::computes_support_vector() const {
	return my_set->computes_support_vector();
}
;

template<typename scalar_type> void sf_unary<scalar_type>::compute_support(
		const math::vdom_vector<Rational>& l, Rational& max_value,
		math::vdom_vector<Rational>& support_vec, bool& is_empty,
		bool& is_bounded) const {
	sf_set<scalar_type>::template compute_support_mapped<Rational> (*my_set, this->get_map(), l, max_value,
			support_vec, is_empty, is_bounded);
}

template<typename scalar_type> void sf_unary<scalar_type>::compute_support(
		const math::vdom_vector<double>& l, double& max_value,
		math::vdom_vector<double>& support_vec, bool& is_empty,
		bool& is_bounded) const {
	this->template compute_support_mapped<double> (*my_set, this->get_map(), l, max_value, support_vec,
			is_empty, is_bounded);
}

}
}

#endif /* SF_UNARY_HPP_ */
