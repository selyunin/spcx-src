/*
 * sf_unary_template.hpp
 *
 *  Created on: apr 5, 2012
 *      Author: frehse
 */

#ifndef SF_UNARY_TEMPLATE_HPP_
#define SF_UNARY_TEMPLATE_HPP_

#include "sf_unary_template.h"

#include "math/vdom/vdom_matrix_operators.h"
#include "math/vdom/vdom_vector_operators.h"
#include "math/vdom/affine_map_utility.h"
#include "math/vdom/positional_vdomain.h"

namespace continuous {
namespace support_function {

template<typename scalar_type, class implementor> sf_unary_template<scalar_type,implementor>::sf_unary_template(
		const implementor_const_ptr& s) :
	my_set(s) {
}
;

template<typename scalar_type, class implementor> sf_unary_template<scalar_type,implementor>::sf_unary_template(
		const implementor_const_ptr& s, const affine_map& M) :
	sf_set<scalar_type> (M), my_set(s) {
}
;

template<typename scalar_type, class implementor> sf_unary_template<scalar_type,implementor>::~sf_unary_template() {
}
;

template<typename scalar_type, class implementor> sf_unary_template<scalar_type,implementor>* sf_unary_template<scalar_type,implementor>::clone() const {
	if (this->get_map())
		return new sf_unary_template<scalar_type,implementor> (my_set, *this->get_map());
	else
		return new sf_unary_template<scalar_type,implementor> (my_set);
}

template<typename scalar_type, class implementor> int sf_unary_template<scalar_type,implementor>::get_memory() const {
	throw std::runtime_error("sf_unary_template : missing implementation get_memory");
	return 0;
}

template<typename scalar_type, class implementor> continuous_set_predicate::ptr sf_unary_template<
		scalar_type,implementor>::get_predicate() const {
	throw std::runtime_error("sf_unary_template : missing implementation get_predicate");
	return continuous_set_predicate::ptr();
}

template<typename scalar_type, class implementor>
math::tribool sf_unary_template<scalar_type,implementor>::is_universe() const {
	return my_set->is_universe();
}

template<typename scalar_type, class implementor>
math::tribool sf_unary_template<scalar_type,implementor>::is_empty() const {
	return my_set->is_empty();
}

template<typename scalar_type, class implementor> void sf_unary_template<scalar_type,implementor>::print(
		std::ostream& os) const {
	os << "support_function( " << my_set;
	if (this->get_map()) {
		os << ", mapped by " << this->get_map();
	}
	os << " )";
}

template<typename scalar_type, class implementor> typename sf_unary_template<
		scalar_type, implementor>::implementor_const_ptr sf_unary_template<
		scalar_type, implementor>::get_implementor() const {
	return my_set;
}

template<typename scalar_type, class implementor> const variable_id_set& sf_unary_template<scalar_type,implementor>::get_variable_ids() const {
	if (!this->get_map() || this->get_map()->is_empty())
		return my_set->get_variable_ids();
	else
		return this->get_map()->codomain().get_variable_ids();
}
template<typename scalar_type, class implementor> void sf_unary_template<scalar_type,implementor>::reassign_primedness(
		unsigned int, unsigned int) {
	throw std::runtime_error(
			"sf_unary_template : missing implementation reassign_primedness");
}
template<typename scalar_type, class implementor> void sf_unary_template<scalar_type,implementor>::increase_primedness(
		unsigned int) {
	throw std::runtime_error(
			"sf_unary_template : missing implementation increase_primedness");
}
template<typename scalar_type, class implementor> void sf_unary_template<scalar_type,implementor>::decrease_primedness(
		unsigned int) {
	throw std::runtime_error(
			"sf_unary_template : missing implementation decrease_primedness");
}

template<typename scalar_type, class implementor> bool sf_unary_template<scalar_type,implementor>::computes_support_vector() const {
	return my_set->computes_support_vector();
}
;

template<typename scalar_type, class implementor> void sf_unary_template<scalar_type,implementor>::compute_support(
		const math::vdom_vector<Rational>& l, Rational& max_value,
		math::vdom_vector<Rational>& support_vec, bool& is_empty,
		bool& is_bounded) const {
	throw std::runtime_error("missing implementation");
}

template<typename scalar_type, class implementor> void sf_unary_template<scalar_type,implementor>::compute_support(
		const math::vdom_vector<double>& l, double& max_value,
		math::vdom_vector<double>& support_vec, bool& is_empty,
		bool& is_bounded) const {
	this->template compute_support_mapped<double> (*my_set, this->get_map(), l, max_value, support_vec,
			is_empty, is_bounded);
}

// Need to specialize implementor as well
//template<class implementor>
//template<> void sf_unary_template<double,implementor>::compute_support_mapped(
//		const implementor& a_set,
//		const affine_map_const_ptr& a_map,
//		const math::vdom_vector<double>& v, double& max_value,
//		math::vdom_vector<double>& support_vec, bool& is_empty,
//		bool& is_bounded) {
//	if (a_map) {
//		math::affine_map<fun_type> M = a_map->template convert_to<
//				fun_type> ();
//		// transform the cost function
//		math::vdom_vector<fun_type> lmapped(v * (M.get_A()));
//		a_set.compute_support(lmapped, max_value, support_vec, is_empty,
//				is_bounded);
//		max_value += scalar_product(v, M.get_b());
//		// transform the support vector if there is one
//		if (support_vec.size() > 0) {
//			support_vec = M.map(support_vec);
//		}
//	} else {
//		a_set.compute_support(v, max_value, support_vec, is_empty, is_bounded);
//	}
//}

template<typename scalar_type, class implementor>
template<typename fun_type> void sf_unary_template<scalar_type,implementor>::compute_support_mapped(
		const implementor& a_set,
		const typename sf_set<scalar_type>::affine_map_const_ptr& a_map,
		const math::vdom_vector<fun_type>& v, fun_type& max_value,
		math::vdom_vector<fun_type>& support_vec, bool& is_empty,
		bool& is_bounded) {
	if (a_map) {
		math::affine_map<fun_type> M = a_map->template convert_to<
				fun_type> ();
		// transform the cost function
		math::vdom_vector<fun_type> lmapped(v * (M.get_A()));
		a_set.compute_support(lmapped, max_value, support_vec, is_empty,
				is_bounded);
		max_value += scalar_product(v, M.get_b());
		// transform the support vector if there is one
		if (support_vec.size() > 0) {
			support_vec = M.map(support_vec);
		}
	} else {
		a_set.compute_support(v, max_value, support_vec, is_empty, is_bounded);
	}
}

}
}

#endif /* SF_UNARY_TEMPLATE_HPP_ */
