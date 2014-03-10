/*
 * sf_lin_con_intersection.hpp
 *
 *  Created on: Apr 27, 2010
 *      Author: frehse
 */

#ifndef SF_LIN_CON_INTERSECTION_HPP_
#define SF_LIN_CON_INTERSECTION_HPP_

#include "sf_lin_con_intersection.h"

#include "math/vdom/vdom_matrix_operators.h"
#include "math/vdom/vdom_vector_operators.h"
#include "math/vdom/affine_map_utility.h"
#include "math/vdom/positional_vdomain.h"
#include "core/continuous/support_function/intersection_by_optimization/sf_plane_intersection.h"

namespace continuous {
namespace support_function {

template<typename scalar_type> sf_lin_con_intersection<scalar_type>::sf_lin_con_intersection(
		const support_function_provider::const_ptr& s,
		const math::lin_constraint<scalar_type>& con,
		const std::string& minbrak_type,
		const double& intersection_error) :
	sf_unary<scalar_type> (s), my_con(con), my_minbrak_type(minbrak_type), my_intersection_error(intersection_error) {
}

template<typename scalar_type> sf_lin_con_intersection<scalar_type>::sf_lin_con_intersection(
		const support_function_provider::const_ptr& s,
		const math::lin_constraint<scalar_type>& con, const affine_map& M,
		const std::string& minbrak_type,
		const double& intersection_error) :
	sf_unary<scalar_type> (s, M), my_con(con), my_minbrak_type(minbrak_type), my_intersection_error(intersection_error) {
}

template<typename scalar_type> sf_lin_con_intersection<scalar_type>::~sf_lin_con_intersection() {
}

template<typename scalar_type> sf_lin_con_intersection<scalar_type>* sf_lin_con_intersection<
		scalar_type>::clone() const {
	support_function_provider::const_ptr new_root(this->my_set->clone());
	if (this->get_map())
		return new sf_lin_con_intersection<scalar_type> (new_root, my_con,
				*this->get_map(), my_minbrak_type, my_intersection_error);
	else
		return new sf_lin_con_intersection<scalar_type> (new_root, my_con, my_minbrak_type, my_intersection_error);
}

template<typename scalar_type> int sf_lin_con_intersection<scalar_type>::get_memory() const {
	throw std::runtime_error(
			"sf_lin_con_intersection : missing implementation get_memory");
	return 0;
}

template<typename scalar_type> continuous_set_predicate::ptr sf_lin_con_intersection<
		scalar_type>::get_predicate() const {
	throw std::runtime_error(
			"sf_lin_con_intersection : missing implementation get_predicate");
	return continuous_set_predicate::ptr();
}

template<typename scalar_type> void sf_lin_con_intersection<scalar_type>::print(
		std::ostream& os) const {
	throw std::runtime_error(
			"sf_lin_con_intersection : missing implementation print");
}

template<typename scalar_type> const variable_id_set& sf_lin_con_intersection<
		scalar_type>::get_variable_ids() const {
	return sf_unary<scalar_type>::get_variable_ids();
}
template<typename scalar_type> void sf_lin_con_intersection<scalar_type>::reassign_primedness(
		unsigned int, unsigned int) {
	throw std::runtime_error(
			"sf_lin_con_intersection : missing implementation reassign_primedness");
}
template<typename scalar_type> void sf_lin_con_intersection<scalar_type>::increase_primedness(
		unsigned int) {
	throw std::runtime_error(
			"sf_lin_con_intersection : missing implementation increase_primedness");
}
template<typename scalar_type> void sf_lin_con_intersection<scalar_type>::decrease_primedness(
		unsigned int) {
	throw std::runtime_error(
			"sf_lin_con_intersection : missing implementation decrease_primedness");
}

template<typename scalar_type> void sf_lin_con_intersection<scalar_type>::compute_support(
		const math::vdom_vector<Rational>& l, Rational& max_value,
		math::vdom_vector<Rational>& support_vec, bool& is_empty,
		bool& is_bounded) const {
	throw std::runtime_error(
			"sf_lin_con_intersection : missing implementation compute_support Rational");
/*	unsigned int count;
	if (!is_intersection_empty(this->my_set,
			this->get_constraint())) {
		// @todo check whether this is bounded
		is_empty = false;
		is_bounded = true;
		compute_support_intersection(this->my_set, l, max_value, is_empty,
			is_bounded, this->get_constraint(), my_minbrak_type, my_intersection_error,count);
		sample_count+=count; // to keep statistics on number of function samplings over the object.
	} else {
		is_empty = true;
		is_bounded = true;
		std::cout << "sf_lin_con_intersection: compute_support: intersection set empty" << std::endl << std::flush;
		sample_count = 0;
	}*/
}

template<typename scalar_type> void sf_lin_con_intersection<scalar_type>::compute_support(
		const math::vdom_vector<double>& l, double& max_value,
		math::vdom_vector<double>& support_vec, bool& is_empty,
		bool& is_bounded) const {
	if (!is_intersection_empty(this->my_set,
			this->get_constraint())) {
		// @todo check whether this is bounded
		is_empty = false;
		is_bounded = true;
		compute_support_intersection(this->my_set, l, max_value, is_empty,
			is_bounded, this->get_constraint(), my_minbrak_type, my_intersection_error);
	} else {
		is_empty = true;
		is_bounded = true;
		std::cout << "sf_lin_con_intersection: compute_support: intersection set empty" << std::endl << std::flush;
	}
}

template<typename scalar_type> void sf_lin_con_intersection<scalar_type>::bounded_compute_support(
		const math::vdom_vector<double>& l, double& max_value,
		math::vdom_vector<double>& support_vec, const unsigned int sample_bound,
		bool& is_empty, bool& is_bounded) const {
	if (!is_intersection_empty(this->my_set,
			this->get_constraint())) {
		// @todo check whether this is bounded
		is_empty = false;
		is_bounded = true;
		compute_support_intersection(this->my_set, l, max_value, is_empty,
			is_bounded, this->get_constraint(), sample_bound, my_minbrak_type, my_intersection_error);
	} else {
		is_empty = true;
		is_bounded = true;
		std::cout << "sf_lin_con_intersection: compute_support: intersection set empty" << std::endl << std::flush;
	}
}

template<typename scalar_type> const math::lin_constraint<scalar_type> & sf_lin_con_intersection<
		scalar_type>::get_constraint() const {
	return my_con;
}

}
}

#endif /* SF_LIN_CON_INTERSECTION_HPP_ */
