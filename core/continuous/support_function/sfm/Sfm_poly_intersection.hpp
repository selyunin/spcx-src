/*
 * Sfm_poly_intersection.cpp
 *
 *  Created on: Dec 30, 2010
 *      Author: ray
 */

#include "Sfm_poly_intersection.h"

namespace continuous {
	namespace support_function {



template<typename scalar_type> Sfm_poly_intersection::Sfm_poly_intersection(
		const typename sfm_cont_set<scalar_type>::const_ptr& s,
		const math::numeric::interval<unsigned int>& my_intv,
		const continuous::polyhedron<scalar_type>& cons_poly) :
	sf_unary<scalar_type> (s), my_interval(my_intv), my_poly(cons_poly){
}

template<typename scalar_type> Sfm_poly_intersection::Sfm_poly_intersection(
		const typename sfm_cont_set<scalar_type>::const_ptr& s,
		const math::numeric::interval<unsigned int>& my_intv,
		const continuous::polyhedron<scalar_type>& cons_poly,
		const affine_map& M) :
	sf_unary<scalar_type> (s, M), my_interval(my_intv), my_poly(cons_poly){
}
template<typename scalar_type> Sfm_poly_intersection<scalar_type>* Sfm_poly_intersection<
		scalar_type>::clone() const {
	support_function_provider::const_ptr new_root(this->my_set->clone());
	/*
	if (this->get_map())
		return new Sfm_poly_intersection<scalar_type> (new_root, my_con,
				*this->get_map());
	else
		return new Sfm_poly_intersection<scalar_type> (new_root, my_con); */
	throw std::runtime_error(
			"Sfm_poly_intersection : missing implementation clone");
}

template<typename scalar_type> int Sfm_poly_intersection<scalar_type>::get_memory() const {
	throw std::runtime_error(
			"Sfm_poly_intersection : missing implementation get_memory");
	return 0;
}

template<typename scalar_type> continuous_set_predicate::ptr Sfm_poly_intersection<
		scalar_type>::get_predicate() const {
	throw std::runtime_error(
			"Sfm_poly_intersection : missing implementation get_predicate");
	return continuous_set_predicate::ptr();
}

template<typename scalar_type> void Sfm_poly_intersection<scalar_type>::print(
		std::ostream& os) const {
	throw std::runtime_error(
			"Sfm_poly_intersection : missing implementation print");
}

template<typename scalar_type> const variable_id_set& Sfm_poly_intersection<
		scalar_type>::get_variable_ids() const {
	return sf_unary<scalar_type>::get_variable_ids();
}
template<typename scalar_type> void Sfm_poly_intersection<scalar_type>::reassign_primedness(
		unsigned int, unsigned int) {
	throw std::runtime_error(
			"Sfm_poly_intersection : missing implementation reassign_primedness");
}
template<typename scalar_type> void Sfm_poly_intersection<scalar_type>::increase_primedness(
		unsigned int) {
	throw std::runtime_error(
			"Sfm_poly_intersection : missing implementation increase_primedness");
}
template<typename scalar_type> void Sfm_poly_intersection<scalar_type>::decrease_primedness(
		unsigned int) {
	throw std::runtime_error(
			"Sfm_poly_intersection : missing implementation decrease_primedness");
}

template<typename scalar_type> void Sfm_poly_intersection<scalar_type>::compute_support(
		const math::vdom_vector<Rational>& l, Rational& max_value,
		math::vdom_vector<Rational>& support_vec, bool& is_empty,
		bool& is_bounded) const {
	throw std::runtime_error(
			"Sfm_poly_intersection : missing implementation compute_support Rational");
}

template<typename scalar_type> void Sfm_poly_intersection<scalar_type>::compute_support(
		const math::vdom_vector<double>& l, double& max_value,
		math::vdom_vector<double>& support_vec, bool& is_empty,
		bool& is_bounded) const {

		//@todo incomplete implementation, to check for emptyness and boundedness
		unsigned int samples; // to catch the # samples in the algo.
		sfm_cont_set<double>::const_ptr sfm_ptr = boost::static_pointer_cast<const sfm_cont_set<double>,const support_function_provider>(this->my_set);
		math::numeric::interval<double> min_intv = guard_intersection_support<double>(
				sfm_ptr ,this->my_interval, l, this->get_constraint(), double(0), samples );
		max_value = min_intv.upper().get_val();
		is_empty = false;
		is_bounded = true;
}

template<typename scalar_type> const continuous::polyhedron<scalar_type> & Sfm_poly_intersection<
		scalar_type>::get_poly() const {
	return my_poly;
}
template<typename scalar_type> Sfm_poly_intersection::~Sfm_poly_intersection() {

}

} // end of namespace support function
} // end of namespace continuous
