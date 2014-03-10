/*
 * support_function_provider_utility.h
 *
 *  Created on: May 18, 2010
 *      Author: frehse
 */

#ifndef SUPPORT_FUNCTION_PROVIDER_UTILITY_H_
#define SUPPORT_FUNCTION_PROVIDER_UTILITY_H_

#include "support_function_provider.h"
#include "math/numeric/container_comp.h"
#include "math/vdom/lin_constraint.h"
#include "core/continuous/polyhedra/hyperbox/finite_hyperbox_utility.h"
#include "core/continuous/support_function/sf_base/sf_unary.h"


namespace continuous {

/**
 * Computes \f$ sup(||x||) over x \in *this\f$, where norm is defined as infinity norm.
 * The support_function_provider must also be an index_to_variable_id_map_provider.
 * Returns true if *this is bounded and false otherwise, in which case the norm is infinity.
 */
template <typename scalar_type> bool get_max_infinity_norm(
		const support_function_provider& prov, scalar_type& mu) {
	mu=scalar_type(0);
	scalar_type x;
	math::vdom_vector<scalar_type> v;
	bool is_empty=false;
	bool is_bounded=true;
	/* Positive and negative direction */
	variable_id_set vars=prov.get_variable_ids();
	for (int i=0; i<2; ++i) {
		for (variable_id_set::const_iterator it=vars.begin(); it!=vars.end()
				&& !is_empty && is_bounded; ++it) {
			math::vdom_vector<scalar_type> l;
			l.set_coeff_with_id(*it, scalar_type(-1+2*i)); // first time with -1, second time with +1
			prov.compute_support(l, x, v, is_empty, is_bounded);
			//std::cout << l << ":" << x << " at " << v << ", empty:" << is_empty << ", bounded:" << is_bounded << std::endl;
			if (!is_empty && abs(x)>mu) {
				mu=abs(x);
			}
		}
	}
	return is_bounded;
}
;

/** Returns the outer polyhedral approximation for a given set of directions.
 */
template<typename scalar_type>
constr_polyhedron<scalar_type> compute_outer_poly(const support_function_provider& s, const std::set<math::vdom_vector<scalar_type>, math::numeric::lex_comp_less<scalar_type,
		math::vdom_vector> >& dirs) {
	constr_polyhedron<scalar_type> res;
	typedef typename std::set<math::vdom_vector<scalar_type>, math::numeric::lex_comp_less<scalar_type,
	math::vdom_vector> > vector_set;
	for (typename vector_set::const_iterator it = dirs.begin(); it
			!= dirs.end(); ++it) {
		double max_value;
		math::vdom_vector<double> support_vec;
		bool is_empty;
		bool is_bounded;
		math::vdom_vector<double> l = it->template convert_to<double> ();
		s.compute_support(l, max_value, support_vec, is_empty, is_bounded);
		if (is_empty) {
			res.add_constraint(
					math::lin_constraint<scalar_type>::zero_dim_false());
			return res;
		} else {
			if (is_bounded) {
				math::lin_constraint<scalar_type> con(*it, convert_element<
						scalar_type> (-max_value), LE);
				res.add_constraint(con);
			}
		}
	}
	return res;
}
;

/**
 *
 * @param G Linear constraint
 * @param S Support function provider object
 *
 * @return True if the linear constraint contains the support function provider object,
 *         false otherwise.
 */

template<typename T>
bool contains(const math::lin_constraint<T>& G, const support_function_provider::const_ptr& S){

	math::numeric::approx_comparator<T> my_comp;
	math::lin_constraint<T> g_canonic = G;

	if (!G.is_canonic()){
		g_canonic = G.get_canonical();
	}
	comparison_operator my_sign = g_canonic.get_sign();

	bool is_empty, is_max_bounded, is_min_bounded;
	T max_val_pos, max_val_neg, G_inh_term;

	math::vdom_vector<T> sv; // support vector (with iimap of S)
	math::vdom_vector<T> le = g_canonic.get_normal();

	G_inh_term = -g_canonic.get_canonic_inh_coeff();

	if(my_sign == LT){
		S->compute_support(le, max_val_pos, sv, is_empty, is_max_bounded);
		if(my_comp.is_definitely_strictly_smaller(max_val_pos, G_inh_term))
			return true;
		else
			return false;
	}
	else if(my_sign == LE){
		S->compute_support(le, max_val_pos, sv, is_empty, is_max_bounded);
		if(my_comp.is_definitely_strictly_smaller(max_val_pos, G_inh_term) || my_comp.is_maybe_equal(max_val_pos, G_inh_term))
			return true;
		else{
//			std::cout << "G_inh_term:" << G_inh_term << std::endl;
//			std::cout << "max_val_pos:" << max_val_pos << std::endl;
//			std::cout << "max_val_neg:" << max_val_neg << std::endl;
			return false;
		}
	}
	else{
		S->compute_support(le, max_val_pos, sv, is_empty, is_max_bounded);
		S->compute_support(-le, max_val_neg, sv, is_empty, is_min_bounded);
		max_val_neg = T(-1) * max_val_neg;
		if(my_comp.is_maybe_equal(max_val_pos, G_inh_term) && my_comp.is_maybe_equal(max_val_neg, G_inh_term)){
			return true;
		}

		else
			return false;
	}
}

/** Returns the center of the bounding box of U.
 *
 * Throws if U is unbounded.
 * */
template<typename scalar_type>
math::vdom_vector<scalar_type> compute_center_of_bounding_box(
		const support_function_provider& U) {
	math::vdom_vector<scalar_type> c;

	finite_hyperbox<scalar_type> U_box;
	try {
		U_box = finite_bounding_box<scalar_type>(U);
		c = U_box.get_c_dom();
	} catch (std::exception& e) {
		std::stringstream ss;
		logger::copyfmt_to(ss);
		throw basic_exception(
				"Can't compute center of the set of values u in transform x'== Ax+ b + u. Did you forget to bound the input variables?\n",
				e);
	}
	return c;
}

/** Checks if the nondeterministic input set U is a point set.
 *
 *
 *  Returns the centered input set in U_centered and the center
 *  coordinates in b.  This implies that if the input set is a point set then b is updated to be the point.
 *  If U is null or empty then b is set to 0.
 *
 *  The domain of b should have all the variables of U.
 * */
template<typename scalar_type>
bool is_input_set_point(const support_function_provider::const_ptr& U,
		support_function_provider::const_ptr& U_centered, math::vdom_vector<scalar_type>& b){


	if (!U || math::definitely(U->is_empty())) {
		// Here are no inputs, so U = {b} and centered U = {0}
		b = math::vdom_vector<scalar_type>(b.domain(), scalar_type(0));
		U_centered = U;
		return true;
	} else {
		math::vdom_vector<scalar_type> U_center_tmp = compute_center_of_bounding_box<scalar_type>(*U);
		U_center_tmp.reorder(b.domain());

		// obtain centered U
		if (!math::numeric::is_MEQ(U_center_tmp.get_vector(), scalar_type(0))) {
			// If nonzero center then shift.
			// In order not to cascade two sf_unary shifts,
			// do the shift on the original.
			math::affine_map<scalar_type> translate_map(b.domain(), - U_center_tmp);

			U_centered = support_function_provider::const_ptr(
					new support_function::sf_unary<scalar_type>(U, translate_map));

		} else {
			U_centered = U;
		}

		b = U_center_tmp;

		scalar_type mu_centered;
		get_max_infinity_norm(*U_centered, mu_centered);
		if (math::numeric::is_MEQ(mu_centered, scalar_type(0))) {
			return true;
		}
		else
			return false;
	}
}
/** Checks if the nondeterministic input set U is a point set.
 *
 *  Returns the centered input set in U_centered and the center
 *  coordinates in b.  This implies that if the input set is a point set then b is updated to be the point.
 *  If U is null or empty then b is set to 0.
 *
 *  The domain of b should have all the variables of U.
 * */
template <typename scalar_type>
bool is_input_set_point(const support_function_provider::const_ptr& U,
		math::vdom_vector<scalar_type>& b) {
	support_function_provider::const_ptr U_centered;
	return is_input_set_point(U, U_centered, b);
}

/** Compute the Euclidean distance of a set to a set of contraints
 *
 * If the distance is negative, one of the constraints is violated.
 */
template<typename scalar_type>
scalar_type distance_to_constraints(const support_function_provider& X,
		const math::lin_constraint_system<scalar_type>& cons) {
	scalar_type mindist;
	bool first = true;
	for (typename math::lin_constraint_system<scalar_type>::const_iterator cit =
			cons.begin(); cit != cons.end(); ++cit) {
		if (cit->is_equality())
			throw std::runtime_error("not implemented");
		math::vdom_vector<double> d, v;
		d = cit->get_normal().template convert_to<double>();
		bool is_empty, is_bounded;
		double s;
		X.compute_support(d, s, v, is_empty, is_bounded);
		double distance = (-cit->get_canonic_inh_coeff() - s)
				/ sqrt(scalar_product(d, d));
		scalar_type sdist(distance);
		if (first || sdist < mindist) {
			first = false;
			mindist = sdist;
		}
	}
	return mindist;
}

}
#endif /* SUPPORT_FUNCTION_PROVIDER_UTILITY_H_ */
