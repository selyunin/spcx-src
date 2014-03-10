#ifndef POLYHEDRON_OPERATORS_H_
#define POLYHEDRON_OPERATORS_H_

#include "polyhedron.h"

namespace continuous {

/** Computes whether the continuous set t1 contains the continuous set t2.
 * @todo what about strict inequalities? */
template<typename scalar_type> inline math::tribool compute_poly_containment(
		const polyhedron<scalar_type>& t1, const support_function_provider& t2) {
	typename math::lin_constraint_system<scalar_type>::const_ptr
			my_poly = t1.get_constraints();
	scalar_type max_val;
	typename math::vdom_vector<scalar_type> support_vec;
	bool is_empty, is_bounded;

//	std::cout << "testing containment of " << t1 << " and " << t2 << std::endl;

	if (t2.is_empty()) {
		return true; // Any non-empty set contains an empty set.
	} else {
		if (t1.is_empty())
			return false; // An empty set contains nothing

		if (t1.is_universe())
			return true; // Universe contains every set.
			//		if (t2.is_universe())
			//			return false; // Universe can be contained only in universe and t1 is not universe.

		/* Containment is true if max_val <= b for all constraints */
		math::tribool contains_res(true);
		for (typename math::lin_constraint_system<scalar_type>::const_iterator
				i = my_poly->begin(); i != my_poly->end() && math::maybe(contains_res); ++i) {
			typename math::vdom_vector<scalar_type> l =
					i->get_normal();
			scalar_type b=-i->get_canonic_inh_coeff();
			// we now have a constraint l^Tx # b, where # is the sign
//std::cout << "testing " << l.get_vector() << l.get_index_to_variable_id_map() << l.get_inh_coeff() << i->get_sign() << std::flush << std::endl;
//std::cout << *i << std::endl;

			// note that l includes the inh. coeff, and the support function value does as well,
			// so we don't have to compare max_val to l.get_inh_coeff()

			using namespace math;
			using namespace math::numeric;
			// # is <, <= or ==
			// false if the max of l^Tx > b
			t2.compute_support(l, max_val, support_vec, is_empty, is_bounded);
			contains_res = contains_res && is_bounded && is_LE(max_val,b);
			if (i->is_equality() && math::maybe(contains_res)) {
				// # is >, >= or ==
				// false if the min of l^Tx = -maxval < b
				t2.compute_support(-l, max_val, support_vec, is_empty,
						is_bounded);
				contains_res = contains_res && is_bounded && is_LE(b,-max_val);
			}
		}
		return math::maybe(contains_res);
	}
}
;
/** Computes the difference of 2 polyhedra. The result is a polyhedron.
 */

template<typename scalar_type>
typename polyhedron<scalar_type>::ptr compute_poly_difference(
		const typename polyhedron<scalar_type>::const_ptr t1,
		const typename polyhedron<scalar_type>::const_ptr t2) {

	typename polyhedron<scalar_type>::ptr diff_set_ptr =
			typename polyhedron<scalar_type>::ptr(new constr_polyhedron<
					scalar_type> ());
	typename math::lin_constraint_system<scalar_type>::const_ptr
			my_poly_1 = t1->get_constraints();
	typename math::lin_constraint_system<scalar_type>::const_ptr
			my_poly_2 = t2->get_constraints();

	diff_set_ptr->add_constraints(my_poly_1); // Copy the constraints of t1 as they are the diff poly.
	//std::cout << "Diff poly after constraints 1: " << diff_set_ptr << std::endl;

	for (typename math::lin_constraint_system<scalar_type>::const_iterator
			i = my_poly_2->begin(); i != my_poly_2->end(); ++i) {
		typename math::lin_constraint<scalar_type> lc =
				math::lin_constraint<scalar_type>();
		lc = *i;
		typename math::lin_constraint<scalar_type>::sign my_sign =
				lc.get_sign();
		if (my_sign == LT) // complement the sign of the constraint.
			lc.set_sign(GE);
		else if (my_sign == LE)
			lc.set_sign(GT);
		else if (my_sign == GE)
			lc.set_sign(LT);
		else if (my_sign == GT)
			lc.set_sign(LE);
		else if (my_sign == EQ) {
			typename math::lin_constraint<scalar_type> lc_1 =
					math::lin_constraint<scalar_type>(lc.get_l(), LT);
			diff_set_ptr->add_constraint(lc_1);
			lc.set_sign(GT);
		} else {
		}
		//std::cout << "Ls before adding" << lc << std::endl;
		diff_set_ptr->add_constraint(lc);
	}
	return diff_set_ptr;
}
;

/** Create the polyhedron resulting from applying the affine map t. */
template<typename scalar_type>
polyhedron<scalar_type>* compute_image(
		const polyhedron<scalar_type>& c, const math::affine_map<
				scalar_type>& t) {
	bool singular;

	// Compute the inverse of the linear transformation matrix
	math::vdom_matrix<scalar_type> A_inverse(
			t.get_A().inverse(singular));
	if (singular) {
		throw std::runtime_error(
				"constr_polyhedron<scalar_type> - compute_affine_transform not implemented yet for singular transforms");
		return c.create_empty();
	} else {
		math::lin_constraint_system<scalar_type> new_poly(*c.get_constraints());

		for (typename math::lin_constraint_system<scalar_type>::iterator
				it = new_poly.begin(); it != new_poly.end(); ++it) {
			math::vdom_vector<scalar_type> row = it->get_l().get_vdom_vec();
			row = row * A_inverse;
			scalar_type transf_inh_coeff = it->get_inh_coeff()
					- scalar_product(row, t.get_b());
			math::lin_expression<scalar_type> transf_le =
					math::lin_expression<scalar_type>(row,
							transf_inh_coeff);
			math::lin_constraint<scalar_type> transf_con =
					math::lin_constraint<scalar_type>(transf_le,
							it->get_sign());
			*it=transf_con;
		}

		polyhedron<scalar_type>* res=c.create_universe();
		res->add_constraints(new_poly);
		return res;
	}
}


}

#endif /*POLYHEDRON_OPERATORS_H_*/
