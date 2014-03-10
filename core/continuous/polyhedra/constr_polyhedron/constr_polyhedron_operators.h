/*
 * constr_polyhedron_operators.h
 *
 *  Created on: Mar 30, 2010
 *      Author: frehse
 */

#ifndef CONSTR_POLYHEDRON_OPERATORS_H_
#define CONSTR_POLYHEDRON_OPERATORS_H_

#include "math/vdom/lin_constraint.h"
#include "constr_polyhedron.h"
#include "core/continuous/support_function_provider_utility.h"

//#include "core/continuous/continuous_set_operator_implementations/compute_transformation.h"
//#include "core/continuous/continuous_set_transforms/reset_affine_transform.h"

namespace continuous {

template<typename scalar_type>
constr_polyhedron<scalar_type> apply_map(
		const constr_polyhedron<scalar_type>& c, const math::affine_map<
				scalar_type>& t) {
	bool singular;

	// Test that all variables in c are in the domain of t
	//variable_id_set vis=t.domain().get_variable_ids();
	variable_id_set vis = c.get_variable_ids();
	set_difference_assign(vis,t.domain().get_variable_ids());
	if (!vis.empty()) {
		std::stringstream ss;
	//	logger::copyfmt_to(ss);
		print_variable_id_set(ss,vis);
		std::stringstream str;
	//	logger::copyfmt_to(str);
		str << ",apply map:" << t << std::endl << "\ndomain: "<< t.domain() << "\ncodomain: "<< t.codomain();
		str << "polyhedron variables: ";
		print_variable_id_set(str,c.get_variable_ids());
		str << "polyhedron constraints: " << c;
		throw basic_exception(
						"Trying to map a polyhedron with a map that doesn't have variables "+ss.str()+" of the polyhedron."+str.str());
	}

	// Compute the inverse of the linear transformation matrix
	math::vdom_matrix<scalar_type> A_inverse(
			t.get_A().inverse(singular));
	if (singular) {
		throw basic_exception(
				"constr_polyhedron<scalar_type> apply_map not implemented for singular transforms");
		return constr_polyhedron<scalar_type> ();
	} else {
		typename constr_polyhedron<scalar_type>::my_poly_ptr new_poly(
				new typename constr_polyhedron<scalar_type>::my_poly_type(*c.get_constraints()));

		for (typename constr_polyhedron<scalar_type>::my_poly_type::iterator
				it = new_poly->begin(); it != new_poly->end(); ++it) {
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

		return constr_polyhedron<scalar_type> (new_poly);
	}
}
/**
 * Takes a mapped polytope and the map and returns the polytope reverse mapped.
 *
 * @ c The mapped Polytope
 * @ t the map
 * @return The reverse mapped polytope
 */
template<typename scalar_type>
constr_polyhedron<scalar_type> reverse_map(
		const constr_polyhedron<scalar_type>& c, const math::affine_map<
		scalar_type>& t){

	// Test that all variables in c are in the codomain of t
	variable_id_set vis = c.get_variable_ids();
	set_difference_assign(vis,t.codomain().get_variable_ids());
	if (!vis.empty()) {
		std::stringstream ss;
	//	logger::copyfmt_to(ss);
		print_variable_id_set(ss,vis);
		std::stringstream str;
	//	logger::copyfmt_to(str);
		str << "reverse map:" << t << std::endl << "\ndomain: "<< t.codomain() << "\ncodomain: "<< t.codomain();
		str << "polyhedron variables: ";
		print_variable_id_set(str,c.get_variable_ids());
		str << "polyhedron constraints: " << c;
		throw basic_exception(
						"Trying to reverse map a polyhedron with a map that doesn't have variables "+ss.str()+" of the polyhedron."+str.str());
	}

	// Compute the new constraint matrix of the reverse mapped poly
	// A_new = A_c * A_map
	// b_new = b_c - A_c*b_map

	using namespace math;

	vdom_matrix<scalar_type> A_t =t.get_A();

	constr_polyhedron<scalar_type> new_poly;

	for (typename constr_polyhedron<scalar_type>::my_poly_type::const_iterator
			it = c.get_constraints()->begin(); it != c.get_constraints()->end(); ++it) {
		vdom_vector<scalar_type> row = it->get_l().get_vdom_vec();
		scalar_type transf_inh_coeff = it->get_inh_coeff()
				+ scalar_product(row, t.get_b());

		row = row * A_t;
		lin_expression<scalar_type> transf_le(row, transf_inh_coeff);
		lin_constraint<scalar_type> transf_con(transf_le, it->get_sign());
		if (!maybe(transf_con.is_always_satisfied())) {
			new_poly.add_constraint(transf_con, false); // false = no redundancy check
		}
	}

	return new_poly;
}

/**
 * Computes the reverse mapped poly on the constraints of cons
 * on which it is possible and returns it.
 * The constraints of cons on which reverse map could not be
 * applied is also returned.
 * Current implementation checks if the entire cons could be
 * reverse mapped. If yes, the reverse mapped poly is returned
 * in rev_cons_poly and dir_cons_poly is empty. If reverse map
 * could not be applied, then rev_cons_map is empty and the
 * dir_cons_poly = cons.
 *
 * @param cons poly to which reverse map is to be applied.
 * @param map The affine map
 * @param U The input set on the map.
 * @param before_cons reversed mapped poly
 * @param after_cons direct mapped poly
 */

template<typename scalar_type>
void reverse_map_constraints(
		typename constr_polyhedron<scalar_type>::const_ptr cons_ptr,
		const math::affine_map<scalar_type>& map,
		const continuous::support_function_provider::const_ptr U,
		constr_polyhedron<scalar_type>& rev_cons_poly,
		constr_polyhedron<scalar_type>& dir_cons_poly) {
	assert(cons_ptr);

	LOGGER(DEBUG6, __FUNCTION__,
			"reverse mapping constraints");

	// get the domain over which U and b are defined
	positional_vdomain U_dom = map.codomain();

	math::vdom_vector<scalar_type> b(U_dom);

	// @todo filter out constraints that have variables not assigned by the map
	variable_id_set excess_vis = cons_ptr->get_variable_ids();
	set_difference_assign(excess_vis, map.codomain().get_variable_ids());

	// quantify away excess variables
	if (!excess_vis.empty()) {
		LOGGER_OS(DEBUG6,__FUNCTION__) << "quantifying over variables " << excess_vis << std::flush;
		print_variable_id_set(std::cout,excess_vis); std::cout << std::endl << std::flush;
		typename constr_polyhedron<scalar_type>::ptr p = typename constr_polyhedron<scalar_type>::ptr(cons_ptr->clone());
		p->existentially_quantify_variables(excess_vis);
		cons_ptr = p;
	}

	if (is_input_set_point(U, b)) {
		b += map.get_b();

		math::affine_map<scalar_type> new_map(map.get_A(), b);
		rev_cons_poly = reverse_map(*cons_ptr, new_map);
		dir_cons_poly = constr_polyhedron<scalar_type> ();
	} else {
		dir_cons_poly = *cons_ptr;
		rev_cons_poly = constr_polyhedron<scalar_type> ();
	}
}

}

#endif /* CONSTR_POLYHEDRON_OPERATORS_H_ */
