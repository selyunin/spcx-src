/*
 * poly_sphere.h
 *
 *  Created on: Apr 27, 2010
 *      Author: frehse
 */

#ifndef POLY_SPHERE_H_
#define POLY_SPHERE_H_

#include "constr_polyhedron.h"
#include "core/continuous/support_function/template_directions/choose_directions.h"

namespace continuous {

/** Construct a polyhedral overapproximation of a sphere.
 *
 * The polyhedron is defined over the domain dom.
 */
template<typename scalar_type>
constr_polyhedron<scalar_type> poly_sphere(const positional_vdomain& dom,
		unsigned int nb) {
	typedef std::list<math::vector<scalar_type> > vdir_type;
	vdir_type vdirs;

	unsigned int dim = dom.size();
	support_function::add_uniform_directions<scalar_type>(dim, nb,
			vdirs);
	assert(vdirs.size()==nb);

	math::vdom_vector<scalar_type> named_v;
	math::lin_constraint<scalar_type> con;
	continuous::constr_polyhedron<scalar_type> cons_poly; // universe set since no constraints.
	for (typename vdir_type::const_iterator it = vdirs.begin(); it != vdirs.end(); ++it) {
		assert(it->size()==dom.size());
		named_v = math::vdom_vector<scalar_type>(dom, *it);
		double norm=std::sqrt(convert_element<double>(scalar_product((*it),(*it))));
		math::lin_expression<scalar_type> lin_exp(named_v, -scalar_type(norm));
		con = math::lin_constraint<scalar_type>(lin_exp, LE);
		cons_poly.add_constraint(con);
	}

	assert(cons_poly.get_constraints()->size()==nb);
	assert(cons_poly.get_dim()==dom.size());

	return cons_poly;
}

}

#endif /* POLY_SPHERE_H_ */
