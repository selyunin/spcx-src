/*
 * polyhedron_utility.h
 *
 *  Created on: Sep 18, 2009
 *      Author: frehse
 */

#ifndef POLYHEDRON_UTILITY_H_
#define POLYHEDRON_UTILITY_H_

#include "core/continuous/polyhedra/polyhedron.h"
#include "math/vdom/lin_expression.h"
#include "math/vdom/lin_constraint.h"

namespace continuous {

/** For every variable x in vis, add the constraint x'==x. */
template<typename scalar_type>
void add_const_relation_constraints(polyhedron<scalar_type>& p, const variable_id_set& vis) {

	for (variable_id_set::const_iterator vit = vis.begin(); vit != vis.end(); ++vit) {
		const variable_id& id = *vit;
		variable_id primed_id = variable::get_id_primedness_increased(id);
		math::lin_expression<scalar_type> e;
		e.set_coeff_with_id(id, scalar_type(1));
		e.set_coeff_with_id(primed_id, scalar_type(-1));
		math::lin_constraint<scalar_type> c(e, EQ);
		p.add_constraint(c);
	}
};

/** Add to p the constraints x==v(x) for every variable x in v. */
template<typename scalar_type>
void add_vertice_constraints(polyhedron<scalar_type>& p,const math::vdom_vector<scalar_type> v) {
	// bring into the same space by embedding variables not in p
	p.embed_variables(v.get_variable_ids());

	for (unsigned int i=0;i<v.get_index_to_variable_id_map()->dimensions();++i) {
		math::lin_expression<scalar_type> e;
		e.set_coeff_with_id(v.get_id(i), scalar_type(-1));
		e.set_inh_coeff(v[i]);
		math::lin_constraint<scalar_type> c(e,EQ);
		p.add_constraint(c);
	}
};

}

#endif /* POLYHEDRON_UTILITY_H_ */
