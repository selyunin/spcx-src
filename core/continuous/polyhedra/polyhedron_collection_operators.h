/*
 * polyhedron_collection_operators.h
 *
 *  Created on: Oct 18, 2010
 *      Author: frehse
 */

#ifndef POLYHEDRON_COLLECTION_OPERATORS_H_
#define POLYHEDRON_COLLECTION_OPERATORS_H_

#include "polyhedron_collection.h"
#include "constr_polyhedron/constr_polyhedron.h"
#include "constr_polyhedron/constr_polyhedron_operators.h"

#include "math/vdom/vdom_vector_operators.h"
#include "utility/logger_stopwatch.h"

namespace continuous {

/** Transforms the polyhedra in the collection using the apply_map
 * operator for const_polyhedra.
 *
 * Conserves the type of polyhedron for each element. */
template<typename scalar_type>
polyhedron_collection<scalar_type> apply_map(const polyhedron_collection<
		scalar_type>& polys, const math::affine_map<scalar_type>& t) {
	polyhedron_collection<scalar_type> res;
	for (typename polyhedron_collection<scalar_type>::const_iterator pt =
			polys.begin(); pt != polys.end(); ++pt) {
		// obtain a constr_poly from the constraints of **pt and map it
		constr_polyhedron<scalar_type> cpoly;
		cpoly.add_constraints(*(*pt)->get_constraints());
		cpoly = apply_map(cpoly, t);

		// create a universe poly of the same type as **pt
		typename polyhedron<scalar_type>::ptr new_poly =
				typename polyhedron<scalar_type>::ptr((*pt)->create_universe());

		// add the constraints of cpoly to new_poly
		new_poly->add_constraints(*cpoly.get_constraints());

		// add the new poly to the result, without redundancy checking
		res.insert(new_poly, false);
	}
	return res;
}
;

}

#endif /* POLYHEDRON_COLLECTION_OPERATORS_H_ */
