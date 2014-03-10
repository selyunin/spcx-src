/*
 * constr_polyhedron_constructors.h
 *
 *  Created on: Dec 29, 2009
 *      Author: frehse
 */

#ifndef CONSTR_POLYHEDRON_CONSTRUCTORS_H_
#define CONSTR_POLYHEDRON_CONSTRUCTORS_H_

#include "core/continuous/polyhedra/constr_polyhedron/constr_polyhedron.h"
#include "math/vdom/convert_to_lin_expression.h"
#include "core/predicates/node_print_visitor.h"

namespace continuous {

/** Construct a constr_polyhedron from a polyhedron with given scalar_type.
 */
template<typename scalar_type, typename orig_type> typename continuous::constr_polyhedron<
		scalar_type>::ptr construct_constr_polyhedron(
		const continuous::polyhedron<orig_type>& c);

/** Construct a constr_polyhedron from a polyhedron with a scalar_type given by
 * a global_types::coeffcient_type.
 */
template<typename before_t>
continuous::continuous_set::ptr construct_constr_polyhedron(
		const continuous::polyhedron<before_t>& c,
		global_types::coefficient_type t);

/** Construct a constr_polyhedron from a textual expression.
 *
 * Looking up only existing symbols within context.
 * Throws if symbol not already present. */
template<typename scalar_type>
typename constr_polyhedron<scalar_type>::ptr construct_constr_polyhedron(
		const std::string& s, const std::string& context, const parser::parse_policy& ppol = parser::parse_policy());

/** Construct a constr_polyhedron from a textual expression.
 *
 * Symbols are created if not already present. */
template<typename scalar_type>
typename constr_polyhedron<scalar_type>::ptr construct_constr_polyhedron(
		const std::string& s);

/** Converts a valuation_function_tree into a constr_polyhedron object of scalar_type. */
template<typename scalar_type>
typename constr_polyhedron<scalar_type>::ptr construct_constr_polyhedron(
		const tree::node::ptr& p);

/** Converts a valuation_function_tree into a lin_expression of scalar_type. */
template<typename scalar_type>
typename constr_polyhedron<scalar_type>::ptr construct_constr_polyhedron(
		const tree::node::const_ptr& p);

/** Construct a constr_polyhedron from a polyhedron with a scalar_type given by
 * a global_types::coeffcient_type.
 */
inline continuous::continuous_set::ptr construct_constr_polyhedron(
		const tree::node::const_ptr& p, global_types::coefficient_type t);

}

#include "constr_polyhedron_constructors.hpp"

#endif /* CONSTR_POLYHEDRON_CONSTRUCTORS_H_ */
