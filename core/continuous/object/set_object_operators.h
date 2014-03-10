/*
 * set_object_operators.h
 *
 *  Created on: Mar 25, 2010
 *      Author: frehse
 */

#ifndef SET_OBJECT_OPERATORS_H_
#define SET_OBJECT_OPERATORS_H_

#include "math/vdom/vdom_matrix.h"
#include "math/vdom/vdom_vector.h"
#include "set_object.h"
#include "core/continuous/support_function/sf_base/sf_set.h"

namespace continuous {
namespace object {

typedef math::vdom_matrix<rational_type> rational_matrix;
typedef math::vdom_matrix<float_type> float_matrix;
typedef math::vdom_vector<rational_type> rational_vector;
typedef math::vdom_vector<float_type> float_vector;

typedef support_function::sf_set<global_types::float_type> sf_float;
typedef sf_float::vector_set float_directions;
typedef support_function::sf_set<global_types::rational_type> sf_rational;
typedef sf_rational::vector_set rational_directions;

/** Compute the image of an affine transform.
 *
 * The returned set is Y=A*X+b. */
void affine_transform_assign(set_object& X, const rational_matrix& A,
		const rational_vector& b);
void affine_transform_assign(set_object& X, const float_matrix& A,
		const float_vector& b);

/** Compute the existential quantification of X over the set of variables vars. */
void existentially_quantify(set_object& X, const std::set<variable>& vars);

/** Compute the intersection of X and Y. */
set_object intersection(const set_object& X, const set_object& Y);

/** Create a support function representing the set X, the default linear constraints
 * of the set_object class. */
set_object support_function(const set_object& X);

/** Construct a support function for the minkowski_sum of X and Y. */
set_object support_minkowski_sum(const set_object& X, const set_object& Y);

/** Construct a support function for the convex hull of X and Y. */
set_object support_convex_hull(const set_object& X, const set_object& Y);

/** Create a polyhedral overapproximation of X using the linear constraints from Y.
 * Throws if Y has no linear contraints. */
set_object outer_poly(const set_object& X, const set_object& Y);

/** Create a polyhedral overapproximation of X using the directions in the set D. */
set_object outer_poly(const set_object& X, const float_directions& dirs);

/** Create n (approximately) uniform directions over a given domain. */
float_directions uni_directions(const positional_domain& dom, unsigned int n);

/** Create (2*n) box directions over a given domain. */
float_directions box_directions(const positional_domain& dom);

/** Create (2*n^2) octagonal directions over a given domain. */
float_directions oct_directions(const positional_domain& dom);

}
}

std::ostream& operator<<(std::ostream& os, const continuous::object::float_directions& dirs);

#endif /* SET_OBJECT_OPERATORS_H_ */
