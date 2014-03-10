#ifndef convert_matrix_math_H_
#define convert_matrix_math_H_

#include <stdexcept>
#include "math/vdom/lin_expression.h"
#include "math/vdom/lin_constraint.h"
//#include "extended_ppl.h"
//#include "fp_interface.h"
#include <ppl.hh>
//#include "rat_linexpression.h"

namespace ppl_polyhedron {

/** Converts a Generator object to a rational vector
 */
math::rational_vector convert_to_rational_vector(const Parma_Polyhedra_Library::Generator& g);

/** Convert a PPL-style Linear_Expression to a lin_expression. */
math::lin_expression<Rational> convert_to_lin_expression(
		const Parma_Polyhedra_Library::Generator& g, const index_to_variable_id_map_ptr& iimap);

/** Convert a PPL-style Linear_Expression to a vdom_vector. */
math::vdom_vector<Rational> convert_to_vdom_vector(
		const Parma_Polyhedra_Library::Generator& g, const index_to_variable_id_map_ptr& iimap);


/** Convert a PPL-style Constraint to a lin_constraint. */
math::lin_constraint<Rational> convert_to_lin_constraint(
		const Parma_Polyhedra_Library::Constraint& c, const index_to_variable_id_map_ptr& iimap);

/**
 * Converts a Generator object to a double vector
 */

math::double_vector convert_to_double_vector(const Parma_Polyhedra_Library::Generator& g);

/** Converts the rational_vector v into a linear expression le and a denominator d such that le[i]/d = v[i]. */
void convert_to_Linear_Expression(const math::rational_vector& v,
		Parma_Polyhedra_Library::Linear_Expression& le, Integer& d);

/** Converts the rational_vector v into a linear expression le and a denominator d such that le[i]/d = v[i]. */
void convert_to_Linear_Expression(const math::double_vector& v,
		Parma_Polyhedra_Library::Linear_Expression& le, Integer& d);

/** Convert a lin_expression to a PPL-style Linear_Expression.
 * Returns true if l has non-zero coefficients for variables that are not in iimap. */
bool convert_to_Linear_Expression(
		const math::lin_expression<Rational>& l, Parma_Polyhedra_Library::Linear_Expression& le,
		Integer& denominator, const index_to_variable_id_map_ptr& iimap);

}

#endif /*convert_matrix_math_H_*/
