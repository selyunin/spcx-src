/*
 * convert_to_lin_expression.h
 *
 *  Created on: Nov 9, 2009
 *      Author: frehse
 */

#ifndef CONVERT_TO_LIN_EXPRESSION_H_
#define CONVERT_TO_LIN_EXPRESSION_H_

#include <sstream>
#include "math/vdom/lin_expression.h"
#include "core/predicates/valuation_function_tree.h"

namespace math {

/** Helper class for conversion */
template<typename scalar_type>
class lin_expression_generator: public valuation_functions::arithmetic_evaluator<
		lin_expression<scalar_type> , scalar_type> {
public:
	//typedef Rational const_type;
	typedef bool eval_type;
	typedef lin_expression<scalar_type> valuation_type;
	typedef bool bool_type;
	typedef lin_expression<scalar_type> LE_type;

	virtual ~lin_expression_generator() {
	}
	;

public:
	LE_type convert(const tree::node::const_ptr& p,
			const index_to_variable_id_map_ptr& iimap);
};

template<typename scalar_type>
typename lin_expression_generator<scalar_type>::LE_type lin_expression_generator<scalar_type>::convert(
		const tree::node::const_ptr& p,
		const index_to_variable_id_map_ptr& iimap) {
	// Create the valuation that maps each variable_id to the corresponding PPL-variable
	typename valuation_functions::variable_valuation<LE_type>::ptr v = typename valuation_functions::variable_valuation<LE_type>::ptr(
			new valuation_functions::variable_valuation<LE_type> ());
	// iterate over the index
	for (index_type i = 0; i < iimap->dimensions(); ++i) {
		LE_type lex(iimap);
		lex[i] = scalar_type(1);
		v->add(iimap->get_id(i), lex);
	}
	// Run the evaluation
	// @todo use an actual const evaluator!!!
	tree::node::ptr pconst = boost::const_pointer_cast<tree::node>(p);
	typename lin_expression_generator<scalar_type>::LE_type result;
//	try {
		result
				= valuation_functions::arithmetic_evaluator<LE_type,
						scalar_type>::arithmetic_eval(pconst, v);
//	} catch (std::exception& e) {
//		std::stringstream s;
//		s << p;
//		throw basic_exception("The following is not a linear expression: "
//				+ s.str(), e);
//	}
	return result;
}

/** Converts a valuation_function_tree into a lin_expression of scalar_type
 * with a given index_to_variable_id_map. */
template<typename scalar_type>
lin_expression<scalar_type> convert_to_lin_expression(
		const tree::node::const_ptr& p,
		const index_to_variable_id_map_ptr& iimap) {

	// get the linear expression corresponding to the function
	lin_expression_generator<scalar_type> gen;
	lin_expression<scalar_type> e = gen.convert(p, iimap);
	return e;
}

/** Converts a valuation_function_tree into a lin_expression of scalar_type. */
template<typename scalar_type>
lin_expression<scalar_type> convert_to_lin_expression(
		const tree::node::const_ptr& p) {

	// Obtain all variables and create an iimap from it
	variable_id_set vis = valuation_functions::get_variable_ids(p);
	//std::cout << vis;

	index_to_variable_id_map_ptr iimap(index_to_variable_id_map::empty_map());
	iimap = iimap->get_map_with_ids_added(vis);
	//std::cout << iimap;

	return convert_to_lin_expression<scalar_type>(p, iimap);
}

}

#endif /* CONVERT_TO_LIN_EXPRESSION_H_ */
