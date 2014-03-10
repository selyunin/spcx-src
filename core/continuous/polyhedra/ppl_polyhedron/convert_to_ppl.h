#ifndef CONVERT_TO_PPL_H_
#define CONVERT_TO_PPL_H_

#include <list>
#include <iostream>

#include "core/continuous/polyhedra/ppl_polyhedron/rat_linexpression.h"
#include "core/predicates/valuation_function_tree.h"

/** Forward declarations */
class index_to_variable_id_map;
typedef boost::shared_ptr<const index_to_variable_id_map> index_to_variable_id_map_ptr;
namespace ppl_polyhedron {
class continuous_set_PPL_NNC;
typedef boost::shared_ptr<continuous_set_PPL_NNC> continuous_set_PPL_NNC_ptr;
typedef boost::shared_ptr<const continuous_set_PPL_NNC> continuous_set_PPL_NNC_const_ptr;
}

namespace ppl_polyhedron {

using namespace valuation_functions;

/** Converts a valuation_function_tree into a Rational_Linear_Expression.
 * \attention The valuation_function_tree must be constructed with constant type const_type. */
class Rational_Linear_Expression_generator: public valuation_functions::arithmetic_evaluator<
		Rational_Linear_Expression, Rational> {
public:
	//typedef Rational const_type;
	typedef bool eval_type;
	typedef Rational_Linear_Expression valuation_type;
	typedef bool bool_type;
	typedef Rational_Linear_Expression scalar_type;
	typedef Rational_Linear_Expression RLE_type;
	typedef Parma_Polyhedra_Library::Linear_Expression LE_type;

	virtual ~Rational_Linear_Expression_generator() {
	}
	;

public:
	Rational_Linear_Expression convert(const tree::node::const_ptr& p,
			const index_to_variable_id_map_ptr& iimap);
};

/** Converts a valuation_function_tree into a continuous_set_PPL_NNC object.
 * \attention The valuation_function_tree must be constructed with constant type const_type. */
class PPL_NNC_generator: public valuation_functions::arithmetic_evaluator<
		Rational_Linear_Expression, Rational> {
public:
	typedef Rational const_type;
	typedef bool eval_type;
	typedef Rational_Linear_Expression valuation_type;
	typedef bool bool_type;
	typedef Rational_Linear_Expression scalar_type;
	typedef Rational_Linear_Expression RLE_type;
	typedef Parma_Polyhedra_Library::Linear_Expression LE_type;

	virtual ~PPL_NNC_generator();

	virtual bool_type boolean_node_eval(boolean_node* p,
			const variable_valuation<valuation_type>::const_ptr& v);

	virtual bool_type comparison_node_eval(comparison_node* p, const variable_valuation<
			valuation_type>::const_ptr& v);

	virtual bool_type boolean_eval(const tree::node::ptr& p, const variable_valuation<
			valuation_type>::const_ptr& v);

	virtual eval_type eval(const tree::node::ptr& p,
			const variable_valuation<valuation_type>::const_ptr& v);

public:
	/** Convert with a given index_to_variable_id_map iimap. */
	continuous_set_PPL_NNC_ptr convert(const tree::node::ptr& p,
			const index_to_variable_id_map_ptr& iimap);

	continuous_set_PPL_NNC_ptr convert(const tree::node::ptr& p);

private:
	continuous_set_PPL_NNC* my_ppl_nnc;
};

Rational_Linear_Expression convert_to_Rational_Linear_Expression(const tree::node::const_ptr& p,
		const index_to_variable_id_map_ptr& iimap);
continuous_set_PPL_NNC_ptr convert_to_continuous_set_PPL_NNC(const tree::node::ptr& p,
		const index_to_variable_id_map_ptr& iimap);
continuous_set_PPL_NNC_ptr convert_to_continuous_set_PPL_NNC(const tree::node::ptr& p);
}

#endif /*CONVERT_TO_PPL_H_*/
