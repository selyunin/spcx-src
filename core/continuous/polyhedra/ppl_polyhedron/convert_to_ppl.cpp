/*
 * convert_to_ppl.cpp
 *
 *  Created on: Sep 7, 2009
 *      Author: frehse
 */

#include "core/continuous/polyhedra/ppl_polyhedron/convert_to_ppl.h"

#include "utility/stl_helper_functions.h"
#include "core/continuous/polyhedra/ppl_polyhedron/continuous_set_PPL_NNC.h"
#include "core/predicates/valuation_function_tree_utility.h"
#include "utility/shared_ptr_output.h"
#include "core/predicates/node_print_visitor.h"
#include "io/common_input/bool_node_creation.h"

namespace ppl_polyhedron {

Rational_Linear_Expression Rational_Linear_Expression_generator::convert(
		const tree::node::const_ptr& p, const index_to_variable_id_map_ptr& iimap) {
	// Create the valuation that maps each variable_id to the corresponding PPL-variable
	variable_valuation<RLE_type>::ptr v = variable_valuation<RLE_type>::ptr(new variable_valuation<
			RLE_type> ());
	// iterate over the index
	for (index_type i = 0; i < iimap->dimensions(); ++i) {
		v->add(iimap->get_id(i), (RLE_type) (Parma_Polyhedra_Library::Variable(i)));
	}
	// Run the evaluation
	// @todo use an actual const evaluator!!!
	tree::node::ptr pconst = boost::const_pointer_cast<tree::node>(p);
	return valuation_functions::arithmetic_evaluator<Rational_Linear_Expression, Rational>::arithmetic_eval(
			pconst, v);
}

Rational_Linear_Expression convert_to_Rational_Linear_Expression(const tree::node::const_ptr& p,
		const index_to_variable_id_map_ptr& iimap) {

	// get the linear expression corresponding to the function
	Rational_Linear_Expression_generator gen;
	Rational_Linear_Expression e = gen.convert(p, iimap);
	return e;
}

// ------------------------------------------------------------

PPL_NNC_generator::~PPL_NNC_generator() {
}

PPL_NNC_generator::bool_type PPL_NNC_generator::boolean_node_eval(boolean_node* p,
		const variable_valuation<valuation_type>::const_ptr& v) {
	boolean_eval(p->child1, v);
	if (p->my_op == NOT) {
		throw std::runtime_error("boolean operation NOT is not allowed in NNC_PPL");
	}
	boolean_eval(p->child2, v);
	if (p->my_op == AND) {
		return true; // Rational(0);
	} else {
		throw std::runtime_error("boolean operation "+valuation_functions::operator_to_string(p)+" is not allow in NNC_PPL");
	}
}

PPL_NNC_generator::bool_type PPL_NNC_generator::comparison_node_eval(comparison_node* p,
		const variable_valuation<valuation_type>::const_ptr& v) {
	scalar_type x1 = valuation_functions::arithmetic_evaluator<Rational_Linear_Expression,
			const_type>::arithmetic_eval(p->child1, v);
	scalar_type x2 = valuation_functions::arithmetic_evaluator<Rational_Linear_Expression,
			const_type>::arithmetic_eval(p->child2, v);
	//		std::cout << " comp:"<< x1 << std::endl;
	//		std::cout << " with:"<< x2 << std::endl;
	Integer den1 = x1.get_denominator();
	Integer den2 = x2.get_denominator();
	Parma_Polyhedra_Library::Constraint c(Parma_Polyhedra_Library::Constraint::zero_dim_false()); // dummy constraint
	if (p->my_op == LT) {
		c = (x1.get_LE() * den2 < x2.get_LE() * den1);
	} else if (p->my_op == LE) {
		c = (x1.get_LE() * den2 <= x2.get_LE() * den1);
	} else if (p->my_op == GT) {
		c = (x1.get_LE() * den2 > x2.get_LE() * den1);
	} else if (p->my_op == GE) {
		c = (x1.get_LE() * den2 >= x2.get_LE() * den1);
	} else if (p->my_op == EQ) {
		c = (x1.get_LE() * den2 == x2.get_LE() * den1);
	} else {
		throw std::runtime_error("unknown comparison operator "+valuation_functions::operator_to_string(p));
	}
	//std::cout << c;
	//ppl_polyhedron::print_constraint(c, "", cout);
	my_ppl_nnc->add_constraint(c);
	return true;
}

PPL_NNC_generator::bool_type PPL_NNC_generator::boolean_eval(const tree::node::ptr& p,
		const variable_valuation<valuation_type>::const_ptr& v) {
	if (boolean_node * q = dynamic_cast<boolean_node*> (p.get())) {
		return boolean_node_eval(q, v);
	} else if (comparison_node * q = dynamic_cast<comparison_node*> (p.get())) {
		return comparison_node_eval(q, v);
		/* } Let's forbid bool constants and variables
		 else if (const_node<const_type>* q
		 = dynamic_cast<const_node<const_type>*>(p.get())) {
		 return const_node_eval(q, v);
		 } else if (variable_node* q = dynamic_cast<variable_node*>(p.get())) {
		 return variable_node_eval(q, v);
		 */
		//Todo: check that is correct
	} else if (p == tree::node::null_node()) {
		return (bool_type) (true);
	} else if (const_node<bool_type> * q = dynamic_cast<const_node<bool_type>*> (p.get())) {
		if (q->my_val) {
			return true; //const_node_eval(q, v);
		} else {
			// the node is false = empty set, so
			// add an unsatisfiable constraint
			Parma_Polyhedra_Library::Constraint c(Parma_Polyhedra_Library::Constraint::zero_dim_false());
			my_ppl_nnc->add_constraint(c);
			return true;
		}
	} else {
		std::stringstream s;
		s << p;
		throw basic_exception("The following boolean expression is not handled (only conjunctions are allowed): " + s.str());
	}
}

PPL_NNC_generator::eval_type PPL_NNC_generator::eval(const tree::node::ptr& p,
		const variable_valuation<valuation_type>::const_ptr& v) {
	if (boolean_node * q = dynamic_cast<boolean_node*> (p.get())) {
		return boolean_node_eval(q, v);
	} else if (comparison_node * q = dynamic_cast<comparison_node*> (p.get())) {
		return comparison_node_eval(q, v);
	} else if (const_node<bool_type> * q = dynamic_cast<const_node<bool_type>*> (p.get())) {
		if (q->my_val) {
			return true; //const_node_eval(q, v);
		} else {
			// the node is false = empty set, so
			// add an unsatisfiable constraint
			Parma_Polyhedra_Library::Constraint c(Parma_Polyhedra_Library::Constraint::zero_dim_false());
			my_ppl_nnc->add_constraint(c);
			return true;
		}

	}
	//Todo: check that is correct
	else if (!p || p == tree::node::null_node()) {
		return (eval_type) (true);
	} else {
		std::stringstream s;
		s << p;
		throw basic_exception("The following boolean expression is not handled (only conjunctions are allowed): " + s.str());
	}
}

continuous_set_PPL_NNC::ptr PPL_NNC_generator::convert(const tree::node::ptr& p,
		const index_to_variable_id_map_ptr& iimap) {
	// Create a universe polyhedron to which we will add constraints
	continuous_set_PPL_NNC::ptr cs = continuous_set_PPL_NNC::ptr(new continuous_set_PPL_NNC(iimap));
	my_ppl_nnc = cs.get();

	// Create the valuation that maps each variable_id to the corresponding PPL-variable
	variable_valuation<RLE_type>::ptr v = variable_valuation<RLE_type>::ptr(new variable_valuation<
			RLE_type> ());
	// iterate over the index
	for (index_type i = 0; i < iimap->dimensions(); ++i) {
		v->add(iimap->get_id(i), (RLE_type) (Parma_Polyhedra_Library::Variable(i)));
	}

	// Run the evaluation
	eval(p, v);

	return cs;
}

continuous_set_PPL_NNC::ptr PPL_NNC_generator::convert(const tree::node::ptr& p) {
	// Obtain all variables and create an iimap from it
	variable_id_set vis = get_variable_ids(p);
	//std::cout << vis;

	index_to_variable_id_map_ptr iimap(index_to_variable_id_map::empty_map());
	iimap = iimap->get_map_with_ids_added(vis);
	//std::cout << iimap;

	continuous_set_PPL_NNC::ptr cs = convert(p, iimap);
	return cs;
}

continuous_set_PPL_NNC::ptr convert_to_continuous_set_PPL_NNC(const tree::node::ptr& p,
		const index_to_variable_id_map_ptr& iimap) {
	PPL_NNC_generator g;
	return g.convert(p, iimap);
}

continuous_set_PPL_NNC::ptr convert_to_continuous_set_PPL_NNC(const tree::node::ptr& p) {
	PPL_NNC_generator g;
	return g.convert(p);
}
}
