/*
 * constr_polyhedron_constructors.hpp
 *
 *  Created on: Dec 29, 2009
 *      Author: frehse
 */

#ifndef CONSTR_POLYHEDRON_CONSTRUCTORS_HPP_
#define CONSTR_POLYHEDRON_CONSTRUCTORS_HPP_

#include "constr_polyhedron_constructors.h"
#include "io/common_input/predicate_parser.h"
#include "io/common_input/bool_node_creation.h"

namespace continuous {

template<typename scalar_type, typename orig_type> typename continuous::constr_polyhedron<
		scalar_type>::ptr construct_constr_polyhedron(
		const continuous::polyhedron<orig_type>& c) {
	// convert the constraints of c to scalar_type
	typename continuous::constr_polyhedron<scalar_type>::my_poly_ptr
			new_cons =
					typename
					continuous::constr_polyhedron<scalar_type>::my_poly_ptr(
							new typename continuous::constr_polyhedron<
									scalar_type>::my_poly_type(
									c.get_constraints()->template convert_to<
											scalar_type> ()));
	typename continuous::constr_polyhedron<scalar_type>::ptr new_c =
			typename continuous::constr_polyhedron<scalar_type>::ptr(
					new continuous::constr_polyhedron<scalar_type>(new_cons));
	return new_c;
}

namespace {
template<typename T>
struct construct_constr_polyhedron_poly_handler {
	template<typename before_t>
	static continuous::continuous_set::ptr implement(
			const continuous::polyhedron<before_t>& c) {
		return construct_constr_polyhedron<T> (c);
	}
	;
};
}

template<typename before_t>
continuous::continuous_set::ptr construct_constr_polyhedron(
		const continuous::polyhedron<before_t>& c,
		global_types::coefficient_type t) {
	return global_types::coefficient_type_caller<
			continuous::continuous_set::ptr, const continuous::polyhedron<
					before_t>&, construct_constr_polyhedron_poly_handler>::call(
			c, t);
}

/** Helper class for the conversion. */
template<typename scalar_type>
class constr_polyhedron_generator: public valuation_functions::arithmetic_evaluator<
		math::lin_expression<scalar_type>, scalar_type> {
public:
	typedef scalar_type const_type;
	typedef bool eval_type;
	typedef scalar_type valuation_type;
	typedef bool bool_type;
	typedef math::lin_expression<scalar_type> LE_type;
	typedef typename valuation_functions::variable_valuation<valuation_type>::const_ptr
			val_const_ptr;

	virtual ~constr_polyhedron_generator() {
	}
	;

	virtual bool_type
	boolean_node_eval(valuation_functions::boolean_node* p,
			const val_const_ptr& v);

	virtual bool_type
	comparison_node_eval(valuation_functions::comparison_node* p,
			const val_const_ptr& v);

	virtual bool_type
	boolean_eval(const tree::node::ptr& p, const val_const_ptr& v);

	virtual eval_type
	eval(const tree::node::ptr& p, const val_const_ptr& v);

public:
	/** Convert with a given index_to_variable_id_map iimap. */
	typename constr_polyhedron<scalar_type>::ptr convert(
			const tree::node::ptr& p);

private:
	typename constr_polyhedron<scalar_type>::ptr my_poly;
};

template<typename scalar_type>
typename constr_polyhedron<scalar_type>::ptr construct_constr_polyhedron(
		const std::string& s, const std::string& context, const parser::parse_policy& ppol) {
	parse_type_chooser::type old_bool=parse_type_chooser::get_bool();
	parse_type_chooser::type old_number=parse_type_chooser::get_number();

	/* Get a predicate representation */
	parse_type_chooser::set_bool(global_types::STD_BOOL);
	parse_type_chooser::set_number(global_types::type_identifier<scalar_type>::coeff);

	// lock the variable creation so that context dependent lookup is done
	bool oldlock=valuation_functions::variable_node_creator::is_locked();
	valuation_functions::variable_node_creator::set_locked(true);

	tree::node_ptr pred_rep = predicate_parser::parse_predicate(s, context, ppol);

	valuation_functions::variable_node_creator::set_locked(oldlock);

	parse_type_chooser::set_bool(old_bool);
	parse_type_chooser::set_number(old_number);

	return construct_constr_polyhedron<scalar_type>(pred_rep);
}


template<typename scalar_type>
typename constr_polyhedron<scalar_type>::ptr construct_constr_polyhedron(
		const std::string& s) {
	parse_type_chooser::type old_bool=parse_type_chooser::get_bool();
	parse_type_chooser::type old_number=parse_type_chooser::get_number();



	/* Get a predicate representation */
	parse_type_chooser::set_bool(global_types::STD_BOOL);
	parse_type_chooser::set_number(global_types::type_identifier<scalar_type>::coeff);

	tree::node_ptr pred_rep = predicate_parser::parse_predicate(s, "");

	parse_type_chooser::set_bool(old_bool);
	parse_type_chooser::set_number(old_number);

	return construct_constr_polyhedron<scalar_type>(pred_rep);
}

template<typename scalar_type>
typename constr_polyhedron<scalar_type>::ptr construct_constr_polyhedron(
		const tree::node::ptr& p) {
	constr_polyhedron_generator<scalar_type> g;
	return g.convert(p);
}

template<typename scalar_type>
typename constr_polyhedron<scalar_type>::ptr construct_constr_polyhedron(
		const tree::node::const_ptr& p) {
	// @todo dirty! clean up constness
	tree::node::ptr pc = boost::const_pointer_cast<tree::node>(p);
	return construct_constr_polyhedron<scalar_type> (pc);
}

namespace {
template<typename T>
struct construct_constr_polyhedron_tree_handler {
	static continuous::continuous_set::ptr implement(
			const tree::node::const_ptr& p) {
		return construct_constr_polyhedron<T> (p);
	}
	;
};
}

inline continuous::continuous_set::ptr construct_constr_polyhedron(
		const tree::node::const_ptr& p, global_types::coefficient_type t) {
	return global_types::coefficient_type_caller<
			continuous::continuous_set::ptr, const tree::node::const_ptr&,
			construct_constr_polyhedron_tree_handler>::call(p, t);
}

template<typename scalar_type>
typename constr_polyhedron_generator<scalar_type>::bool_type constr_polyhedron_generator<
		scalar_type>::boolean_node_eval(valuation_functions::boolean_node* p,
		const val_const_ptr& v) {
	boolean_eval(p->child1, v);
	if (p->my_op == NOT) {
		std::stringstream s;
		s << p;
		throw basic_exception(
				"boolean operation NOT is not allowed: "+s.str());
	}
	boolean_eval(p->child2, v);
	if (p->my_op == AND) {
		return true; // Rational(0);
	} else {
		std::stringstream s;
		valuation_functions::print_node(s,p);
		throw basic_exception("boolean operation "
				+ valuation_functions::operator_to_string(p)
				+ " is not allowed: "+s.str());
	}
}

template<typename scalar_type>
typename constr_polyhedron_generator<scalar_type>::bool_type constr_polyhedron_generator<
		scalar_type>::comparison_node_eval(
		valuation_functions::comparison_node* p, const val_const_ptr& v) {

	math::lin_constraint<scalar_type> con;
	try {
		LE_type x1 = math::convert_to_lin_expression<scalar_type>(p->child1);
		LE_type x2 = math::convert_to_lin_expression<scalar_type>(p->child2);
		con = math::lin_constraint<scalar_type>(x1 - x2, p->my_op);
	} catch (std::exception& e) {
		std::stringstream s;
		valuation_functions::print_node(s,p);
		throw basic_exception("The following is not a linear constraint:\n"+s.str(),e);
	}

	my_poly->add_constraint(con);
	return true;
}

template<typename scalar_type>
typename constr_polyhedron_generator<scalar_type>::bool_type constr_polyhedron_generator<
		scalar_type>::boolean_eval(const tree::node::ptr& p,
		const val_const_ptr& v) {
	if (valuation_functions::boolean_node * q = dynamic_cast<valuation_functions::boolean_node*> (p.get())) {
		return boolean_node_eval(q, v);
	} else if (valuation_functions::comparison_node * q = dynamic_cast<valuation_functions::comparison_node*> (p.get())) {
		return comparison_node_eval(q, v);
	} else if (p == tree::node::null_node()) {
		return (bool_type) (true);
	} else if (valuation_functions::const_node<bool_type> * q
			= dynamic_cast<valuation_functions::const_node<bool_type>*> (p.get())) {
		if (q->my_val) {
			return true; //const_node_eval(q, v);
		} else {
			// the node is false = empty set, so
			// add an unsatisfiable constraint
			math::lin_constraint<scalar_type> c =
					math::lin_constraint<scalar_type>::zero_dim_false();
			my_poly->add_constraint(c);
			return true;
		}
	} else {
		std::stringstream s;
		s << p;
		throw basic_exception("The following boolean expression is not handled (only conjunctions are allowed): " + s.str());
	}
}

template<typename scalar_type>
typename constr_polyhedron_generator<scalar_type>::eval_type constr_polyhedron_generator<
		scalar_type>::eval(const tree::node::ptr& p, const val_const_ptr& v) {
	if (valuation_functions::boolean_node * q = dynamic_cast<valuation_functions::boolean_node*> (p.get())) {
		return boolean_node_eval(q, v);
	} else if (valuation_functions::comparison_node * q = dynamic_cast<valuation_functions::comparison_node*> (p.get())) {
		return comparison_node_eval(q, v);
	} else if (valuation_functions::const_node<bool_type> * q
			= dynamic_cast<valuation_functions::const_node<bool_type>*> (p.get())) {
		if (q->my_val) {
			return true; //const_node_eval(q, v);
		} else {
			// the node is false = empty set, so
			// add an unsatisfiable constraint
			math::lin_constraint<scalar_type> c =
					math::lin_constraint<scalar_type>::zero_dim_false();
			my_poly->add_constraint(c);
			return true;
		}

	}
	//Todo: check that is correct
	else if (!p || p == tree::node::null_node()) {
		return (eval_type) (true);
	} else {
		std::stringstream s;
		s << p;
		throw basic_exception("The following expression is not handled (only conjunctions are allowed): " + s.str());
	}
}

template<typename scalar_type>
typename constr_polyhedron<scalar_type>::ptr constr_polyhedron_generator<
		scalar_type>::convert(const tree::node::ptr& p) {
	// Create a universe polyhedron to which we will add constraints
	my_poly = typename constr_polyhedron<scalar_type>::ptr(
			new constr_polyhedron<scalar_type> ());

	// Create the valuation that maps each variable_id to the corresponding PPL-variable
	typename valuation_functions::variable_valuation<valuation_type>::ptr
			v =
					typename
					valuation_functions::variable_valuation<valuation_type>::ptr(
							new valuation_functions::variable_valuation<
									valuation_type>());

	// Run the evaluation
	try {
		eval(p, v);
	} catch (std::exception& e) {
		throw basic_exception("The expression is not a conjunction of linear constraints.",e);
	}

	return my_poly;
}

}

#endif /* CONSTR_POLYHEDRON_CONSTRUCTORS_HPP_ */
