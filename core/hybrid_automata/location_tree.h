/*
 * location_tree.h
 *
 *  Created on: Jul 25, 2009
 *      Author: gvincent
 */

#ifndef LOCATION_TREE_H_
#define LOCATION_TREE_H_

#include "core/predicates/valuation_function.h"
#include "core/predicates/valuation_function_tree_nodes.h"
#include "core/hybrid_automata/hybrid_automaton.h"
#include "core/hybrid_automata/location_eq_node.h"

namespace hybrid_automata {


template<typename T> class automaton_valuation : public valuation_functions::valuation<T, automaton_id> {
public:

	automaton_valuation(const valuation_functions::valuation<T, automaton_id>& v) : valuation_functions::valuation<T, automaton_id>(v){
	}
	;

	automaton_valuation(){
	}
	;

	virtual ~automaton_valuation(){
	}
	;

	/** Attribute the value val to the variable with name n.
	 * \todo define add (depend of automate)*/
	virtual void add(const std::string& n, const T& val) {
		//valuation<T, variable_id>::add(variable::get_or_add_variable_id(n), val);
	}
	;

	/** Returns the value of the variable with name n.
	 * Throws an exception if the variable has no value attributed to it.
	 * \todo define add (depend of automate)*/
	virtual const T& get_value(const std::string& n) const {
		//return valuation<T, variable_id>::get_value(variable::get_or_add_variable_id(n));
		return valuation_functions::valuation<T, automaton_id>::get_value(0);
	}
	;
};


template <typename bool_type> class location_evaluator {
public:
	virtual ~location_evaluator() {
	}
	;

	virtual bool_type const_bool_node_eval(valuation_functions::const_node<bool_type>* p,
				const typename automaton_valuation<location_id>::const_ptr& v) {
			return (bool_type)(p->my_val);
		}
		;

		virtual bool_type boolean_node_eval(valuation_functions::boolean_node* p,
				const typename automaton_valuation<location_id>::const_ptr& v) {
			bool_type x1=location_eval(p->child1, v);
			if (p->my_op==NOT) {
				return !x1;
			}
			bool_type x2=location_eval(p->child2, v);
			if (p->my_op==AND) {
				return x1 && x2;
			} else if (p->my_op==OR) {
				return x1 || x2;
			} else {
				throw std::runtime_error("unknown boolean operation");
			}
		}
		;

		virtual bool_type location_eq_node_eval(location_eq_node* p,
				const typename automaton_valuation<location_id>::const_ptr& v) {
			if(p->get_equal())
				return (bool_type)(p->get_location_id() == v->get_value(p->get_automaton_id()));
			else
							return (bool_type)(p->get_location_id() != v->get_value(p->get_automaton_id()));
		}
		;

	virtual bool_type location_eval(const tree::node::ptr& p,
			const typename automaton_valuation<location_id>::const_ptr& v) {
		if (valuation_functions::boolean_node* q = dynamic_cast<valuation_functions::boolean_node*>(p.get())) {
					return boolean_node_eval(q, v);
				} else if (location_eq_node* q = dynamic_cast<location_eq_node*>(p.get())) {
					return location_eq_node_eval(q, v);
				} else if (valuation_functions::const_node<bool_type>* q
						= dynamic_cast<valuation_functions::const_node<bool_type>*>(p.get())) {
					return const_bool_node_eval(q, v);
				} else
			throw std::runtime_error("unmatched node in location_eval");
	}
	;

};

template<typename bool_type> class location_tree :
	public shared_ptr_user<location_tree<bool_type> > {
public:
	/** We want to use \p ptr from shared_ptr_user to refer to the base class,
	 * so here we define new pointers for the tree. */
	typedef boost::shared_ptr<location_tree<bool_type> > tree_ptr;
	typedef boost::shared_ptr<const location_tree<bool_type> >
			tree_const_ptr;
	typedef location_evaluator<bool_type> evaluator_type;
	typedef boost::shared_ptr<evaluator_type> evaluator_ptr;

	location_tree(evaluator_ptr eval, const tree::node::ptr& root) :
		my_root(root), my_eval(eval) {
	}
	;
		location_tree(const tree::node::ptr& root) :
		my_root(root) {
		evaluator_ptr p_eval= evaluator_ptr(new evaluator_type);
		my_eval=p_eval;
	}
	;
	virtual ~location_tree() {
	}
	;
	virtual const tree::node::ptr& get_root() {
		return my_root;
	}
	;
	virtual void set_root(const tree::node::ptr& root) {
		my_root=root;
	}
	;
	bool_type eval(const typename automaton_valuation<location_id>::const_ptr& v) const {
		return my_eval->location_eval(my_root, v);
	}
	;
	bool_type print(const typename automaton_valuation<location_id>::const_ptr& v,
			std::ostream& os) const {
		return my_eval->location_eval(my_root, v);
	}
	;
	variable_id_set get_variable_ids() const {
		throw std::runtime_error("get_variables is not implemented");
	}
	;

private:
	tree::node::ptr my_root;
	evaluator_ptr my_eval;
};

}

#endif /* LOCATION_TREE_H_ */
