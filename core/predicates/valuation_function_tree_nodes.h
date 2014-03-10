#ifndef VALUATION_FUNCTION_TREE_NODES_H_
#define VALUATION_FUNCTION_TREE_NODES_H_

#include <iostream>
#include <sstream>
#include "utility/operator_enums.h"
#include "utility/stl_helper_functions.h"
#include "utility/tree_node.h"
#include "math/vdom/variable.h"
#include "io/common_input/parse_policy.h"
//#include "valuation_function.h"

#include "core/predicates/dot_context_lookup.h"

namespace valuation_functions {

template<typename ret_type> class const_node: public tree::node {
public:
	const_node(const ret_type& x) :
		my_val(x) {
	}
	;
	virtual ~const_node() {
	}
	;
	virtual void accept(tree::node_visitor& v) const {
		v.visit(this);
		;
	}
	;
	//private:
	ret_type my_val;
};

class variable_node: public tree::node {
public:
	variable_node(const std::string& name, unsigned int dim = 0) {
		my_id = variable::get_or_add_variable_id(name, dim);
		dim1 = dim;
		dim2 = 0;
		if (dim1 > 0)
			dim2 = 1;
	}
	;
	variable_node(const std::string& name, unsigned int d1, unsigned int d2) {
		my_id = variable::get_or_add_variable_id(name, d1 * d2);
		dim1 = d1;
		dim2 = d2;
	}
	;
	variable_node(variable_id id, unsigned int d1, unsigned int d2) {
		my_id = id;
		dim1 = d1;
		dim2 = d2;
	}
	;
	variable_node(variable_id id) {
		my_id = id;
		dim1 = 0;
		dim2 = 0;
		variable var(id);
		if (var.get_dimension() > 0) {
			dim1 = var.get_dimension();
			dim2 = 1;
		}
	}
	;
	variable_node(const variable& var) {
		my_id = var.get_id();
		dim1 = 0;
		dim2 = 0;
		if (var.get_dimension() > 0) {
			dim1 = var.get_dimension();
			dim2 = 1;
		}
	}
	virtual ~variable_node() {
	}
	;
	virtual void accept(tree::node_visitor& v) const {
		v.visit(this);
	}
	;
	const unsigned int get_dim1() const {
		return dim1;
	}
	;
	const unsigned int get_dim2() const {
		return dim2;
	}
	;
	//private:
	variable_id my_id;
private:
	unsigned int dim1;
	unsigned int dim2;
};

/** A class that controls the creation/lookup of variables.
 *
 * Having it as a class gives us the ability to add
 * static memory, e.g., for locking the creation of variables.
 *
 * When parsing a model file, variables are already created
 * when the symbol table is created, so we need lookup only.
 * A variable that is not found should generate an exception.
 *
 * However, when we want to parse a simple string inside
 * a tester or elsewhere, we might not want to bother
 * with a symbol table etc. For this the creation
 * should be unlocked.
 * */
class variable_node_creator {
public:
	/** Create a scalar or vector variable */
	static variable_node* create(const std::string& name,
			const std::string& context, unsigned int dim1 = 0, const parser::parse_policy& ppol = parser::parse_policy()) {
		unsigned int dim2 = 0;
		if (dim1 >= 1)
			dim2 = 1;
		return create(name,context,dim1,dim2,ppol);
	}
	;
	/** Create a variable of given dimension */
	static variable_node* create(const std::string& name,
			const std::string& context, unsigned int d1, unsigned int d2, const parser::parse_policy& ppol = parser::parse_policy()) {
		std::string namec = name_dot_context(name, context);
		if (locked) {
			namec = context_lookup(name, context);
		}
		unsigned int dim1 = d1;
		unsigned int dim2 = d2;
		if (dim1 <= 1 && dim2 <= 1) {
			dim1 = ppol.scalar_dim;
			dim2 = ppol.scalar_dim;
		}
		variable_id id = variable::get_or_add_variable_id(namec, dim1 * dim2);
		return new variable_node(id, dim1, dim2);
	}
	;
	static bool is_locked() {
		return locked;
	}
	static void set_locked(bool l) {
		locked = l;
	}
	/** Return the name plus context */
	static std::string name_dot_context(const std::string& name,
			const std::string& context) {
		if (context.empty()) {
			return name;
		} else {
			return context + "." + name;
		}
	}
	static std::string context_lookup(std::string name,
			const std::string& context) {
		// get the inside context name to treat cases where name already includes the context
		name = dot_context::in_context_name(name, context);
		std::string namec = name_dot_context(name, context);
		// try lookup inside context
		if (!name.empty() && !variable::has_variable(namec)) {
			std::set<std::string> candidates = dot_context::lookup(name,
					context, variable::get_names().begin(),
					variable::get_names().end());
			if (candidates.size() == 1) {
				namec = name_dot_context(*candidates.begin(), context);
			} else if (candidates.size() == 0) {
				throw basic_exception("Could not find variable " + name
						+ " in component " + context + ".");
			} else {
				std::string s;
				for (std::set<std::string>::const_iterator it =
						candidates.begin(); it != candidates.end(); ++it) {
					if (it != candidates.begin()) {
						s += ", ";
					}
					s += *it;
				}
				throw basic_exception("Could not find variable " + name
						+ " in component " + context + ". Found instead " + s
						+ ". Please use a more explicit identifier.");
			}
		}
		return namec;
	}
private:

	static bool locked;
};

class arithmetic_node: public tree::binary_node {
public:
	arithmetic_node(arithmetic_operator op, const tree::node::ptr& child_one,
			const tree::node::ptr& child_two) :
		binary_node(child_one, child_two), my_op(op) {
	}
	;
	arithmetic_node(arithmetic_operator op, const tree::node::ptr& child_one) :
		binary_node(child_one, tree::node::null_node()), my_op(op) {
	}
	;
	virtual void accept(tree::node_visitor& v) const {
		v.visit(this);
	}
	;
	arithmetic_operator my_op;
};

class boolean_node: public tree::binary_node {
public:
	boolean_node(boolean_operator op, const tree::node::ptr& child_one,
			const tree::node::ptr& child_two) :
		binary_node(child_one, child_two), my_op(op) {
	}
	;
	boolean_node(boolean_operator op, const tree::node::ptr& child_one) :
		binary_node(child_one, tree::node::null_node()), my_op(op) {
	}
	;
	virtual void accept(tree::node_visitor& v) const {
		v.visit(this);
	}
	;
	boolean_operator my_op;
};

class comparison_node: public tree::binary_node {
public:
	comparison_node(comparison_operator op, const tree::node::ptr& child_one,
			const tree::node::ptr& child_two) :
		binary_node(child_one, child_two), my_op(op) {
	}
	;
	virtual void accept(tree::node_visitor& v) const {
		v.visit(this);
	}
	;
	comparison_operator my_op;
};

}

#endif /*VALUATION_FUNCTION_TREE_NODES_H_*/

