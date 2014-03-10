#include "core/hybrid_automata/assignment_traitment.h"

#include "utility/tree_traverser.h"

namespace assignment_traitment {

class has_var: public tree::tree_traverser {
public:
	has_var() {
		var = false;
	}

	bool eval(tree::node::ptr p) {
		if (dynamic_cast<valuation_functions::variable_node*> (p.get())) {
			var = true;
			return false;
		}
		return true;
	}
	bool var;
};

bool has_variable(const tree::node::ptr& p) {
	has_var g;
	g.traverse_preorder(p);
	return g.var;
}

class all_define: public tree::tree_traverser {
public:
	all_define(const index_to_variable_id_map_ptr& index_id_map) :
		iimap(index_id_map), define(index_id_map->dimensions(), false) {
	}

	bool eval(tree::node::ptr p) {
		if (valuation_functions::comparison_node* c = dynamic_cast<valuation_functions::comparison_node*>(p.get()))
			if (valuation_functions::variable_node* v = dynamic_cast<valuation_functions::variable_node*>(c->child1.get()))
				define[iimap->get_index(
						variable::get_id_primedness_decreased(v->my_id))]
						= true;
		return true;
	}

	bool test() {
		bool all_v_define = true;
		for (unsigned int i = 0; i < define.size(); ++i)
			if (define[i] == false)
				all_v_define = false;
		return all_v_define;
	}
	index_to_variable_id_map_ptr iimap;
	math::vector<bool> define;
};

bool all_variable_define(const tree::node::ptr& p) {
	index_to_variable_id_map_ptr iimap(index_to_variable_id_map::empty_map());
	iimap = iimap->get_map_with_ids_added(
			valuation_functions::get_unprimed_variable_ids(p));
	all_define g(iimap);
	g.traverse_preorder(p);
	return g.test();
}

class add_Ode: public tree::tree_traverser {
public:
	add_Ode() {
	}

	bool eval(tree::node::ptr p) {
		if (valuation_functions::comparison_node* q = dynamic_cast<valuation_functions::comparison_node*>(p.get())) {
			if (valuation_functions::variable_node* v
					= dynamic_cast<valuation_functions::variable_node*>(q->child1.get())) {
				add.insert(std::make_pair(v->my_id, q->child2));
			} else {
				std::stringstream s;
				logger::copyfmt_to(s);
				s << p;
				throw basic_exception("Cannot handle "+s.str()+".");
			}
		}
		return true;
	}
	std::map<variable_id, tree::node::ptr> add;
};

std::map<variable_id, tree::node::ptr> get_all_Ode(const tree::node::ptr& p) {
	add_Ode g;
	g.traverse_preorder(p);
	return g.add;
}

class all_constraint: public tree::tree_traverser {
public:
	all_constraint() {
	}

	bool eval(tree::node::ptr p) {
		if (dynamic_cast<valuation_functions::comparison_node*> (p.get())) {
			all.push_back(p);
		}
		return true;
	}
	std::vector<tree::node::ptr> all;
};

std::vector<tree::node::ptr> get_all_constraint(const tree::node::ptr& p) {
	all_constraint g;
	g.traverse_preorder(p);
	return g.all;
}

class all_primed: public tree::tree_traverser {
public:
	all_primed() {
		all = true;
	}

	bool eval(tree::node::ptr p) {
		if (valuation_functions::variable_node* q = dynamic_cast<valuation_functions::variable_node*>(p.get())) {
			if (variable::get_prime_count(q->my_id) < 1) {
				all = false;
				return false;
			}
		}
		return true;
	}
	bool all;
};

bool is_LHA_dynamics(const tree::node::ptr& p) {
	all_primed g;
	g.traverse_preorder(p);
	return g.all;
}

class is_all_Ode: public tree::tree_traverser {
public:
	is_all_Ode() {
		all = true;
	}

	bool eval(tree::node::ptr p) {
		bool is_ode = true;
		if (valuation_functions::comparison_node* q = dynamic_cast<valuation_functions::comparison_node*>(p.get())) {
			// Comparisons need to be equalities
			if (q->my_op != EQ) {
				is_ode = false;
			} else {
				if (valuation_functions::variable_node* v
					= dynamic_cast<valuation_functions::variable_node*>(q->child1.get())) {
					// left side is a variable,
					// if it is primed then right side must not be primed
					if (variable::get_prime_count(v->my_id) > 0
							&& valuation_functions::has_primed_variable(
									q->child2)) {
						is_ode = false;
					}
				} else if (valuation_functions::variable_node* v
						= dynamic_cast<valuation_functions::variable_node*>(q->child2.get())) {
					// right side is a variable,
					// if it is primed then left side must not be primed
					if (variable::get_prime_count(v->my_id) > 0
							&& valuation_functions::has_primed_variable(
									q->child1)) {
						is_ode = false;
					}
				} else
					is_ode = false;
			}
		}
		all = all && is_ode;
		return is_ode;
	}
	bool all;
};

bool is_Ode_dynamics(const tree::node::ptr& p) {
	is_all_Ode g;
	g.traverse_preorder(p);
	return g.all;
}

/** Checks whether an expression is affine
 *
 * @remark The tree traverser stops when eval is false.
 * The result of the affine test is stored in affine.
 */
class is_affine: public tree::tree_traverser {
public:
	is_affine() {
		affine = true;
	}

	bool eval(tree::node::ptr p) {
		if (valuation_functions::arithmetic_node* q =
				dynamic_cast<valuation_functions::arithmetic_node*>(p.get())) {
			if (q->my_op == ADD || q->my_op == SUB || q->my_op == NEG) {
				// don't do anything, children will be checked
			} else if (q->my_op == DIV) {
				affine = affine && !(has_variable(q->child2));
			} else if (q->my_op == MUL) {
				affine = affine && !(has_variable(q->child1) && has_variable(q->child2));
			} else {
				// in all other cases, consider as nonlinear is variables are present
				// if there are children then they must have no variables
				affine = affine && (!q->child1 || !has_variable(q->child1));
				affine = affine && (!q->child2 || !has_variable(q->child2));
			}
		}
		return affine; // continue traversing if still affine (check children)
	}
	bool affine;
};

bool is_affine_tree(const tree::node::ptr& p) {
	is_affine g;
	g.traverse_preorder(p);
	return g.affine;
}

bool is_false_cheap(const tree::node::ptr& p) {
	if (valuation_functions::const_node<bool>* q = dynamic_cast<valuation_functions::const_node<bool>*>(p.get())) {
		return !q->my_val;
	}
	return false;
}

}
