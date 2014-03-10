#ifndef ASSIGNMENT_TRAITMENT_H_
#define ASSIGNMENT_TRAITMENT_H_

#include "core/predicates/valuation.h"
#include "core/predicates/valuation_function_tree_utility.h"
#include "core/predicates/valuation_function_tree_nodes.h"
#include "math/vdom/variable.h"
#include "utility/share_tree.h"
#include "math/vdom/index_to_variable_id_map.h"
#include "math/scalar_types/rational.h"
#include "math/vdom/vdom_matrix.h"
#include "math/vdom/vdom_vector.h"
#include "math/type_conversion.h"

namespace assignment_traitment {

/**
 * Use class has_var to determine if the tree contains a variable
 */
bool has_variable(const tree::node::ptr& p);

/**
 * Use class all_variable_define to determine if all the variables of the tree
 * are define
 */
bool all_variable_define(const tree::node::ptr& p);

/**
 * \return a std::map<variable_id, tree::node::ptr> with all Ode equation in the tree
 */
std::map<variable_id, tree::node::ptr> get_all_Ode(const tree::node::ptr& p);

/**
 * \return a vector<tree::node::ptr> with all constraint in the tree
 */
std::vector<tree::node::ptr> get_all_constraint(const tree::node::ptr& p);

/**
 * Use class all_primed to check that all the variable of the tree are primed
 */
bool is_LHA_dynamics(const tree::node::ptr& p);

bool is_Ode_dynamics(const tree::node::ptr& p);

bool is_affine_tree(const tree::node::ptr& p);

template<typename scalar_type> bool is_affine_expression(
		const tree::node::ptr& p) {
	if (!is_affine_tree(p))
		return false;

	if (!dynamic_cast<valuation_functions::arithmetic_node*> (p.get()))
		if (!dynamic_cast<valuation_functions::const_node<scalar_type>*> (p.get()))
			if (!dynamic_cast<valuation_functions::variable_node*> (p.get()))
				return false;
	return true;
}

/** Cheap check whether the predicate is false */
bool is_false_cheap(const tree::node::ptr& p);

/**
 * count the coefficient of each variable in a tree and stock them in a map <variable_id, Rational>.
 * count the constant too
 * the map must be initialized
 */
template<typename scalar_type> class variable_evaluator {
public:

	tree::node::ptr my_tree;
	valuation_functions::variable_valuation<scalar_type> my_coeffs;
	scalar_type my_const;

	variable_evaluator(
			const tree::node::ptr& p,
			const valuation_functions::variable_valuation<scalar_type>& coefficients,
			const scalar_type& constant) {
		my_tree = p;
		my_coeffs = coefficients;
		my_const = constant;

		typedef typename std::map<variable_id, scalar_type>::const_iterator
				const_map_iterator;
		const_map_iterator it;

		if (valuation_functions::arithmetic_node * a
				= dynamic_cast<valuation_functions::arithmetic_node*> (p.get())) {
			variable_evaluator<scalar_type> child1(a->child1, coefficients,
					constant);
			if (a->my_op == NEG) {
				for (it = child1.my_coeffs.get_valuation_map().begin(); it
						!= child1.my_coeffs.get_valuation_map().end(); ++it) {
					my_coeffs.add(it->first, -it->second);
				}
				my_const = -child1.my_const;
				my_tree = tree::node::ptr(
						new valuation_functions::arithmetic_node(NEG,
								child1.my_tree));
			} else {
				variable_evaluator<scalar_type> child2(a->child2, coefficients,
						constant);
				if (a->my_op == ADD) {
					for (it = child1.my_coeffs.get_valuation_map().begin(); it
							!= child1.my_coeffs.get_valuation_map().end(); ++it) {
						my_coeffs.add(it->first, child1.my_coeffs.get_value(
								it->first) + child2.my_coeffs.get_value(
								it->first));
					}
					my_const = child1.my_const + child2.my_const;
				} else if (a->my_op == SUB) {
					for (it = child1.my_coeffs.get_valuation_map().begin(); it
							!= child1.my_coeffs.get_valuation_map().end(); ++it) {
						my_coeffs.add(it->first, child1.my_coeffs.get_value(
								it->first) - child2.my_coeffs.get_value(
								it->first));
					}
					my_const = child1.my_const - child2.my_const;
				} else if (a->my_op == MUL) {
					if (has_variable(child1.my_tree)) {
						for (it = child1.my_coeffs.get_valuation_map().begin(); it
								!= child1.my_coeffs.get_valuation_map().end(); ++it) {
							my_coeffs.add(it->first, it->second
									* child2.my_const);
						}
						my_const = child1.my_const * child2.my_const;
					} else {
						for (it = child2.my_coeffs.get_valuation_map().begin(); it
								!= child2.my_coeffs.get_valuation_map().end(); ++it) {
							my_coeffs.add(it->first, it->second
									* child1.my_const);
						}
						my_const = child1.my_const * child2.my_const;
					}
				} else if (a->my_op == DIV) {
					for (it = child1.my_coeffs.get_valuation_map().begin(); it
							!= child1.my_coeffs.get_valuation_map().end(); ++it) {
						my_coeffs.add(it->first, it->second / child2.my_const);
					}
					my_const = child1.my_const / child2.my_const;
				}
				my_tree = tree::node::ptr(
						new valuation_functions::arithmetic_node(a->my_op,
								child1.my_tree, child2.my_tree));
			}
		} else if (valuation_functions::variable_node * v
				= dynamic_cast<valuation_functions::variable_node*> (p.get())) {
			my_coeffs.add(v->my_id, scalar_type(1));
		} else if (valuation_functions::const_node<scalar_type> * c
				= dynamic_cast<valuation_functions::const_node<scalar_type>*> (p.get())) {
			my_const = convert_element<scalar_type> (c->my_val);
		} else {
			throw std::runtime_error(" Tree is not a arithmetic expression .");
		}
	}
	;
};

/**
 * Initialize a std::map<variable_id, Rational> with class index_var then count coefficient
 * of all variable of the tree with function count_var
 */
template<typename scalar_type> std::pair<
		valuation_functions::variable_valuation<scalar_type>, scalar_type> calc_aff_coef(
		const tree::node::ptr& p) {
	if (!is_affine_expression<scalar_type> (p))
		throw std::runtime_error("Tree is not a affine expression.");

	variable_id_set ids = valuation_functions::get_variable_ids(p);
	valuation_functions::variable_valuation<scalar_type> coeffs;
	for (variable_id_set::const_iterator it = ids.begin(); it != ids.end(); ++it)
		coeffs.add(*it, scalar_type(0));

	variable_evaluator<scalar_type> v(p, coeffs, scalar_type(0));
	//valuation_functions::variable_valuation<scalar_type> valut()
	std::pair<valuation_functions::variable_valuation<scalar_type>, scalar_type>
			new_pair(std::make_pair<valuation_functions::variable_valuation<
					scalar_type>, scalar_type>(v.my_coeffs, v.my_const));

	return new_pair;
}

/**
 * Calcul coefficients of each variables in p with the function calc_aff_coef(tree::node::ptr p).
 * Next, put this coefficients in matrix A and constants in vector b
 */
template<typename scalar_type> void insert_affine_coeff(math::vdom_matrix<
		scalar_type>& A, math::vdom_vector<scalar_type>& b, variable_id vid,
		const tree::node::ptr& p) {
	//create a std::map<variable_id, Rational> with coefficient of each variables of p,
	//with the function valuation_functions::calc_aff_coef (in parser/assignment_traitment.h)

	std::pair<valuation_functions::variable_valuation<scalar_type>, scalar_type>
			coefs_n_const(calc_aff_coef<scalar_type> (p));
	valuation_functions::variable_valuation<scalar_type> coefficients(
			coefs_n_const.first);

	unsigned int row_index = A.codomain().pos(variable(variable::get_primed_id(vid,0)));
	for (typename std::map<variable_id, scalar_type>::const_iterator it =
			coefficients.get_valuation_map().begin(); it
			!= coefficients.get_valuation_map().end(); ++it) {
		unsigned int column_index = A.domain().pos(variable(it->first));
		A(row_index, column_index) = it->second;
	}
	b[row_index] = coefs_n_const.second;
}

/**
 * Calcul coefficients of each variables in p with the function calc_aff_coef(tree::node::ptr p).
 * Next, put this coefficients in vector a and constants in b, ordered by iimap
 */
template<typename scalar_type> void insert_linear_affine_coeff(
		math::vector<scalar_type>& a, scalar_type& b,
		const index_to_variable_id_map_ptr& iimap, const tree::node::ptr& p) {

	//create a std::map<variable_id, Rational> with coefficient of each variables of p,
	//with the function valuation_functions::calc_aff_coef (in parser/assignment_traitment.h)
	std::pair<valuation_functions::variable_valuation<scalar_type>, scalar_type>
			coefs_n_const(calc_aff_coef<scalar_type> (p));
	valuation_functions::variable_valuation<scalar_type> coefficients(
			coefs_n_const.first);

	for (typename std::map<variable_id, scalar_type>::const_iterator it =
			coefficients.get_valuation_map().begin(); it
			!= coefficients.get_valuation_map().end(); ++it) {
		a[iimap->get_index(it->first)] = it->second;
	}
	b = coefs_n_const.second;
}

}
#endif /*ASSIGNMENT_TRAITMENT_H_*/
