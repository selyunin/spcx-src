#ifndef CONVERT_PREDICATE_H_
#define CONVERT_PREDICATE_H_

//#include "ppl_polyhedron/convert_to_ppl.h"
#include "core/predicates/valuation_function_tree_nodes.h"
#include "core/hybrid_automata/location_constraint_set.h"
#include "math/vdom/lin_constraint.h"
#include "core/continuous/polyhedra/constr_polyhedron/constr_polyhedron.h"
#include "core/continuous/polyhedra/constr_polyhedron/constr_polyhedron_operators.h"
#include "core/hybrid_automata/assignment_traitment.h"

#include "core/continuous/predicate_continuous_set_constructors.h"

/** Forward declarations */
namespace discrete {
class discrete_set;
typedef boost::shared_ptr<discrete_set> discrete_set_ptr;
}
namespace continuous {
class continuous_set;
typedef boost::shared_ptr<continuous_set> continuous_set_ptr;
class continuous_set_transform;
typedef boost::shared_ptr<continuous_set_transform> continuous_set_transform_ptr;
class continuous_dynamics;
typedef boost::shared_ptr<continuous_dynamics> continuous_dynamics_ptr;
}
namespace hybrid_automata {
class symbolic_state_collection;
}

namespace continuous {

/** Converts a predicate to an affine map.
 *
 * Returns true if the conversion is successful,
 * and false otherwise.
 *
 * If include_codomain_in_domain=true, the variables
 * in the codomain are included in the domain, even
 * if they are not used.
 *
 * @todo check if every variable is assigned only once
 */
template<typename scalar_type>
bool convert_predicate_to_affine_map(tree::node::ptr p, math::affine_map<
		scalar_type>& M, bool include_codomain_in_domain = false) {
	if (assignment_traitment::is_affine_tree(p)
			&& assignment_traitment::is_Ode_dynamics(p)) {
		std::map<variable_id, tree::node::ptr> all_Ode(
				assignment_traitment::get_all_Ode(p));

		using namespace valuation_functions;
		// Tho codomain (target) are the primed variables
		positional_vdomain codom(get_primed_variable_ids(p, 1));
		// The domain are the unprimed variables
		positional_vdomain dom(get_primed_variable_ids(p, 0));

		if (include_codomain_in_domain) {
			// to avoid unnecessary remapping, create the dom from the codom
			// and then add additional variables.
			// that way, if there are no additional variables then dom=codom
			std::set<variable> vset=dom.get_variables();
			dom=codom;
			dom.add_variables(vset);
//			std::cout << "codomain extended with domain: " << codom << std::endl;
		}

		math::vdom_matrix<scalar_type> A(codom, dom);
		math::vdom_vector<scalar_type> b(codom);

		for (std::map<variable_id, tree::node::ptr>::const_iterator it =
				all_Ode.begin(); it != all_Ode.end(); ++it) {
			// for each ode in all_Ode, insert coefficients in A and constant in b
			variable_id target_var = it->first;
			assignment_traitment::insert_affine_coeff<scalar_type>(A, b,
					target_var, it->second);
		}
		if (!all_Ode.empty()) {
			M = math::affine_map<scalar_type>(A, b);
		} else {
			// check whether the predicate is void (false) or universe
			if ( assignment_traitment::is_false_cheap(p) ) {
				M = math::affine_map<scalar_type>::void_map();
			} else {
				M = math::affine_map<scalar_type>::universe_map();
			}
		}

		//std::cout << "converted " << p << " to " << M << std::endl;
		return true;
	} else {
		return false;
	}
};

/** Converts a predicate to a vector of ODE functions
 *
 * Deterministic ODE (Ordinary Differential Equation) system, i.e.,
 *  \f$ \dot x_i = f_i(x_i) \f$ over some \f$i\f$.
 *
 * The codomain is the domain of the derivatives \dot x_i to which
 * the function funcs[i] is attributed.
 *
 * Note that x_i are the unprimed variables.
 *
 * Returns true if the conversion is successful,
 * and false otherwise.
 *
 * @todo check if every variable is assigned only once
 */
template<typename scalar_type>
bool convert_predicate_to_function_vector(tree::node::ptr p, positional_vdomain& codom,
		std::vector<tree::node::ptr>& funcs) {
	// check if the predicate is in ODE form (only equalities with exactly one primed
	// variable on each side
	if (assignment_traitment::is_Ode_dynamics(p)) {
		std::map<variable_id, tree::node::ptr> all_Ode(
				assignment_traitment::get_all_Ode(p));

		using namespace valuation_functions;
		// Tho codomain (target) are the primed variables
		codom=positional_vdomain(get_primed_variable_ids(p, 1));

		// The domain are the unprimed variables
		positional_vdomain dom(get_primed_variable_ids(p, 0));

		// assign the functions
		funcs = std::vector<tree::node::ptr>(codom.size());
		for (std::map<variable_id, tree::node::ptr>::const_iterator it =
				all_Ode.begin(); it != all_Ode.end(); ++it) {
			// for each ode in all_Ode, insert coefficients in A and constant in b
			// remove prime from variable
			variable_id target_var = variable::get_primed_id(it->first,0);
			funcs[codom.pos(variable(target_var))] = it->second;
		}

		//std::cout << "converted " << p << " to " << M << std::endl;
		return true;
	} else {
		return false;
	}
};

template<typename scalar_type> math::lin_expression<Rational> convert_affine_tree_to_lin_expression(
		const tree::node::ptr& p) {
	if (!assignment_traitment::is_affine_expression<scalar_type>(p))
		throw std::runtime_error("Tree is not a affine expression.");

	index_to_variable_id_map_ptr iimap(index_to_variable_id_map::empty_map());
	iimap = iimap->get_map_with_ids_added(valuation_functions::get_variable_ids(p));
	Rational default_value = Rational(0);
	math::rational_vector a(iimap->dimensions(), default_value);
	Rational b = default_value;
	// insert coefficients in A and constant in b
	assignment_traitment::insert_linear_affine_coeff<scalar_type>(a, b, iimap, p);
	math::lin_expression<Rational> l(a, b, iimap);

	return l;
}

template<typename scalar_type> math::lin_constraint<Rational> convert_affine_predicate_to_lin_constraint(
		const tree::node::ptr& p) {
	if (!assignment_traitment::is_affine_tree(p))
		throw std::runtime_error("Tree is not a affine predicate.");

	if (valuation_functions::comparison_node * q
			= dynamic_cast<valuation_functions::comparison_node*> (p.get())) {
		math::lin_constraint<Rational> l(
				convert_affine_tree_to_lin_expression<scalar_type> (q->child1)
						- convert_affine_tree_to_lin_expression<scalar_type> (q->child2), q->my_op);
		return l;
	} else
		throw std::runtime_error("Tree is not a affine predicate.");
}

template<typename scalar_type>
continuous::constr_polyhedron<Rational> convert_affine_predicates_to_constr_polyhedron(
		const tree::node::ptr& p) {

	continuous::constr_polyhedron<Rational> poly;
	if (assignment_traitment::is_affine_tree(p)) {
		std::vector<tree::node::ptr> all_constr(assignment_traitment::get_all_constraint( p));

		for (std::vector<tree::node::ptr>::size_type i = 0; i < all_constr.size(); ++i) {
			poly.add_constraint(convert_affine_predicate_to_lin_constraint<scalar_type> (
					all_constr[i]));
		}
	} else
		throw std::runtime_error("Tree is not a affine predicates.");
	return poly;
}

template<typename scalar_type>
valuation_functions::variable_valuation<math::lin_expression<Rational> > create_lin_expression_valuation(
		const tree::node::ptr& p) {
	if (assignment_traitment::is_affine_tree(p) && assignment_traitment::is_Ode_dynamics(p)
			&& assignment_traitment::all_variable_define(p)) {

		valuation_functions::variable_valuation<math::lin_expression<Rational> > v;
		std::map<variable_id, tree::node::ptr> all_Ode(assignment_traitment::get_all_Ode( p));
		for (std::map<variable_id, tree::node::ptr>::const_iterator it = all_Ode.begin(); it
				!= all_Ode.end(); ++it) {
			// for each ode in all_Ode, convert tree to lin_expression, then insert it in the valuation
			math::lin_expression<Rational> le = convert_affine_tree_to_lin_expression<
					scalar_type> (it->second);
			variable_id vid = variable::get_id_primedness_decreased(it->first);
			v.add(variable(vid).get_name(), le);
		}
		return v;
	} else
		throw std::runtime_error("Tree is not ode affine.");
}

}

namespace hybrid_automata {

/** Convert a predicate to a set of symbolic states and adds it to an existing
 * collection.
 */
void convert_and_add_predicate_to_symbolic_state_collection(const tree::node::ptr& p,
		hybrid_automata::symbolic_state_collection& sstates);

}

#endif /*CONVERT_PREDICATE_H_*/
