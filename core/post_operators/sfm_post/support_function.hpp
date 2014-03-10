#include "support_function.h"

#include <sstream>
#include "utility/basic_exception.h"
#include "utility/basic_warning.h"
#include "utility/logger_stopwatch.h"

#include "math/numeric/container_comp.h"
//#include "core/continuous/support_function/sf_unary.h"
//#include "core/continuous/support_function/sfm/sfm_cont_set.h"
//#include "core/continuous/support_function/sf_ode_deterministic.hpp"
#include "core/continuous/support_function_provider_utility.h"
#include "core/continuous/polyhedra/hyperbox/bounding_box_utility.h"
#include "sf_evaluator.h"

namespace continuous {
namespace support_function {

template<class scalar_type> typename sfm_cont_set<scalar_type>::ptr support_function_time_elapse(
		const support_function_provider::const_ptr& csup,
		const math::affine_map<scalar_type>& M_orig, support_function_provider::const_ptr U, const typename polyhedron<
				scalar_type>::const_ptr& inv, double time_horizon,
		double time_step, double step_tolerance) {

	if (time_horizon < 0.0) {
		throw std::runtime_error(
				"support_function_time_elapse : cannot handle negative time_horizon "
						+ to_string(time_horizon));
	}
	if (time_step <= 0.0) {
		throw std::runtime_error(
				"support_function_time_elapse : cannot handle negative time_step "
						+ to_string(time_step));
	}

	if (csup->is_empty()) {
		typename sfm_cont_set<scalar_type>::ptr empty_set(sfm_cont_set<
				scalar_type>::empty_set());
		return empty_set;
	}

	// call implementation of iterative computation of time elapse set using support function.

	//	typedef Rational scalar_type;
	typedef math::matrix<scalar_type> matrix_type;
	typedef math::vector<scalar_type> vector_type;

	//	dimension_t dim = csup->get_dim();

	// If there are const variables, add their initial values to the inputs
	typename polyhedron<scalar_type>::const_ptr new_inv=inv;

	math::affine_map<scalar_type> M = M_orig;
	// If the domain and codomain are not the same,
	// we need to reorder and seperate state from input variables
	if (M.domain() != M.codomain()) {
		// split A into square A1 and A2
		math::affine_map<scalar_type> Mstate;
		math::affine_map<scalar_type> Minput;
		math::separate_states_from_inputs(M, Mstate, Minput);

		M = Mstate;
		LOGGER_OS(DEBUG5,"support_function_time_elapse") << "split dynamics into " << Mstate << " and " << Minput << " with input domain " << Minput.domain();

		// Define the set of inputs using the invariant.
		// This is the original invariant with both state and nonstate variables
		if (Minput.domain().size() > 0) {
			if (U) {
				// don't change U
			} else if (new_inv) {
				U = new_inv;
			} else {
				// make a universe set for this domain
				U	= support_function_provider::ptr(new hyperbox<
								scalar_type> ());
			}
			U = support_function_provider::ptr(new support_function::sf_unary<
					scalar_type>(U, Minput));
		}
	}
	matrix_type A = M.get_A().get_matrix();
	math::vdom_vector<scalar_type> b_named = M.get_b();
	positional_vdomain dom = M.codomain();
	index_to_variable_id_map_ptr iimap = dom.get_index_to_variable_id_map(); // iimap that tells which variables are at which index
	assert(b_named.domain().get_index_to_variable_id_map() == iimap);

	std::list<vector_type> directions;
	choose_directions(dom, directions); // R is the number of directions.

//		std::cout << "supp time elapse with dynamics A=" << A << ", b=" << b_named << " and U={" << U
//				<< "}" << " domain " << dom << std::endl;

	matrix_type symbolic_reach_set;
	//cout << csup << endl;
	//cout << U << endl;

	// Add the inv directions to the directions_list. Also construct a vector of the inhomogenous terms of the invariants.
	// The order of the inhomogenous terms in the vector should be the same as the corresponding order of directions in the directions list.
	// For directions which are not invariant directions, the corresponding inhomogenous term is set to positive infinity.

	typename math::lin_constraint_system<scalar_type>::const_ptr
			inv_constr_ptr = new_inv->get_constraints();

	std::list<scalar_with_infinity<scalar_type> > inh_terms;
	scalar_with_infinity<scalar_type> my_scalar;

	/* This code below checks if a direction is one of the negative invariant directions. If it is not, then pos_infiity
	 * is stored for its corresponding entry in the inh_coeffs list. If the direction is one of the invariant
	 * direction, then the inhomogenous term of the invariant constraint is stored in the corresponding entry
	 * of inh_coeff list.
	 */
	for (typename std::list<vector_type>::const_iterator it =
			directions.begin(); it != directions.end(); it++) {
		bool flag = false;
		for (typename math::lin_constraint_system<scalar_type>::const_iterator
				iter = inv_constr_ptr->begin(); iter != inv_constr_ptr->end(); ++iter) {
			math::vdom_vector<scalar_type> rvec =
					iter->get_canonic_l().get_vdom_vec();
			// if the direction has all the rvec's variables
			if (rvec.remap(dom)) {
				vector_type inv_direction = rvec.get_vector();

				//			std::cout << "inv_direction:" << inv_direction << std::endl;
				//			std::cout << "*it:" << *it << std::endl;

				if (math::numeric::is_MEQ(*it, -inv_direction)) {
					inh_terms.push_back(
							-scalar_with_infinity<scalar_type> (
									iter->get_canonic_inh_coeff()));
					flag = true; // The assumption is the inv set contains no redundant constraints
					break;
				}
			}
		}
		if (!flag)
			inh_terms.push_back(my_scalar.pos_infty());
	}

	// Variant with negative inv directions @author Rajarshi
	/*
	 * This code below checks if a negative invariant direction is already present in existing list of directions.
	 * Only the new negative invariant directions are added to the directions list.
	 */
	// GF 2011-12-03 : Also added positive invariant directions (was a side effect of previous algo)
	bool add_positive_inv_directions = true;
	for (typename math::lin_constraint_system<scalar_type>::const_iterator it =
			inv_constr_ptr->begin(); it != inv_constr_ptr->end(); ++it) {
		bool red_flag = false; // neg inv direction already present
		bool green_flag = false; // pos inv direction already present
		math::vdom_vector<scalar_type> rvec =
				it->get_canonic_l().get_vdom_vec();
		// don't use any constraints that are over nonstate variables
		if (rvec.remap(dom)) {
			vector_type inv_direction = rvec.get_vector();
			//		std::cout << "invariant direction is:" << inv_direction << std::endl;

			for (typename std::list<vector_type>::const_iterator iter =
					directions.begin(); iter != directions.end() && (!red_flag || !green_flag); iter++) {
				if (math::numeric::is_MEQ(*iter, -inv_direction)) {
					red_flag = true; // the neg invariant is already present as a direction
				}
				if (add_positive_inv_directions && math::numeric::is_MEQ(*iter, inv_direction)) {
					green_flag = true; // the neg invariant is already present as a direction
				}
			}
			if (!red_flag) {
				directions.push_back(-inv_direction); // Add the non-redundant negative inv direction to the list of directions
				inh_terms.push_back(-scalar_with_infinity<scalar_type> (
						it->get_canonic_inh_coeff()));
			}
			if (add_positive_inv_directions && !green_flag) {
				directions.push_back(inv_direction); // Add the non-redundant negative inv direction to the list of directions
				inh_terms.push_back(my_scalar.pos_infty());
			}
		}
	}

	//	testing
	/*	std::cout << "directions after adding inv directions : "<< std::endl;
	 for(typename std::list<vector_type>::const_iterator it = directions.begin(); it != directions.end(); it++){
	 std::cout << "direction:" << *it << std::endl;
	 }
	 */
	std::vector<scalar_type> delta_vec;
	/*
	 * Since the new inv directions are added to the list of directions at the end with the same order
	 * in which the corresponding inv inh_terms are added to the list inh_coeffs, making the vector is just
	 * a direct elementwise copy.
	 */
	//make a vector from the list
	std::vector<scalar_with_infinity<scalar_type> >
			inh_coeffs(inh_terms.size());
	//	std::cout << "the inh coeffs are: " << std::endl;
	unsigned int i = 0;
	for (typename std::list<scalar_with_infinity<scalar_type> >::const_iterator
			it = inh_terms.begin(); it != inh_terms.end(); ++it, i++) {
		inh_coeffs[i] = *it;
		//		std::cout << *it << std::endl;
	}
	symbolic_reach_set = sf_evaluator<scalar_type> (A, b_named.get_vector(),
			iimap, csup, U).compute_matrix_tol(directions, inh_coeffs,
			scalar_type(time_step), time_horizon, delta_vec, step_tolerance);
	//cout << "Reach Set Matrix: " << symbolic_reach_set;

	/**
	 * Create a list of intervals Ij such that each omega
	 */
	// create the postc  parameter structure .
	postc_params<scalar_type> my_postc_params;
	my_postc_params.dynamics_A = A;
	my_postc_params.dynamics_b = b_named.get_vector();
	my_postc_params.initial_set_ptr = csup;
	my_postc_params.invariant_set_ptr = new_inv;
	my_postc_params.input_set_ptr = U;
	my_postc_params.delta = scalar_type(time_step);
	my_postc_params.time_horizon = scalar_type(time_horizon);
	my_postc_params.delta_vec = delta_vec;

	typename sfm_cont_set<scalar_type>::direction_store dir_store = convert_to_direction_store<scalar_type>(directions);
	typename sfm_cont_set<scalar_type>::ptr sfm_ptr = typename sfm_cont_set<
			scalar_type>::ptr(new sfm_cont_set<scalar_type> (my_postc_params,
			symbolic_reach_set, dir_store, iimap));

	// @todo What is this for?
	//sfm_ptr->embed_variables(csup->get_variable_ids());
	//testing..
	//cout << "IIMAP :" << iimap << std::endl;

	return sfm_ptr;
}
;

}
}

