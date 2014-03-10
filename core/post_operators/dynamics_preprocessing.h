/*
 * dynamics_preprocessing.h
 *
 *  Created on: Jun 19, 2013
 *      Author: notroot
 */

#ifndef DYNAMICS_PREPROCESSING_H_
#define DYNAMICS_PREPROCESSING_H_

#include "utility/logger_stopwatch.h"
#include "core/continuous/support_function_provider.h"

#include "math/vdom/vdom_matrix.h"
#include "core/continuous/support_function/sf_base/sf_unary.h"
#include "core/continuous/support_function_provider_utility.h"

// eliminating variables
#include "math/vdom/alg_var_elim.h"

/** Obtain the constraints projected to the state variables given by state_dom */
template<typename scalar_type>
math::lin_constraint_system<double> compute_state_variable_constraints(
		const positional_vdomain& state_dom,
		math::lin_constraint_system<double> cons) {

	// erase constraints on nonstate variables
	variable_id_set excess_variables = cons.get_variable_ids();
	set_difference_assign(excess_variables, state_dom.get_variable_ids());
	//std::cout << "state vars " << state_dom << ", excess: ";
	//print_variable_id_set(std::cout,excess_variables);

	if (!excess_variables.empty()) {
		typedef continuous::constr_polyhedron<scalar_type> poly_type;
		poly_type poly;
		poly.add_constraints(cons);

		// Use existential quantification to eliminate variables
		LOGGER(DEBUG5, __FUNCTION__,
				"quantifying over "+to_string(excess_variables.size())+" non-state variable(s)");

		//std::cout << "poly: " << poly << std::endl;
		poly.existentially_quantify_variables(excess_variables);
		//std::cout << "poly exist: " << poly << std::endl;
		poly.remove_redundant_constraints();
		//std::cout << "poly red: " << poly << std::endl;

		cons = *poly.get_constraints();
	}
	LOGGER_OS(DEBUG7,__FUNCTION__) << "State variable constraints:" << cons;
	return cons;
}
;

/** This file contains routines for bringing
 *  dynamics from the form
 *
 *  x' == Ax + Bu + b
 *  Cx + Du + Ev >= 0,
 *
 *  onto the form
 *
 *  x' == Ax + b + u,
 *  u \in U.
 *
 */

/** Convert dynamics and invariant to nondeterministic ODE. */
template<typename scalar_type>
continuous::ode_affine_dynamics<scalar_type> convert_to_ODE(const math::affine_map<scalar_type>& M_orig,
		const typename continuous::polyhedron<scalar_type>::const_ptr& Inv, typename continuous::polyhedron<scalar_type>::ptr& new_Inv) {
	using namespace continuous;
	using namespace math;

	support_function_provider::const_ptr U;
	affine_map<scalar_type> M = M_orig;
	// divide the dynamics into state equation and additive input disturbance,
	// given the invariant (over both state and input sets)

	lin_constraint_system<scalar_type> con_sys, red_con_sys;
	if (Inv) {
		con_sys = *Inv->get_constraints();
	}

	compute_centered_input_set(M, con_sys, U, red_con_sys);

	// return the new invariant
	new_Inv = typename continuous::polyhedron<scalar_type>::ptr();
	if (Inv) {
		new_Inv = typename continuous::polyhedron<scalar_type>::ptr(Inv->create_universe());

		// remove constraints not over the state variables
		positional_vdomain xdom =  M.codomain();
		red_con_sys = compute_state_variable_constraints<scalar_type>(xdom,red_con_sys);
		new_Inv->add_constraints(red_con_sys);
		new_Inv->remove_redundant_constraints();
	}

	return continuous::ode_affine_dynamics<scalar_type>(M,U);
}

/** Convert the given dynamics and invariant into dynamics and centered inputs
 *
 * Takes M = { x' = A1 x + A2 v + b } and invariant {(x,v) in V} and returns dynamics of the form
 * M = { x' = A1 x + b' + u }, u in U = { A2 V - c }
 * where b' = b + c, c = center(A2 V).
 *
 * The resulting U may be null.
 * The constraints that remain in charge of the invariant are returned in remaining_cons.
 */
template<typename scalar_type>
void compute_centered_input_set(math::affine_map<scalar_type>& M,
		const typename math::lin_constraint_system<scalar_type>& con_sys,
		continuous::support_function_provider::const_ptr& U_centered, math::lin_constraint_system<scalar_type>& remaining_cons) {
	using namespace continuous;
	using namespace math;
	typedef affine_map<scalar_type> affine_map;

	/** @attention The state variables are given by the codomain of A!
	 * The domain contains also the input variables.
	 */
	positional_vdomain dom = M.codomain();

	// start with null = no inputs
	U_centered = support_function_provider::const_ptr();

	/** First, split the domain into state and input variables.
	 * Put the input variables into U. */
	bool use_old_cons = true;
	if (M.domain() != M.codomain()) {
		// split A into square A1 and A2
		affine_map Mstate;
		affine_map Minput;
		separate_states_from_inputs(M, Mstate, Minput);

		M = Mstate;
		LOGGER_OS(DEBUG7, __FUNCTION__) << "split dynamics into "
				<< Mstate << " and " << Minput << " with input domain "
				<< Minput.domain();

		// Define the set of inputs using the invariant.
		// This is the original invariant with both state and nonstate variables
		if (Minput.domain().size() > 0) {
			if (con_sys.size() > 0) {
				use_old_cons = false;

				// try to eliminate algebraic variables using the invariant
				// state variables
				positional_vdomain xdom = Mstate.codomain();
				// algebraic variables
				positional_vdomain udom = Minput.domain();
				// the rest of the variables
				positional_vdomain vdom(con_sys.get_variable_ids());
				vdom.remove_variables(xdom.get_variables());
				vdom.remove_variables(udom.get_variables());

				// try to simplify the constraints

				LOGGER_OS(DEBUG6, __FUNCTION__)
						<< "Eliminating algebraic variables using constraints "
						<< con_sys;

				alg_var_elim<double> my_eliminator = alg_var_elim<double>(xdom,
						udom, vdom, Mstate.get_A().get_matrix(),
						Minput.get_A().get_matrix(),
						Mstate.get_b().get_vector(), con_sys);
				matrix<scalar_type> A_elim, B_elim;
				vector<scalar_type> b0_elim;
				//my_eliminator.print();

				positional_vdomain B_dom;
				my_eliminator.eliminate(A_elim, B_elim, B_dom, b0_elim, remaining_cons);

//				std::cout << "A elim:" << A_elim << "B elim:" << B_elim
//						<< "b0 elim:" << b0_elim;

				Mstate = affine_map(xdom, A_elim, b0_elim);
				M = Mstate;
				// assume no inputs are needed (B_elim is zero)
				Minput = affine_map();

				if (!B_elim.is_zero()) {
					//std::cout << "centering inputs" << std::endl;
					vdom_matrix<scalar_type> B(xdom, B_dom, B_elim);
					Minput = affine_map(B);

					// here we could filter out constraints not on U,
					// but probably not worth the effort
					typename constr_polyhedron<scalar_type>::ptr V(
							new constr_polyhedron<scalar_type>());
					V->add_constraints(remaining_cons);
					typename support_function::sf_unary<scalar_type>::ptr U(
							new support_function::sf_unary<scalar_type>(V,
									Minput));

					vdom_vector<scalar_type> b_new = Mstate.get_b();
					vdom_vector<scalar_type> c(xdom);
					bool is_point = is_input_set_point<scalar_type>(U,
							U_centered, c);
					if (is_point) {
						// we don't need U_centered
						U_centered = support_function_provider::const_ptr();
					} else {
						// remove the constraints not related to the inputs
						lin_constraint_system<scalar_type> U_cons =
								compute_state_variable_constraints<scalar_type>(B_dom,
										remaining_cons);
						V = typename constr_polyhedron<scalar_type>::ptr(
								new constr_polyhedron<scalar_type>());
						V->add_constraints(U_cons);

						// translate the constraints by -c
						vdom_vector<scalar_type> c_U =
								compute_center_of_bounding_box<scalar_type>(*V);
						math::affine_map<scalar_type> translate_map(B_dom,
								-c_U);
//						std::cout << "c_U: " << c_U << ", domain: " << con_dom
//								<< ", translating by " << translate_map
//								<< std::endl;
						typename constr_polyhedron<scalar_type>::ptr V_centered =
								typename constr_polyhedron<scalar_type>::ptr(
										new constr_polyhedron<scalar_type>(
												apply_map(*V, translate_map)));
						U_centered = typename support_function::sf_unary<
								scalar_type>::ptr(
								new support_function::sf_unary<scalar_type>(V_centered,
										Minput));
					}
					M = affine_map(Mstate.get_A(), Mstate.get_b() + c);
					//std::cout << "b = " << Mstate.get_b() << ", c = " << c << ", bfinal = " << M.get_b() << std::endl;
				}
			} else {
				throw basic_exception(
						"Unbounded algebraic variables "
								+ logger::formatted_to_string(
										Minput.domain().get_variables()));
			}
		}
	}
	if (use_old_cons) {
		// no need to change constraints
		remaining_cons = con_sys;
	}

	// report eigenvalues
	IFLOGGER(NEVER)
	{
		vector<scalar_type> v_real_d, v_imag_d;
		matrix<scalar_type> V_d;
		compute_eigenvalues(M.get_A().get_matrix(), v_real_d, v_imag_d, V_d);
		LOGGER_OS(DEBUG4, __FUNCTION__) << "real parts of Eigenvalues:"
				<< v_real_d << ", imag parts of Eigenvalues:" << v_imag_d;
		LOGGER_OS(DEBUG4, __FUNCTION__) << " Eigenvectors:" << V_d;
	}
}

#endif /* DYNAMICS_PREPROCESSING_H_ */
