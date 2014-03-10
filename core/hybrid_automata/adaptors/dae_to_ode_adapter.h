/*
 * dae_to_ode_adapter.h
 *
 *  Created on: Jun 24, 2013
 *      Author: frehse
 */

#ifndef DAE_TO_ODE_ADAPTER_H_
#define DAE_TO_ODE_ADAPTER_H_

#include "core/hybrid_automata/location.h"
#include "core/hybrid_automata/transition.h"
#include "core/hybrid_automata/hybrid_automaton_on_the_fly_adapter.h"
#include "core/post_operators/dynamics_preprocessing.h"
#include "core/continuous/continuous_dynamics/ode_affine_dynamics.h"
#include "core/continuous/support_function/sf_derived/sf_sum.h"

// for transitions
#include "core/continuous/polyhedra/constr_polyhedron/constr_polyhedron_operators.h"


namespace hybrid_automata {

/** Converts a location from DAE to ODE form */
template<typename scalar_type>
class dae_to_ode_location_adapter : public location_adapter {
public:
	/** Default behavior is to do return the location without any modification. */
	location_ptr accept(hybrid_automaton::ptr aut, location_ptr l) {
		LOGGERSWOC(DEBUG6,__FUNCTION__,"DAE to ODE conversion for locations");

		using namespace continuous;

		// get the data structures from l:
		const time_constraints& tcons = l->get_time_constraints();
		// dynamics, invariant
		typename polyhedron<scalar_type>::const_ptr inv_poly_ptr;
		if (tcons.get_invariant()) {
			inv_poly_ptr = boost::dynamic_pointer_cast<
					const polyhedron<scalar_type> >(tcons.get_invariant());
			if (!inv_poly_ptr) {
				std::string tname = typeid(*tcons.get_invariant().get()).name();
				throw basic_exception(
						std::string(__FUNCTION__) + " cannot handle invariant of type "
								+ tname);
			}
		}

		typedef ode_affine_dynamics<scalar_type> ode_aff_dyn;
		typename ode_aff_dyn::const_ptr dp = boost::dynamic_pointer_cast<
				const ode_aff_dyn>(tcons.get_dynamics());
		typename typed_dynamics<scalar_type>::const_ptr dp_nonlin =
				boost::dynamic_pointer_cast<const typed_dynamics<scalar_type> >(
						tcons.get_dynamics());
		if (dp) {
			// first, we treat degenerated cases
			if ((dp && dp->is_universe())) {
				// there is no restriction on the dynamics,
				/** Save original invariant */
				l->add_tag("INV",inv_poly_ptr);
				l->add_tag("STATE_DOM",dp->codomain());

			} else if ((dp && dp->is_zero())) {
				// there is no change in the variables,
				// since the derivatives are identical zero.
				/** Save original invariant */
				l->add_tag("INV",inv_poly_ptr);
				l->add_tag("STATE_DOM",dp->codomain());

			} else {
				if (dp->domain() != dp->codomain()) {
					typedef typename ode_aff_dyn::offset_type offset_type;
					typedef typename offset_type::const_ptr offset_const_ptr;

					// keep the original offset
					offset_const_ptr previous_U = dp->get_U();
					const math::affine_map<scalar_type>& M_orig(*dp);
					typename polyhedron<scalar_type>::ptr new_Inv;
					ode_affine_dynamics<scalar_type> ode_dyn = convert_to_ODE(
							M_orig, inv_poly_ptr, new_Inv);

					typename ode_affine_dynamics<scalar_type>::ptr ode_ptr(
							new ode_affine_dynamics<scalar_type>(ode_dyn));

					// add previous offset to new offset
					if (previous_U) {
						offset_const_ptr new_U = ode_dyn.get_U();
						offset_const_ptr sum_U(
								new support_function::sf_sum<scalar_type>(previous_U,
										new_U));
						ode_ptr->set_U(sum_U);
					}

					/** Save original invariant */
					l->add_tag("INV",inv_poly_ptr);
					l->add_tag("STATE_DOM",ode_ptr->codomain());

					time_constraints new_tcon(new_Inv, ode_ptr);
					l->set_time_constraints(new_tcon);

					LOGGER_OS(DEBUG7,__FUNCTION__) << "took location with inv " << inv_poly_ptr
							<< " and dynamics " << *dp << std::endl;
					LOGGER_OS(DEBUG6,__FUNCTION__) << "produced location with inv " << new_Inv
							<< " and dynamics " << ode_ptr << std::endl;
					LOGGER_OS(DEBUG6,__FUNCTION__) << "state variables " << l->get_tag<positional_vdomain>("STATE_DOM");
				} else {
					if (!l->has_tag("STATE_DOM")) {
						l->add_tag("STATE_DOM", dp->codomain());
					}
				}
			}
		}
		return l;
	}
	;
};

/** Transform transitions to added input form
 *
 * Adds the pre-image of the target invariant to the guard.
 * Brings affine transforms to the form
 *    x'=Ax+b+U
 *
 */
template<typename scalar_type>
class dae_to_ode_transition_adapter : public transition_adapter {
public:
	typedef continuous::polyhedron<scalar_type> poly_type;
	typedef typename poly_type::const_ptr poly_const_ptr;

	typedef continuous::reset_affine_transform<scalar_type>  reset_trans_type;
	typedef boost::shared_ptr<reset_trans_type> reset_trans_ptr;

	struct problem_instance {
		poly_const_ptr sinv_ptr;
		poly_const_ptr tinv_ptr;
		poly_const_ptr sinv_orig_ptr; // contains algebraic constraints
		poly_const_ptr guard_poly_ptr;
		const reset_trans_type* reset_trans_ptr;
		positional_vdomain source_state_dom; // state variables of the source location
		positional_vdomain target_state_dom; // state variables of the target location
	};

	problem_instance get_instance(hybrid_automaton_ptr aut, transition_ptr trans_ptr) const {
		problem_instance inst;

		if (!trans_ptr) return inst;

		const location_id& s_id = trans_ptr->get_source();
		const location_id& t_id = trans_ptr->get_target();
		poly_const_ptr source_inv;
		poly_const_ptr target_inv;
		const jump_constraints& j = trans_ptr->get_jump_constraints();
		try {
			LOGGER_OS(DEBUG5,__FUNCTION__) << "adapting transition with label " <<  named_label::get_name(trans_ptr->get_label()) << " from location " << canonicalize_location_constraint(aut, s_id) << " to location " << canonicalize_location_constraint(aut, t_id);

			location_const_ptr sl = aut->get_location(s_id);
			location_const_ptr tl = aut->get_location(t_id);
			// Get the invariant of the source location
			source_inv = boost::static_pointer_cast<const poly_type>(sl->get_time_constraints().get_invariant());
			// GF-2013-10-11: Don't know where here this used to be the OLD invariant (INV)
			if (sl->has_tag("INV")) {
				inst.sinv_orig_ptr = sl->get_tag<poly_const_ptr>("INV");
			} else {
				inst.sinv_orig_ptr = source_inv;
			}

			// Get the invariant of the target location
			if (false && tl->has_tag("INV")) {
				target_inv = tl->get_tag<poly_const_ptr>("INV");
			} else {
				target_inv = boost::static_pointer_cast<const poly_type>(tl->get_time_constraints().get_invariant());
			}

			if (sl->has_tag("STATE_DOM")) {
				inst.source_state_dom = sl->get_tag<positional_vdomain>("STATE_DOM");
			} else {
				throw basic_exception("missing source state variables");
			}
			if (tl->has_tag("STATE_DOM")) {
				inst.target_state_dom = tl->get_tag<positional_vdomain>("STATE_DOM");
			} else {
				throw basic_exception("missing target state variables");
			}
		} catch (std::exception& e) {
			std::stringstream s;
			logger::copyfmt_to(s);
			s << canonicalize_location_constraint(aut, s_id);
			s << " to location ";
			s << canonicalize_location_constraint(aut, t_id);
			throw basic_exception(
					"could not adapt transition with label "
							+ named_label::get_name(trans_ptr->get_label())
							+ ", from location " + s.str() + ".", e);
		}

		inst.sinv_ptr = source_inv;
		inst.tinv_ptr = target_inv;

		inst.guard_poly_ptr = poly_const_ptr();
		if (j.get_guard()) {
			inst.guard_poly_ptr = boost::dynamic_pointer_cast<
					const continuous::polyhedron<scalar_type> >(j.get_guard());
			if (!inst.guard_poly_ptr) {
				throw std::runtime_error(
						"guard type not supported");
			}
		}
		inst.reset_trans_ptr = 0;
		if (j.get_transform()) {
			inst.reset_trans_ptr =
					dynamic_cast<const continuous::reset_affine_transform<
							scalar_type>*>(j.get_transform().get());
			if (!inst.reset_trans_ptr) {
				std::stringstream ss;
				logger::copyfmt_to(ss);
				ss << "offending transform:" << j.get_transform();
				throw std::runtime_error(
						"discrete_post_stc: transform type not supported\n"
								+ ss.str());
			}
		}
		return inst;
	}
	;

	/** Adapt a transition */
	virtual transition_ptr accept(hybrid_automaton_ptr aut, transition_ptr t) {
		LOGGERSWOC(DEBUG6,__FUNCTION__,"DAE to ODE conversion for transitions");
		using namespace continuous;
		using namespace math;
		using namespace support_function;

		typedef typename polyhedron<scalar_type>::ptr poly_ptr;
		typedef typename constr_polyhedron<scalar_type>::ptr con_poly_ptr;
		typedef lin_constraint_system<scalar_type> cons_type;

		if (!t) {
			return t;
		}
		problem_instance inst = get_instance(aut,t);
		jump_constraints tcons = t->get_jump_constraints();

		/***************************************************************************
		 *  Compute the effective guard set
		 *
		 *  Add (reverse mapped) invariants of source and target location to guard,
		 *  then project to state variables
		 ***************************************************************************/

		// add the pre-transformed invariant constraints to the guard constraint
		con_poly_ptr guard_poly_ptr(new constr_polyhedron<scalar_type>());
		constr_polyhedron<scalar_type>& guard_poly(*guard_poly_ptr);
		if (inst.guard_poly_ptr) {
			// add guard constraints (duh)
			guard_poly.add_constraints(inst.guard_poly_ptr->get_constraints());
		}
		if (inst.sinv_ptr) {
			// add source invariant constraints
			guard_poly.add_constraints(inst.sinv_ptr->get_constraints());
		}
		if (inst.tinv_ptr) {
			constr_polyhedron<scalar_type> exit_set_poly;
			typename constr_polyhedron<scalar_type>::const_ptr tgt_inv_poly_ptr =
					boost::dynamic_pointer_cast<
							const constr_polyhedron<scalar_type> >(
							inst.tinv_ptr);
			if (tgt_inv_poly_ptr) {
				affine_map<scalar_type> M = *inst.reset_trans_ptr;

				support_function_provider::const_ptr U =
						boost::dynamic_pointer_cast<
								const support_function_provider>(
								inst.reset_trans_ptr->get_input_set());

				constr_polyhedron<scalar_type> dummy;

				LOGGER_OS(DEBUG7,__FUNCTION__)
						<< "target invariant: "
						<< tgt_inv_poly_ptr << ", map: "<< M;

				/** reverse map target invariant */
				reverse_map_constraints(tgt_inv_poly_ptr, M, U, exit_set_poly,
						dummy);
				LOGGER_OS(DEBUG7,__FUNCTION__)
						<< "invariant poly reverse mapped to: "
						<< exit_set_poly;
				guard_poly.add_constraints(exit_set_poly.get_constraints());
			} else
				throw std::runtime_error(
						"Target location invariant set type not supported\n");
		}

		/** Compute state variables */
		variable_id_set source_state_variables = inst.source_state_dom.get_variable_ids();
		variable_id_set source_nonstate_variables = guard_poly_ptr->get_variable_ids();
		set_difference_assign(source_nonstate_variables, source_state_variables);

		/**************************************************************
		 * Compute the input set of the transform
		 *
		 * Project the guard to the input set of the transform
		 *
		 * @note By definition, the source invariant should suffice,
		 * but the guard might be more accurate in some cases (when)?
		 * ************************************************************/
		if (inst.reset_trans_ptr) {
			affine_map<scalar_type> M = *inst.reset_trans_ptr;

			// remove nonstate variables
			positional_vdomain dom=M.domain();
			positional_vdomain codom=M.codomain();
			variable_id_set target_state_variables = inst.target_state_dom.get_variable_ids();
			//variable_id_set target_nonstate_variables = inst.tinv_ptr->get_variable_ids();
			variable_id_set target_nonstate_variables = codom.get_variable_ids();
			set_difference_assign(target_nonstate_variables, target_state_variables);

			std::set<variable> target_nonstate = create_variable_set(target_nonstate_variables);
			if (!target_nonstate.empty()) {
				LOGGER(DEBUG7, __FUNCTION__,
						"removing nonstate variables "
								+ logger::formatted_to_string(target_nonstate));
				math::affine_map<scalar_type> M_proj = math::affine_map<
						scalar_type>::projection_map(inst.target_state_dom, M.codomain());
				M = concatenate(M_proj, M);
			}

			if (M.codomain() != inst.target_state_dom) {
				if (target_nonstate.empty()) {
					// we just need to reorder
					M.reorder(inst.target_state_dom,M.domain());
				} else {
					LOGGER_OS(LOW,__FUNCTION__) << "Map = " << M << std::endl;
					LOGGER_OS(LOW,__FUNCTION__) << "nonstate: "
							<< target_nonstate << std::endl << ", codom: "
							<< M.codomain() << ", target: "
							<< inst.target_state_dom;
					throw basic_exception(
							"Codomain of map does not correspond to target states.");
				}
			}

			support_function_provider::const_ptr mapped_U;

			// If the domain and codomain are not the same,
			// we need to reorder and seperate state from input variables
			if (M.domain() != inst.source_state_dom) {
				// split A into square A1 and A2
				affine_map<scalar_type> Mstate;
				affine_map<scalar_type> Minput;
				decompose_on_domain(inst.source_state_dom, M, Mstate, Minput);

				M = Mstate;
				//									std::cout << "split dynamics into " << Mstate << " and " << Minput
				//											<< " with input domain " << Minput.domain() << std::endl;

				// Define the set of inputs using the invariant.
				// This is the original invariant with both state and nonstate variables
				if (Minput.domain().size() > 0 && !Minput.is_zero()) {
					LOGGER(DEBUG7,__FUNCTION__,
							"nondeterministic inputs " + logger::formatted_to_string(Minput.domain()));
					if (guard_poly_ptr) {
						LOGGER(DEBUG7,__FUNCTION__,
								"using guard " + logger::formatted_to_string(*guard_poly_ptr));
						//					positional_vdomain inv_vars(sinv_ptr->get_variable_ids());
						//					if (!Minput.domain().contains_variables(inv_vars)) {
						//						affine_map<scalar_type> M_proj = affine_map<scalar_type>::projection_map(M.domain(),domt);
						//						M = concatenate(M,M_proj);
						//					} else {
						con_poly_ptr U_poly = con_poly_ptr(guard_poly_ptr->clone());
						if (inst.sinv_orig_ptr) {
							U_poly->add_constraints(inst.sinv_orig_ptr->get_constraints());
						}
						/** Project away unused variables
						 * @note If this is too costly, one could project using a map */
						variable_id_set input_variables = Minput.domain().get_variable_ids();
						variable_id_set noninput_variables = U_poly->get_variable_ids();
						set_difference_assign(noninput_variables, input_variables);
						U_poly->existentially_quantify_variables(noninput_variables);
						mapped_U = U_poly;
						//					}
					} else {
						// make a universe set for this domain
						mapped_U = support_function_provider::const_ptr(
								new hyperbox<scalar_type>());
					}
					mapped_U = support_function_provider::const_ptr(
							new sf_unary<scalar_type>(mapped_U, Minput));
				}

				if (M.domain() != inst.source_state_dom) {
					LOGGER_OS(LOW,__FUNCTION__) << "Map = " << M << ", domain: " << M.domain() << " should be: " << inst.source_state_dom;
					throw basic_exception(
							"Domain of map does not correspond to source states.");
				}

				reset_trans_ptr new_trans(new reset_trans_type(M,mapped_U));

				LOGGER(DEBUG7,__FUNCTION__,
						"final map is " + logger::formatted_to_string(new_trans) + " with domain " + logger::formatted_to_string(new_trans->domain()));

				tcons.set_transform(new_trans);
			}
		}

		/** Project the guard to the state variables
		 *
		 * @note: This is done after computing the input set U,
		 * for which the non-projected guard is used. */
		guard_poly_ptr->existentially_quantify_variables(source_nonstate_variables);
//		LOGGER_OS(DEBUG5,__FUNCTION__) << "exist" << guard_poly;
		/** Simplify the guard constraint representation */
//		LOGGER_OS(DEBUG5,"discrete_post_stc::post") << "guard cons after rev map:" << guard_poly;
		guard_poly.remove_redundant_constraints();
//		LOGGER_OS(DEBUG5,__FUNCTION__) << "remove_red" << guard_poly;
//		LOGGER_OS(DEBUG5,"discrete_post_stc::post") << "guard constraints after redundancy check: " << guard_poly;
		//collapse guard constraints
		guard_poly.collapse_inequalities();
		LOGGER_OS(DEBUG7,__FUNCTION__)
				<< "guard constraints list after collapsing : " << guard_poly << " empty: " << std::boolalpha << guard_poly.is_empty();

		/** Set the new guard */
		tcons.set_guard(guard_poly_ptr);

		// Define the new jump constraints
		t->set_jump_constraints(tcons);
		return t;
	}
	;
};

/** Converts a hybrid automaton from DAE to ODE form
 * (when possible)
 *
 * Conversion is done on the fly.
 * Currently no caching.
 */
template<typename scalar_type>
hybrid_automaton::ptr dae_to_ode_adapter(hybrid_automaton::ptr impl) {
	location_adapter::ptr loc_ad(
			new dae_to_ode_location_adapter<scalar_type>());
	transition_adapter::ptr trans_ad(new dae_to_ode_transition_adapter<scalar_type>());

	hybrid_automaton_wrapper::ptr adapt(new hybrid_automaton_on_the_fly_adapter(loc_ad, trans_ad));
	adapt->set_impl(impl);
	return adapt;
}

}

#endif /* DAE_TO_ODE_ADAPTER_H_ */
