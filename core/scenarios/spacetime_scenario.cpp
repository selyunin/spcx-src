/*
 * spacetime_scenario.cpp
 *
 *  Created on: Nov 3, 2012
 *      Author: notroot
 */

#include "spacetime_scenario.h"

#include "core/post_operators/spacetime_post/continuous_post_stc.h"
#include "core/post_operators/spacetime_post/discrete_post_stc.h"
#include "core/pwl/PLWL_sstate_collection_stl_list.h"
#include "core/symbolic_states/symbolic_state_collection_stl_list.h"
#include "core/hybrid_automata/explicit_automaton.h"
#include "core/hybrid_automata/pairwise_hybrid_automaton_network.h"
#include "core/hybrid_automata/adaptors/support_function_adaptors.h"
#include "io/common_input/input_options.h"
#include "io/common_input/parse_policy.h"

namespace hybrid_automata {

static const global_types::coefficient_type stc_number_type =
		global_types::STD_DOUBLE;

void spacetime_scenario_set_options(reachability_scenario* scp,
		const continuous_post_stc<double>::time_elapse_parameters& params, const continuous::spacetime_flowpipe<double>::cut_point_method& m) {
	if (params.affine_error.abs() < 0.0) {
		throw basic_exception("Flowpipe-tolerance needs to be >= 0 in scenario STC. Value is "+to_string(params.affine_error.abs()));
	}

	using namespace continuous;
	using namespace spacetime;
	using namespace math::numeric;

	//scp->set_intersection_error(m.approx_error.abs());
	typedef global_types::type_selector<stc_number_type>::type scalar_type;

	scp->set_continuous_post_operator(
			continuous_post::ptr(new continuous_post_stc<scalar_type>(params)));

	scp->set_discrete_post_operator(
			discrete_post::ptr(
					new discrete_post_stc<
							global_types::type_selector<stc_number_type>::type>(
							m)));

	// set the error for numeric computations in spacetime_plif
	// for now, use the same error as for comparing numbers
	spacetime::spacetime_plif::error_type num_err(
			approx_comparator<double>::get_rel_error_bound(),
			approx_comparator<double>::get_abs_error_bound());
	spacetime::spacetime_plif::set_numeric_error(num_err);
}

void spacetime_scenario_apply_options(reachability_scenario* sp,
		options::options_processor::variables_map& vmap) {
	/** @todo This should somehow not be necessary at this point. Global policy? */
	parser::parse_policy ppol = parser::parse_policy::SX_policy();
	// accept unknown vars so that subcomponent variables can be specified too
	ppol.add_unknown_vars = true;

	typedef global_types::type_selector<stc_number_type>::type scalar_type;
	continuous_post_stc<scalar_type>::time_elapse_parameters params;

	std::string s;
	std::vector<std::string> svec;
	if (options::options_processor::get_string_vector_option(vmap, "directions",
			svec)) {
		try {
			continuous::support_function::direction_chooser::reset();
			for (unsigned int i = 0; i < svec.size(); ++i) {
				s = svec[i];
				continuous::support_function::direction_chooser::add_directions(
						s, options::get_system_name(vmap), ppol);
			}
		} catch (std::exception& e) {
			std::string msg;
			if (s.size() >= 3 && s.substr(0, 3) == "uni") {
				msg =
						" For uniform directions, please specify the number of directions.";
			}
			throw basic_exception(
					"Could not parse the direction option \"" + s + "\"." + msg,
					e);
		}
	}
	if (options::options_processor::get_string_option(vmap, "time-horizon",
			s)) {
		sp->set_time_horizon(from_string<double>(s));
	}
	double inters_error_abs = -1.0; // default is same as flowpipe
	double inters_error_rel = 0.0; // default is no rel error
	double flowpipe_error_abs = 1.0; // default 1
	double flowpipe_error_rel = 0.0; // default is no rel error
	if (options::options_processor::get_string_option(vmap,
			"flowpipe-tolerance", s)) {
		flowpipe_error_abs = from_string<double>(s);
	}
	if (options::options_processor::get_string_option(vmap,
			"flowpipe-tolerance-rel", s)) {
		flowpipe_error_rel = from_string<double>(s);
	}
	if (options::options_processor::get_string_option(vmap,
			"intersection-error", s)) {
		inters_error_abs = from_string<double>(s);
	}
	if (options::options_processor::get_string_option(vmap,
			"intersection-error-rel", s)) {
		inters_error_rel = from_string<double>(s);
	}
	continuous::spacetime_flowpipe<double>::error_type flowpipe_error(flowpipe_error_rel,flowpipe_error_abs);

	typedef continuous::spacetime_flowpipe<double> flowpipe;
	flowpipe::cut_point_method m;
	m.type = flowpipe::cut_point_method::MIN_CONCAVE_PIECES;
	//m.type = spacetime_plif::cut_point_method::FIXED_COUNT_SAME_SIZE_PIECES;
	//m.piece_count = 1;
//						 m.type = spacetime_plif::cut_point_method::TIME_STEP;
//						 m.time_step = 0.5;
	m.approx_error = flowpipe_error; // default: same as flowpipe

	// by default, use the same error for intersection as for flowpipe
	if (inters_error_abs == 0.0) {
		m.type = flowpipe::cut_point_method::ALL_PIECES;
	} else if (inters_error_abs >= 0.0) {
		m.type = flowpipe::cut_point_method::MIN_CONCAVE_PIECES;
		m.approx_error = flowpipe::error_type(inters_error_rel,
				inters_error_abs);
	}

	if (options::options_processor::get_string_option(vmap,"set-aggregation",s)) {
		if (s == "chull") {
			m.type = flowpipe::cut_point_method::FIXED_COUNT_SAME_SIZE_PIECES;
			m.piece_count = 1;
			m.lower_is_error_reference = false;
		} else if (s.substr(0, 4) == "none") {
			// don't do anything
		} else
			throw basic_exception("Illegal set-aggregation \"" + s + "\"");
	}

	// MIN_CONCAVE_PIECES is used, the flowpipe error needs to be smaller than the intersection error!
	// here: let it be at most 50%
	if (m.type == flowpipe::cut_point_method::MIN_CONCAVE_PIECES) {
		flowpipe_error = continuous::spacetime_flowpipe<double>::error_type(
				0.125*flowpipe_error.rel(),std::min(flowpipe_error.abs(),0.125*m.approx_error.abs()));
		LOGGER(DEBUG, __FUNCTION__,
				"Corrected flowpipe error to "+to_string(flowpipe_error)+" in order to make clustering feasible for intersection error "+to_string(m.approx_error));
	}

	LOGGER(DEBUG, __FUNCTION__,
			"Set flowpipe error to "+to_string(flowpipe_error)+", intersection error to "+to_string(m.approx_error));

	if (options::options_processor::get_string_option(vmap,"flowpipe-filename",s)) {
		continuous::spacetime_flowpipe<double>::temp_filename = s;
	}

	if (options::options_processor::get_string_option(vmap,"flowpipe-simplify-concave",s)) {
		continuous::spacetime_flowpipe<double>::simplify_concave = from_string<bool>(s);
	}

	if (options::options_processor::get_string_option(vmap,"flowpipe-simplify-convex",s)) {
		continuous::spacetime_flowpipe<double>::simplify_convex = from_string<bool>(s);
	}

	if (options::options_processor::get_string_option(vmap,"nonlinear.flowcut",s)) {
		params.nonlinear_flowcut = from_string<bool>(s);
	}
	if (options::options_processor::get_string_option(vmap,"nonlinear.linearization-error-factor",s)) {
		params.linearization_error_factor = from_string<double>(s);
	}
	if (options::options_processor::get_string_option(vmap,"nonlinear.linearization-upscale-factor",s)) {
		params.linearization_upscale_factor = from_string<double>(s);
	}
	if (options::options_processor::get_string_option(vmap,"nonlinear.domain-distance-factor",s)) {
		params.domain_distance_factor = from_string<double>(s);
	}
	if (options::options_processor::get_string_option(vmap,"nonlinear.refine-linearization-error",s)) {
		params.refine_linearization_error = from_string<bool>(s);
	}
	if (options::options_processor::get_string_option(vmap,"nonlinear.refine-linearization-with-reach",s)) {
		params.refine_linearization_with_reach = from_string<bool>(s);
	}
	if (options::options_processor::get_string_option(vmap,"nonlinear.default-dwell-time-factor",s)) {
		params.default_dwell_time_factor = from_string<double>(s);
	}
	if (options::options_processor::get_string_option(vmap,"nonlinear.split-linearization",s)) {
		params.split_linearization = from_string<bool>(s);
	}

	/** Note: flowpipe error is the error made y flowpipe approximation, while
	 * m is the error made by convexification, clustering (before discrete map)
	 */
	params.affine_error = flowpipe_error;
	spacetime_scenario_set_options(sp, params, m);
}

reachability_scenario get_spacetime_scenario() {

	reachability_scenario s;
	s.set_name("stc");

	passed_and_waiting_list::ptr PLWL_ptr(new PLWL_sstate_collection_stl_list);
	PLWL_ptr->set_option("MERGEP_ON");
	s.set_passed_and_waiting_list(PLWL_ptr);
	s.set_symbolic_state_collection(symbolic_state_collection::ptr(
			new symbolic_state_collection_stl_list));
	s.set_hybrid_automaton(hybrid_automaton::ptr(new explicit_automaton));
	s.set_hybrid_automaton_network(hybrid_automaton_network::ptr(
			new pairwise_hybrid_automaton_network));

	// obtain default options using empty options map
	options::options_processor::variables_map empty_map;
	spacetime_scenario_apply_options(&s,empty_map);

	automaton_to_supp_f_adaptor_ptr aut_ad = automaton_to_supp_f_adaptor_ptr(
			new automaton_to_supp_f_adaptor(global_types::STD_BOOL,
					stc_number_type,false));
	aut_ad->init();
	s.set_adapt_automaton_visitor(aut_ad);
	s.set_option_handler(&spacetime_scenario_apply_options);
	return s;
}

}

