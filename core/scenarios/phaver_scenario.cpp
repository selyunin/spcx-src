/*
 * phaver_scenario.cpp
 *
 *  Created on: Aug 31, 2009
 *      Author: frehse
 */

#include "core/scenarios/phaver_scenario.h"

#include "core/post_operators/pcd_post/constant_bound_time_elapse.h"
#include "core/post_operators/direct_discrete_post.h"
#include "core/pwl/PLWL_sstate_collection_stl_list.h"
#include "core/symbolic_states/symbolic_state_collection_stl_list.h"
#include "core/hybrid_automata/explicit_automaton.h"
#include "core/hybrid_automata/pairwise_hybrid_automaton_network.h"
#include "core/hybrid_automata/adaptors/ppl_adaptors.h"
#include "core/hybrid_automata/adaptors/constr_poly_adaptors.h"

namespace hybrid_automata {

void phaver_scenario_apply_options(reachability_scenario* sp,
		options::options_processor::variables_map& vmap) {

	std::string s;
	std::vector<std::string> svec;
	if (options::options_processor::get_string_option(vmap,"refine-max-iter",s)) {
		unsigned int r = from_string<double> (s);
		constant_bound_time_elapse_post::max_refine_count = r;
	}
}

reachability_scenario get_phaver_scenario() {
	reachability_scenario s;
	s.set_name("phaver");
	s.set_continuous_post_operator(continuous_post::ptr(new constant_bound_time_elapse_post));
	s.set_discrete_post_operator(discrete_post::ptr(new direct_discrete_post));
	passed_and_waiting_list::ptr PLWL_ptr(new PLWL_sstate_collection_stl_list);
	PLWL_ptr->set_option("MERGEP_ON");
	s.set_passed_and_waiting_list(PLWL_ptr);
	s.set_symbolic_state_collection(symbolic_state_collection::ptr(
			new symbolic_state_collection_stl_list));
	s.set_hybrid_automaton(hybrid_automaton::ptr(new explicit_automaton));
	s.set_hybrid_automaton_network(hybrid_automaton_network::ptr(
			new pairwise_hybrid_automaton_network));
	automaton_to_PPL_adaptor_ptr aut_ad=automaton_to_PPL_adaptor_ptr(new automaton_to_PPL_adaptor());
//	automaton_to_constr_poly_adaptor_ptr aut_ad = automaton_to_constr_poly_adaptor_ptr(
//			new automaton_to_constr_poly_adaptor(global_types::STD_BOOL,global_types::GMP_RATIONAL));
	aut_ad->init();
	s.set_adapt_automaton_visitor(aut_ad);

	s.set_option_handler(&phaver_scenario_apply_options);

	return s;
}

}
