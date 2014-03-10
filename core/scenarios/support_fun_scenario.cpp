/*
 * support_fun_scenario.cpp
 *
 *  Created on: Aug 31, 2009
 *      Author: frehse
 */

#include "core/scenarios/support_fun_scenario.h"

#include "core/post_operators/sfm_post/continuous_post_sfm.h"
#include "core/post_operators/sfm_post/continuous_post_inters.h"
#include "core/post_operators/sfm_post/continuous_post_sfm_inters.h"
#include "core/post_operators/sfm_post/discrete_post_sfm_inters.h"
#include "core/post_operators/sfm_post/discrete_post_branch_bound.h"
#include "core/post_operators/sfm_post/discrete_post_simult.h"
#include "core/post_operators/sfm_post/discrete_post_intv_split.h"
#include "core/post_operators/sfm_post/discrete_post_sf_poly.h"
#include "core/post_operators/sfm_post/omega_model_factory.h"
#include "core/pwl/PLWL_sstate_collection_stl_list.h"
#include "core/symbolic_states/symbolic_state_collection_stl_list.h"
#include "core/hybrid_automata/explicit_automaton.h"
#include "core/hybrid_automata/pairwise_hybrid_automaton_network.h"
#include "core/hybrid_automata/adaptors/support_function_adaptors.h"
#include "core/scenarios/phaver_scenario.h"
#include "io/common_input/input_options.h"
#include "io/common_input/parse_policy.h"

namespace hybrid_automata {

static const global_types::coefficient_type supp_f_number_type =
		global_types::STD_DOUBLE;

void support_fun_scenario_set_options(reachability_scenario* scp, double e =
		-1.0, double aggr_choice = 0, double clustering = 0,
		double flowpipe_tolerance = -1) {
	scp->set_intersection_error(e);

	if (e >= 0.0) {

		scp->set_continuous_post_operator(continuous_post::ptr(
				new continuous_post_sfm<global_types::type_selector<
						supp_f_number_type>::type> (flowpipe_tolerance)));

		std::string inters_type = scp->get_intersection_type();
		std::string minbrak_type = scp->get_minbrak_type();
		size_t split_size = scp->get_split_size();

		if(inters_type.compare("lb_chull_old")==0){
			if(split_size == 0){
				scp->set_discrete_post_operator(discrete_post::ptr(
						new discrete_post_branch_bound<global_types::type_selector<
							supp_f_number_type>::type> (aggr_choice, clustering, minbrak_type, e)));
			}
			else{
				scp->set_discrete_post_operator(discrete_post::ptr(
						new discrete_post_intv_split<global_types::type_selector<
							supp_f_number_type>::type> (aggr_choice, clustering, minbrak_type, split_size, e)));
			}
		}
		else if(inters_type.compare("lb_chull")==0){
			if(split_size == 0){
				// GF 2011-12-04 : New version shares code with simult function
				scp->set_discrete_post_operator(discrete_post::ptr(
						new discrete_post_simult<global_types::type_selector<
							supp_f_number_type>::type>(aggr_choice, clustering, minbrak_type, e, true)));
			}
			else{
				scp->set_discrete_post_operator(discrete_post::ptr(
						new discrete_post_intv_split<global_types::type_selector<
							supp_f_number_type>::type> (aggr_choice, clustering, minbrak_type, split_size, e)));
			}
		}
		else if(inters_type.compare("lb")==0){
			scp->set_discrete_post_operator(discrete_post::ptr(
					new discrete_post_sfm_inters<global_types::type_selector<
							supp_f_number_type>::type> (aggr_choice,clustering, minbrak_type, e)));
		}
		else if(inters_type.compare("lb_simult")==0){
			scp->set_discrete_post_operator(discrete_post::ptr(

					new discrete_post_simult<global_types::type_selector<
						supp_f_number_type>::type>(aggr_choice, clustering, minbrak_type, e, false)));
		}
		else{
			throw std::runtime_error("Intersection Type:"+to_string(inters_type)+" not known");
		}
	} else {
		scp->set_continuous_post_operator(continuous_post::ptr(
				new continuous_post_sfm<global_types::type_selector<
						supp_f_number_type>::type> (flowpipe_tolerance)));
		// old operator: discrete_post_sfm
		scp->set_discrete_post_operator(discrete_post::ptr(
				new discrete_post_sf_poly<global_types::type_selector<
						supp_f_number_type>::type> (aggr_choice, clustering)));
	}
}

void support_fun_scenario_apply_options(reachability_scenario* sp,
		options::options_processor::variables_map& vmap) {
	/** @todo This should somehow not be necessary at this point. Global policy? */
	parser::parse_policy ppol = parser::parse_policy::SX_policy();
	// accept unknown vars so that subcomponent variables can be specified too
	ppol.add_unknown_vars = true;

	std::string s;
	std::vector<std::string> svec;
	if (options::options_processor::get_string_vector_option(vmap,"directions",svec)) {
		try {
			continuous::support_function::direction_chooser::reset();
			for (unsigned int i = 0; i < svec.size(); ++i) {
				s = svec[i];
				continuous::support_function::direction_chooser::add_directions(
						s,options::get_system_name(vmap),ppol);
			}
		} catch (std::exception& e) {
			std::string msg;
			if (s.size() >= 3 && s.substr(0, 3) == "uni") {
				msg
						= " For uniform directions, please specify the number of directions.";
			}
			throw basic_exception("Could not parse the direction option \"" + s
					+ "\"." + msg, e);
		}
	}
	if (options::options_processor::get_string_option(vmap,"time-horizon",s)) {
		sp->set_time_horizon(from_string<double> (s));
	}
	if (options::options_processor::get_string_option(vmap,"sampling-time",s)) {
		sp->set_sampling_time(from_string<double> (s));
	}
	double inters_error = -1.0; // default is no precise intersection
	double aggr_choice = 0; // default is convex hull
	double clustering = 0.3; // default is clustering up to 30% error
	double flowpipe_tolerance = -1; // default is fixed step
	if (options::options_processor::get_string_option(vmap,"set-aggregation",s)) {
		if (s == "chull")
			aggr_choice = 0;
		else if (s.substr(0, 5) == "thull") {
			aggr_choice = 1;
			if (s.size() > 5)
				aggr_choice = from_string<double> (s.substr(5));
		} else if (s.substr(0, 4) == "none") {
			aggr_choice = -1;
		} else
			throw basic_exception("Illegal set-aggregation \"" + s + "\"");
	}
	if (options::options_processor::get_string_option(vmap,"intersection-error",s)) {
		inters_error = from_string<double> (s);
		sp->set_intersection_type("lb_chull"); // lb_search branch and bound set as default.
		sp->set_minbrak_type("gold_desc"); // Golden descend method set as default for min bracketing
		sp->set_split_size(0); // split size set to 0 by default, meaning complete chull.
	}
	if (options::options_processor::get_string_option(vmap,"intersection-method",s)) {
		sp->set_intersection_type(s);
	}
	if (options::options_processor::get_string_option(vmap,"minbrak",s)) {
		sp->set_minbrak_type(s);
	}
	// default value for split size
	sp->set_split_size(0);
	if (options::options_processor::get_string_option(vmap,"split",s)) {
		sp->set_split_size(from_string<size_t>(s));
	}

	if (options::options_processor::get_string_option(vmap,"clustering",s)) {
		clustering = from_string<double> (s) / 100;
		if (clustering < 0.0 || clustering > 1.0001) {
			throw basic_exception(
					"Clustering percentage needs to be between 0 and 100.");
		}
	}
	if (options::options_processor::get_string_option(vmap,"flowpipe-tolerance",s)) {
		flowpipe_tolerance = from_string<double> (s);
	}
	if (options::options_processor::get_string_option(vmap,"error-model",s)) {
		continuous::support_function::omega_model_factory::set_model(s);
	}

	if (options::options_processor::get_string_option(vmap,"compact-sections",s)) {
		continuous::support_function::sfm_section<double>::compact_sfms=from_string<bool>(s);
	}

	if (options::options_processor::get_string_option(vmap,"use-all-dirs",s)) {
		hybrid_automata::discrete_post_simult<double>::use_all_dirs=from_string<bool>(s);
	}

	if (options::options_processor::get_string_option(vmap,
			"quantify-unused-jump-variables", s)) {
		if (s == "exact") {
			hybrid_automata::discrete_post_sf_poly<double>::quantify_unused_jump_variables
					= true;
		} else if (s == "outer") {
			hybrid_automata::discrete_post_sf_poly<double>::quantify_unused_jump_variables
					= false;
		} else {
			throw basic_exception(
					"Unknown value for option quantify-unused-jump-variables: \""
							+ s + "\".");

		}
	}


	support_fun_scenario_set_options(sp, inters_error, aggr_choice, clustering,
			flowpipe_tolerance);

}

reachability_scenario get_support_fun_scenario() {

	reachability_scenario s = get_phaver_scenario();
	s.set_name("supp");

	//	s.set_continuous_post_operator(continuous_post::ptr(
	//			new continuous_post_sfm<global_types::type_selector<
	//					supp_f_number_type>::type> ));
	//	s.set_discrete_post_operator(discrete_post::ptr(
	//			new discrete_post_sfm_inters<global_types::type_selector<
	//					supp_f_number_type>::type> ));

	// let continuous and discrete post be chosen according to default options
	support_fun_scenario_set_options(&s);

	automaton_to_supp_f_adaptor_ptr aut_ad = automaton_to_supp_f_adaptor_ptr(
			new automaton_to_supp_f_adaptor(global_types::STD_BOOL,
					supp_f_number_type,false));
	aut_ad->init();
	s.set_adapt_automaton_visitor(aut_ad);
	s.set_option_handler(&support_fun_scenario_apply_options);
	return s;
}

}
