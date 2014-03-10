/*
 * scenario_options.cpp
 *
 *  Created on: Nov 12, 2009
 *      Author: frehse
 */

#include "core/scenarios/scenario_options.h"

#include <iostream>
#include <boost/algorithm/string/trim.hpp>

#include "utility/stl_helper_functions.h"
//#include "../abstract_framework/hybrid_automata/automaton_cache.h"
#include "core/hybrid_automata/automaton_cache.h"
#include "core/hybrid_automata/hybrid_automaton.h"
#include "core/hybrid_automata/hybrid_automaton_utility.h"
#include "core/symbolic_states/symbolic_state_collection.h"
#include "core/scenarios/reachability_scenario.h"
#include "core/scenarios/scenario_chooser.h"
#include "core/scenarios/parse_with_scenario.h"
#include "core/analysis_algorithms/reachability_algorithm.h"
#include "core/continuous/support_function/template_directions/choose_directions.h"

#include "application/options.h"
//#include "application/output_options.h"

namespace options {

void add_scenario_options() {
	options::options_processor::config.add_options()(
			"scenario,c",
			boost::program_options::value<std::string>(),
			"Set the active scenario:\n- phaver : PHAVer (exact LHA)\n- stc : Space-Time Approximation with Clustering\n- supp : Support Functions (floating point, affine dynamics)\n- simu : Trajectory Simulation");
	options::options_processor::config.add_options()(
			"directions",
			boost::program_options::value<std::vector<std::string> >()->composing(),
			"Set the directions for template polyhedra:\n- box : box constraints\n- oct : octagonal constraints\n- uniXXX : XXX uniform constraints\n- a conjunction of linear constraints (in {})\n If the directions option is given more than once, the directions are added to one another.");
	options::options_processor::config.add_options()(
			"time-horizon",
			boost::program_options::value<std::string>(),
			"Set the time horizon for the analysis. A negative value denotes infinite time.");
	options::options_processor::config.add_options()("sampling-time",
			boost::program_options::value<std::string>(),
			"Set the sampling time for the analysis. A negative value denotes no sampling.");
	options::options_processor::config.add_options()(
			"intersection-error",
			boost::program_options::value<std::string>(),
			"(default -1) Set a bound on the absolute error bound for intersection, i.e., for image computation of discrete transitions (supp and stc scenario only). \nsupp scenario: If negative or unspecified, the outer approximation is used. If positive or zero, an optimization algorithm tries to bring the error below the given boudn (buggy).\nstc scenario: If negative or unspecified, the flowpipe error is used and the min. possible number of pieces of the outer approximation. If zero, the flowpipe error is used and all pieces of the outer approximation.");
	options::options_processor::config.add_options()(
			"intersection-error-rel",
			boost::program_options::value<std::string>(),
			"(default 0) Set a bound on the relative error bound for intersection, whenever the absolute error is exceeded (stc scenario only).");
	options::options_processor::config.add_options()(
			"intersection-method",
			boost::program_options::value<std::string>(),
			"(default lb_chull) Sets the type of intersection algorithm (supp scenario only):lb for naive lower bound search, lb_chull for lower bound search"
			" with branch and bound method, lb_simult for simultaneous solving of multiple lower bound search problems for optimization."
			"If unspecified, the chull algorithm is used for guard intersection. The split option with lb_chull defines the # sfm sections "
			"which intersects with the guard to chull before computing the final guard intersection set. If split option is not defined, then all the intersecting"
			" sfm section members are chull-ed before computing the intersection with the guard constraints."
			);
	options::options_processor::config.add_options()(
			"split",
			boost::program_options::value<std::string>(),
			"(default 0) Defines the sfm interval size for computing convex hull and intersecting with the guard constraints "
			"(supp scenario only with positive intersection error). split is set to 0 by default meaning that all the members of the sfm which "
			" intersects with the guard will be chull-ed before computing the intersection with the guard constraint in the discrete jump step."
			" split k would mean that at most k sfm members should be be chulled before intersecting with the guard constraint. "
			" Higher values of split should mean lesser precision but faster computation of discrete jump step."
			);
	options::options_processor::config.add_options()(
			"minbrak",
			boost::program_options::value<std::string>(),
			"(default gold_desc) Sets the method for bracketing a function minima (supp scenario only with positive intersection-error): "
			"Set [gold_desc] for naive magnification of the descend step "
			"with the factor of 1 + Golden ratio, [parab_desc] for using parabolic extrapolation based downhill decend.");
	options::options_processor::config.add_options()(
			"set-aggregation",
			boost::program_options::value<std::string>(),
			"(default chull) How the convex sets that constitute flowpipe are aggregated after clustering (supp and stc scenario only):\n- chull : convex hull\n- none : treat each convex set separately, which means multiplying the number of successor sets after each iteration.");
	options::options_processor::config.add_options()(
			"clustering",
			boost::program_options::value<std::string>(),
			"(default 30) supp scenario: Cluster the convex sets that constitute the flowpipe up to relative error XXX in percent before computing transition successors. 0 means no clustering. 100 means clustering into one single set. A value in between leads to clustering such that the relative distance (Hausdorff) to the original is below the given value (smaller values indicate higher accuracy). Clustering is applied before aggregation.");
	options::options_processor::config.add_options()(
			"flowpipe-tolerance",
			boost::program_options::value<std::string>(),
			"(default -1) If x >= 0, the flowpipe is computed such that the error estimate is below that value. The semantics of the error estimate depend on the scenario.");
	options::options_processor::config.add_options()(
			"flowpipe-tolerance-rel",
			boost::program_options::value<std::string>(),
			"(default 0) If x >= 0, the flowpipe is computed such that the relative error is guaranteed to lie below that value whenever the absolute error is exceeded. Only valid in the stc scenario. ");
	options::options_processor::config.add_options()(
			"error-model",
			boost::program_options::value<std::string>(),
			"(default forw) Determines how the overapproximation error is obtained.\n- forw: simple forward error\n- interp: interpolated error\n- interpfb: forward/backward interpolated error");
	options::options_processor::config.add_options()(
			"refine-max-iter",
			boost::program_options::value<std::string>(),
			"(default 1) Set the max number of local refinement loops for time elapse on affine dynamics (phaver scenario only). Zero means derivatives are given by quantifying over flow and invariant. Higher values iteratively quantify over the flow and the reachable states of the previous iteration.");

	options::options_processor::config.add_options()(
				"global-time-horizon",
				boost::program_options::value<std::string>(),
				"(default -1) Sets the global time horizon, in the sense of the integration. A negative number means infinity");
	options::options_processor::config.add_options()(
				"simu-intern-algo",
				boost::program_options::value<std::string>(),
				"(default 0) Sets the internal simulation algorithm (0='safe' & 1 = 'adapt).");
	options::options_processor::config.add_options()(
				"simu-init-sampling-points",
				boost::program_options::value<std::string>(),
				"(default 0) Sets the number of points taken for sampling the initial states, counted per convex set. If zero, takes the center of the bounding box.");
	options::options_processor::config.add_options()(
			"compact-sections",
			boost::program_options::value<std::string>(),
			"(default true) Let sfm sections be compacted by moving the initial states");
	options::options_processor::config.add_options()(
			"use-all-dirs",
			boost::program_options::value<std::string>(),
			"(default false) Propagate all constraints from pre to post in discrete post");

	options::options_processor::config.add_options()(
			"quantify-unused-jump-variables",
			boost::program_options::value<std::string>(),
			"(default exact) Specifies how unused jump variables are eliminated:\n- exact: quantifier elimination\n- outer: outer approximation");

	options::options_processor::config.add_options()(
			"flowpipe-filename",
			boost::program_options::value<std::string>(),
			"(default none) If specified, flowpipe plots are written to files of this name, appended with a number (stc scenario only).");
	options::options_processor::config.add_options()(
			"flowpipe-simplify-concave",
			boost::program_options::value<std::string>(),
			"(default true) If true, flowpipe evolutions are simplified to reduce their complexity. (stc scenario only)");
	options::options_processor::config.add_options()(
			"flowpipe-simplify-convex",
			boost::program_options::value<std::string>(),
			"(default true) If true, flowpipe evolutions are simplified to reduce their complexity. (stc scenario only)");
	options::options_processor::config.add_options()(
			"nonlinear.flowcut",
			boost::program_options::value<std::string>(),
			"(default false) If true, boundary of linearization domain will be intersected. (stc scenario only)");
	options::options_processor::config.add_options()(
			"nonlinear.linearization-error-factor",
			boost::program_options::value<std::string>(),
			"(default 10) The linearization error is obtained by multiplying the affine_error with this factor. Attention: Currently, only the absolute error is used. (stc scenario only)");
	options::options_processor::config.add_options()(
			"nonlinear.linearization-upscale-factor",
			boost::program_options::value<std::string>(),
			"(default 2) If the linearization fails, augment error bound by this factor. (stc scenario only)");
	options::options_processor::config.add_options()(
			"nonlinear.domain-distance-factor",
			boost::program_options::value<std::string>(),
			"(default 10) If the initial set is closer than this distance*affine_error.abs() to the domain, increase linearization error. (stc scenario only)");
	options::options_processor::config.add_options()(
			"nonlinear.refine-linearization-error",
			boost::program_options::value<std::string>(),
			"(default true) Refine linearization error with a-posteriori measurements. (stc scenario only)");
	options::options_processor::config.add_options()(
			"nonlinear.refine-linearization-with-reach",
			boost::program_options::value<std::string>(),
			"(default false) Refine linearization error with reachable states and repeat reach. (stc scenario only)");
	options::options_processor::config.add_options()(
			"nonlinear.default-dwell-time-factor",
			boost::program_options::value<std::string>(),
			"(default 0.1) If the dwell time estimation doesn't give a result, use time_horizon*factor. (stc scenario only)");
	options::options_processor::config.add_options()(
			"nonlinear.split-linearization",
			boost::program_options::value<std::string>(),
			"(default true) Split set if linearization error too large. If false, the error bound is increased instead. (stc scenario only)");
}

bool check_scenario_options(options::options_processor::variables_map& vmap) {
	return true;
}

bool apply_scenario_options_wo_system(options::options_processor::variables_map& vmap) {
	std::string s;
	if (options::options_processor::get_string_option(vmap,"scenario",s)) {
		hybrid_automata::scenario_chooser::set_scenario(s);
	}
	// @todo Check why this is here twice
	if (options::options_processor::get_string_option(vmap,"scenario",s)) {
			hybrid_automata::scenario_chooser::set_scenario(s);
	}
	return true;
}

bool apply_scenario_options_with_system(options::options_processor::variables_map& vmap) {
	// Hand the rest of the options over to the scenario itself
	hybrid_automata::reachability_scenario scen =
			hybrid_automata::scenario_chooser::get_scenario();
	scen.apply_options(vmap);

	hybrid_automata::scenario_chooser::set_scenario(scen);
	return true;
}

}
