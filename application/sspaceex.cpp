//============================================================================
// Name        : sspaceex.cpp
// Author      : Goran Frehse
// Version     :
// Copyright   : (C) 2008
// Description : Hello World in C++, Ansi-style
//============================================================================

//unsigned int VERBOSE_LEVEL = 600011;
//unsigned int line_number = 0;

//#define NDEBUG


#include <iostream>
#include <fstream>

#include <glpk.h>

#include "application/options.h"
#include "core/hybrid_automata/hybrid_automaton.h"
#include "core/hybrid_automata/hybrid_automaton_utility.h"
#include "core/hybrid_automata/automaton_cache.h"
#include "core/symbolic_states/symbolic_state_collection.h"
#include "core/scenarios/scenario_chooser.h"
#include "core/scenarios/reachability_scenario.h"
#include "core/scenarios/parse_with_scenario.h"

#include "core/analysis_algorithms/monitoring/monitoring_options.h"
#include "io/common_input/input_options.h"
#include "core/scenarios/scenario_options.h"
#include "core/analysis_algorithms/reachability_options.h"
#include "io/common_output/output_options.h"
#include "base_options.h"
#include "core/hybrid_automata/automaton_name_formatter.h"
#include "math/vdom/context_variable_formatter.h"

#include "utility/logger.h"
#include "utility/logger_stopwatch.h"
#include "utility/basic_warning.h"

#include "application/version_info.h"

using namespace std;

int main(int ac, char* av[]) {
#ifdef NDEBUG
	try {
#endif
	hybrid_automata::context_automaton_name_formatter
			global_automaton_formatter("", "system");
	// set variable output to be context dependent
	context_variable_formatter global_variable_formatter("system");

	logger::set_stream(std::cout);

	// Register options across components
	options::add_base_options();
	options::add_scenario_options();
	options::add_input_options();
	options::add_output_options();
	options::add_reachability_options();
	options::add_monitoring_options();

	// Define options (visible, hidden, ...)
	options::options_processor::define_options();

	// Check if additional-options is set
    for (int i = 1; i + 1 < ac; i++) {
		if (std::string(av[i]) == "--additional-options") {
			options::options_processor::read_command_line_from_string(
					std::string(av[i + 1]));
		}
	}

	// Read options from command line
	int command_options_exit = options::options_processor::read_command_line(
			ac, av);

//	if (options::options_processor::get_string_option(vmap,
//			"additional-options", str)) {
//		options::options_processor::read_command_line_from_string(str);
//	}


	if (command_options_exit != 0) {
		cerr
				<< "Error in command line options.\nTry `--help' for more information."
				<< std::endl;
	} else {
		options::options_processor::variables_map& vmap =
				options::options_processor::get_vmap();

		options::set_verbosity(vmap);

		LOGGER(DEBUG, "main", "Processing options");

		std::string str;
		if (options::options_processor::get_string_option(vmap, "config", str)) {
			options::options_processor::read_config_file(str);
		}

		if (vmap.size() ==0) {
			std::cout << "Missing arguments." << std::endl << "Try `--help' for more information." << std::endl;
		} else if (vmap.count("help")) {
			std::cout << options::options_processor::get_visible() << "\n";
		} else if (vmap.count("version")) {
			std::cout << "SpaceEx State Space Explorer, v"
					<< sspaceex::get_version() << ", compiled "
					<< sspaceex::get_compilation_time() << ", "
					<< global_types::coefficient_type_summary() << "\n";
		} else {
			bool exec = true;
			// Check options for consistency etc.
			if (vmap.count("show-option-values")) {
				LOGGER_OS(LOW,"main") << "Option values:" << std::endl;
				IFLOGGER(LOW) {
					options::options_processor::print_options(vmap,logger(logger_level::LOW,"main").buffer());
				}
			}

			LOGGER(DEBUG, "main", "Checking options");
			exec &= options::check_base_options(vmap);
			exec &= options::check_scenario_options(vmap);
			exec &= options::check_input_options(vmap);
			exec &= options::check_reachability_options(vmap);
			exec &= options::check_monitoring_options(vmap);

			// Apply options
			if (exec) {
				LOGGER(DEBUG, "main", "Applying options");
				// apply the scenario options that don't need to know about the system
				exec &= options::apply_scenario_options_wo_system(vmap);
				exec &= options::apply_base_options(vmap);
				// now the system gets defined
				exec &= options::apply_input_options(vmap);
				if (exec) {
					/** Apply logger formatter */
					global_automaton_formatter.set_context(options::get_system_name(vmap));
					global_variable_formatter.set_context(options::get_system_name(vmap));
					*logger::get_stream() << global_automaton_formatter;
					*logger::get_stream() << global_variable_formatter;

					// checked after system is defined
					exec &= options::apply_scenario_options_with_system(vmap);
					exec &= options::check_output_options(vmap);
					exec &= options::apply_output_options(vmap);
					exec &= options::apply_reachability_options(vmap);
					exec &= options::apply_monitoring_options(vmap);
				}
			}

			// Execute application

			if (exec) {
#ifndef NDEBUG
		LOGGER(LOW, "main", "Running with assertion checking on.");
#endif

				// choose whether to output
				std::vector<io::output_formatter*> of =
						options::create_output_formatter(vmap);
				for (unsigned int i = 0; i < of.size(); ++i) {
					of[i]->prologue();
				}

				// choose which system to analyze
				LOGGER(DEBUG, "main", "Adapting system to scenario");
				hybrid_automata::hybrid_automaton::ptr aut_net =
						options::adapt_system(vmap);

				// choose the initial states
				LOGGER(DEBUG, "main", "Defining initial states");
				hybrid_automata::symbolic_state_collection_ptr ini_states =
						options::define_initial_states(vmap, aut_net);
				// intersect with the invariant
				ini_states = hybrid_automata::intersect_with_invariants(
						*ini_states, *aut_net);

				// choose which algorithm to apply
				if (options::use_monitoring(vmap)) {
					options::monitor(aut_net, ini_states, vmap, of);
				} else {
					options::execute_reachability(aut_net, ini_states, vmap, of);
				}

				// output system automaton etc.
				options::apply_post_analysis_output_options(vmap);

				// Finish output and terminate application
				for (unsigned int i = 0; i < of.size(); ++i) {
					of[i]->epilogue();
				}
				options::destroy_output_formatter(of);
			}
			options::options_processor::clean_up();
		}
	}

	basic_warning::output_statistics();
	LOGGER(DEBUG, "main", "Done.");

	logger_stopwatch::report_all(logger_level::DEBUG);

	logger::terminate();

#ifdef NDEBUG
} catch (exception& e) {
	std::cerr << std::endl << "ERROR: "<< e.what() << "\n";
	return 1;
}
#endif
	return 0;
}
