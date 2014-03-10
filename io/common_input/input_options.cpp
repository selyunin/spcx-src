/*
 * input_options.cpp
 *
 *  Created on: Sep 28, 2009
 *      Author: frehse
 */

#include "io/common_input/input_options.h"

#include <iostream>
#include <fstream>
#include <sstream>

#include "application/options.h"
#include "utility/logger.h"
#include "utility/stl_helper_functions.h"
#include "core/hybrid_automata/location.h"
#include "core/hybrid_automata/hybrid_automaton.h"
#include "core/hybrid_automata/hybrid_automaton_utility.h"
#include "core/hybrid_automata/hybrid_automaton_network.h"
#include "core/hybrid_automata/automaton_cache.h"
#include "core/hybrid_automata/adapt_automaton_visitor.h"
#include "core/hybrid_automata/automaton_name_formatter.h"
#include "core/symbolic_states/symbolic_state_collection.h"
#include "core/scenarios/parse_with_scenario.h"
#include "core/scenarios/scenario_chooser.h"
#include "math/vdom/context_variable_formatter.h"
#include "io/TXT_format/TXT_formatter.h"

#include "core/hybrid_automata/adaptors/dae_to_ode_adapter.h"

namespace options {

std::string get_system_name(options::options_processor::variables_map& vmap) {
	std::string sys_name = "system";
	options::options_processor::get_string_option(vmap,"system",sys_name);
	return sys_name;
}

hybrid_automata::hybrid_automaton::ptr get_system(
		options::options_processor::variables_map& vmap) {
	std::string sys_name = get_system_name(vmap);

	hybrid_automata::hybrid_automaton::ptr aut_net =
			hybrid_automata::hybrid_automaton_cache::get_automaton(sys_name);

	if (!aut_net) {
		std::string msg="Trying to analyze component \""+sys_name+"\", but could not find a component with that name.";
		if (sys_name=="system") {
			msg+=" Did you forget to specify the name of the component you want to analyze?";
		}
		msg+=" Use option --list-components to obtain information about known components.";
		throw basic_exception(msg);
	}

	return aut_net;
}

hybrid_automata::hybrid_automaton::ptr adapt_system(
		options::options_processor::variables_map& vmap) {
	hybrid_automata::hybrid_automaton::ptr aut_net = get_system(vmap);

	const hybrid_automata::reachability_scenario& scen =
			hybrid_automata::scenario_chooser::get_scenario();

	// get an adaptor for adapting the system
	hybrid_automata::adapt_automaton_visitor_ptr adaptor =
			scen.create_adapt_automaton_visitor();

	try {
		//parse_type_chooser::set_bool(adaptor->get_bool_type());
		//parse_type_chooser::set_number(adaptor->get_number_type());
		adaptor->init();

		aut_net->accept(*adaptor);
	} catch (std::exception& e) {
		throw basic_exception("Failed to adapt automaton " + get_system_name(
				vmap) + " to scenario " + scen.get_name() + ".", e);
	}

	if (!adaptor->get_success()) {
		throw basic_exception("Failed to adapt automaton " + get_system_name(
				vmap) + " to scenario " + scen.get_name() + ".");
	}

	/** Temporary adapter hack for dae to ode adapting */
	if (scen.get_name()!="phaver") {
	aut_net = hybrid_automata::dae_to_ode_adapter<double>(aut_net);
	hybrid_automata::hybrid_automaton_cache::add_automaton(aut_net);
	}

	return aut_net;
}

hybrid_automata::symbolic_state_collection_ptr define_initial_states(
		options::options_processor::variables_map& vmap,
		hybrid_automata::hybrid_automaton::ptr sys) {
	hybrid_automata::symbolic_state_collection_ptr ini =
			hybrid_automata::symbolic_state_collection_ptr();
	hybrid_automata::symbolic_state_collection_ptr model_ini =
			sys->get_initial_states();
	if (model_ini) {
		hybrid_automata::canonicalize_location_constraints(*model_ini);
	}
	std::string ini_states_string;
	if (options::options_processor::get_string_option(vmap,"initially",ini_states_string)) {
		const hybrid_automata::reachability_scenario& scen =
				hybrid_automata::scenario_chooser::get_scenario();

		std::string sys_name = get_system_name(vmap);

		parser::parse_policy ppol = parser::parse_policy::SX_policy();
		// accept unknown vars so that subcomponent variables can be specified too
		ppol.add_unknown_vars = true;
		try {
		ini = hybrid_automata::parse_symbolic_state_collection_with_scenario(
				scen, ini_states_string, sys_name,ppol);
		} catch ( std::exception& e ) {
			throw basic_exception("Failed to define initial states.",e);
		}

		if (ini) {
			hybrid_automata::canonicalize_location_constraints(*ini);
			std::string s;
			bool override = options::options_processor::get_string_option(vmap,
					"override-model-initial-states", s) && (s == "no");
			if (model_ini && override) {
				ini->intersection_assign(model_ini);
			}
		}
	} else
		ini = model_ini;
	if (!ini) {
		throw basic_exception(
				"No initial states specified. Please specify in options or in model file.");
	}

	// Check for variables
	variable_id_set vis = ini->get_variable_ids();
	throw_if_not_contains(variable_id_set(), get_primed_variables(vis),
			"Primed variable(s) not permitted in initial states, but found ",
			".");
	throw_if_not_contains(sys->get_variable_ids(), vis, "Unknown variable(s) ",
			" in initial states.");

	return ini;
}

std::string format_ident(const std::string& ident, const std::string& h) {
	std::string new_str = dot_context::in_context_name(ident,h);
	return new_str;
}

std::string format_var(const std::string& ident, const std::string& h) {
	return context_variable_formatter(h).format(ident);
}

std::string format_aut(const std::string& ident, const std::string& h) {
	std::string new_str = dot_context::in_context_name(ident, h);
	//	replace_first(new_str, h + ".", "");
	if (new_str.find(".") != std::string::npos) {
		std::set<std::string> known =
				hybrid_automata::hybrid_automaton_cache::get_automaton_names();
		new_str = dot_context::smallest_unique_name(new_str, h,
				known.begin(), known.end());
	}
	return new_str;
}

void add_base_components(std::list<std::string>& component_list,
		hybrid_automata::hybrid_automaton::ptr h) {
	// List the leaf components of a network
	if (hybrid_automata::hybrid_automaton_network::ptr aut_net=boost::dynamic_pointer_cast<hybrid_automata::hybrid_automaton_network>(h)) {
		hybrid_automata::automaton_id_set aut_ids = aut_net->get_automata();
		for (hybrid_automata::automaton_id_set::const_iterator it =
				aut_ids.begin(); it != aut_ids.end(); ++it) {
			hybrid_automata::hybrid_automaton::ptr a =
					hybrid_automata::hybrid_automaton_cache::get_automaton(*it);
			add_base_components(component_list, a);
		}
	} else {
		component_list.push_back(h->get_name());
	}
}

std::list<std::string> location_names(hybrid_automata::hybrid_automaton::ptr h) {
	std::list<std::string> loc_names;
	std::pair<hybrid_automata::hybrid_automaton::location_const_iterator,
			hybrid_automata::hybrid_automaton::location_const_iterator> I =
			h->get_locations();
	for (hybrid_automata::hybrid_automaton::location_const_iterator i = I.first; i
			!= I.second; ++i) {
		hybrid_automata::location_ptr l = h->get_location(*i);
		loc_names.push_back(l->get_name());
	}
	return loc_names;
}

void list_automaton_base_components(std::ostream& os,
		hybrid_automata::hybrid_automaton::ptr h, const std::string& context) {
	std::list<std::string> component_list;
	add_base_components(component_list, h);

	for (std::list<std::string>::const_iterator it = component_list.begin(); it
			!= component_list.end(); ++it) {
		if (it != component_list.begin())
			os << ",";
		os << format_aut(*it, context);
	}
}

void list_automaton_info(std::ostream& os, const std::string& aut_name, int show_subcomponents,
		const std::string& parent = "") {
	// Output variable+label info about the automaton to os.
	// Resolve identifiers w.r.t. to context aut_name
	// (e.g., output "x" instead of "sys.x" if context is "sys").
	// List the variables

	hybrid_automata::hybrid_automaton::ptr h =
			hybrid_automata::hybrid_automaton_cache::get_automaton(aut_name);

	// set the variable output to context lookup
	context_variable_formatter var_form(aut_name);
	os << var_form;

	variable_id_set vars = h->get_variable_ids();
	variable_id_set contrvars = h->get_variable_ids();
	variable_id_set inpvars = h->get_input_variables();
	variable_id_set constvars = h->get_const_variables();
	set_difference_assign(contrvars, inpvars);
	set_difference_assign(contrvars, constvars);

	if (contrvars.begin() != contrvars.end()) {
		std::cout << "Controlled = ";
		for (variable_id_set::const_iterator it = contrvars.begin(); it
				!= contrvars.end(); ++it) {
			if (it != contrvars.begin())
				os << ",";
			os << variable(*it);
		}
		os << std::endl;
	}

	if (inpvars.begin() != inpvars.end()) {
		std::cout << "Uncontrolled = ";
		for (variable_id_set::const_iterator it = inpvars.begin(); it
				!= inpvars.end(); ++it) {
			if (it != inpvars.begin())
				os << ",";
			os << variable(*it);
		}
		os << std::endl;
	}

	if (constvars.begin() != constvars.end()) {
		std::cout << "Constant-Dynamics = ";
		for (variable_id_set::const_iterator it = constvars.begin(); it
				!= constvars.end(); ++it) {
			if (it != constvars.begin())
				os << ",";
			os << variable(*it);
		}
		os << std::endl;
	}

	if (false) {
		// List labels
		hybrid_automata::label_id_set labels = h->get_labels();
		std::cout << "Synchronization Labels = ";
		for (hybrid_automata::label_id_set::const_iterator it = labels.begin(); it
				!= labels.end(); ++it) {
			if (it != labels.begin())
				os << ",";
			os << format_ident(hybrid_automata::named_label::get_name(*it),
					aut_name);
		}
		os << std::endl;
	}

	if (!parent.empty()) {
		os << "Parent = " << parent << std::endl;
	}

	// List Components if it's a network
	if (hybrid_automata::hybrid_automaton_network::ptr aut_net=boost::dynamic_pointer_cast<hybrid_automata::hybrid_automaton_network>(h)) {
		os << "Base-components = ";
		list_automaton_base_components(os, aut_net, aut_name);
		os << std::endl;
		// Define a default initial state (for editing by user)
		os << "Default-Ini = \"";
		if (h->get_initial_states()) {
			std::stringstream ss;
			io::TXT_formatter formatter(ss);
			formatter.set_context(aut_name);
			formatter.output(*h->get_initial_states());
//			hybrid_automata::context_automaton_name_formatter form("",aut_name);
//			ss << form;
//			ss << h->get_initial_states();
			// remove accolades
			std::string s = ss.str();
			trim_quotes(s, '{', '}');
			os << s;
		} else {
			variable_id_set vis = h->get_variable_ids();
			bool first_var = true;
			if (!vis.empty()) {
				for (variable_id_set::const_iterator it = vis.begin(); it
						!= vis.end(); ++it) {
					if (it != vis.begin())
						os << " & ";
					os << variable(*it) << "==0";
					first_var = false;
				}
			}
			std::list<std::string> component_list;
			add_base_components(component_list, h);
			if (!component_list.empty()) {
				bool first_loc = true;
				for (std::list<std::string>::const_iterator it =
						component_list.begin(); it != component_list.end(); ++it) {
					hybrid_automata::hybrid_automaton::ptr
							g =
									hybrid_automata::hybrid_automaton_cache::get_automaton(
											*it);
					std::list<std::string> loc_names = location_names(g);
					// only include a restriction on the location if there is more than one
					if (loc_names.size() > 1) {
						if (!first_var) {
							// if there's been a var constraint, and
							// if it's the very first location constraint,
							// then we need to connect it to the variable constraints
							os << " & ";
						}
						first_loc = false;
						os << "loc(";
						os << format_aut(*it, aut_name) << ")=="
								<< *loc_names.begin();
					}
				}
			}
		}
		os << "\"" << std::endl;

		os << std::endl;
		// Show subcomponents
		if (show_subcomponents!=0) {
			hybrid_automata::automaton_id_set aut_ids = aut_net->get_automata();
			for (hybrid_automata::automaton_id_set::const_iterator it =
					aut_ids.begin(); it != aut_ids.end(); ++it) {
				hybrid_automata::hybrid_automaton::ptr a =
						hybrid_automata::hybrid_automaton_cache::get_automaton(
								*it);
				std::cout << "[" << a->get_name() << "]" << std::endl;
				list_automaton_info(os, a->get_name(), show_subcomponents-1, aut_name);
			}
		}
	} else {
		// otherwise list its locations
		os << "Locations = ";
		std::list<std::string> loc_names = location_names(h);
		for (std::list<std::string>::const_iterator i = loc_names.begin(); i
				!= loc_names.end(); ++i) {
			if (i != loc_names.begin())
				os << ",";
			os << *i;
		}
		os << std::endl;
		// Define a default initial state (for editing by user)
		os << "Default-Ini = \"";
		if (h->get_initial_states()) {
			std::stringstream ss;
//			hybrid_automata::context_automaton_name_formatter form("",aut_name);
//			ss << form;
//			ss << h->get_initial_states();

			io::TXT_formatter formatter(ss);
			formatter.set_context(aut_name);
			formatter.output(*h->get_initial_states());

			// remove accolades
			std::string s = ss.str();
			trim_quotes(s,'{','}');
			os << s;
		} else {
			variable_id_set vis = h->get_variable_ids();
			if (!vis.empty()) {
				for (variable_id_set::const_iterator it = vis.begin(); it
						!= vis.end(); ++it) {
					if (it != vis.begin())
						os << " & ";
					os << variable(*it) << "==0";
				}
				os << " & ";
			}
			if (!loc_names.empty()) {
				os << "loc()==" << *loc_names.begin();
			}
		}
		os << "\"" << std::endl;
		os << std::endl;
	}

}

void add_input_options() {
	options::options_processor::config.add_options()("system",
			boost::program_options::value<std::string>(),
			"name of the component to be analyzed");
	options::options_processor::config.add_options()("model-file,m",
			boost::program_options::value<std::vector<std::string> >()->composing(),
			"model file(s)");
	options::options_processor::config.add_options()("initially,i",
			boost::program_options::value<std::string>(),
			"Set the initial states.");
	options::options_processor::config.add_options()(
			"override-model-initial-states",
			boost::program_options::value<std::string>()->implicit_value("no"),
			"(yes/no) By default, initial states passed in the options are combined (intersected) with the initial states passed in the model. If override is activated, the model initial state are not taken into account initial states provided in the model.");
	options::options_processor::config.add_options()(
			"list-components",
			boost::program_options::value<std::string>()->implicit_value(""),
			"Output all components and their variables. Optionally, the name of a component can be given, and only that component and its subcomponents are shown.");
}

bool check_input_options(options::options_processor::variables_map& vmap) {
	std::vector<std::string> model_files;
	if (options::options_processor::get_string_vector_option(vmap,"model-file",model_files)) {
		//std::cout << "Model files are: " << model_files << "\n";
		for (unsigned int j = 0; j < model_files.size(); ++j) {
			std::ifstream ifs(model_files[j].c_str(), std::ifstream::in);
			if (!ifs.good()) {
				throw basic_exception("Error opening model file "
						+ model_files[j] + " for reading.");
			}
			ifs.close();
		}
	} else {
		throw basic_exception("No model file specified.");
	}
	return true;
}

bool apply_input_options(options::options_processor::variables_map& vmap) {
	const hybrid_automata::reachability_scenario& scen =
			hybrid_automata::scenario_chooser::get_scenario();

	bool list_components = false;

	// Check if component list is requested. If so, output list and suspend execution.
	// Note: We don't list all components in the cache because this includes
	//       wrapper automata etc.
	std::string str;
	std::string list_components_context;
	if (options::options_processor::get_string_option(vmap,"list-components",list_components_context)) {
		list_components = true; // interrupt further execution
	}

	std::vector<std::string> model_files;
	if (options::options_processor::get_string_vector_option(vmap,"model-file",model_files)) {
		std::string instantiate_component = get_system_name(vmap);

		// If listing components, only parse requested
		if (list_components)
			instantiate_component = list_components_context;

		//std::cout << "Model files are: " << model_files << "\n";
		for (unsigned int j = 0; j < model_files.size(); ++j) {
			std::vector<hybrid_automata::hybrid_automaton_ptr> list_automata;
			list_automata = parse_models_in_file_with_scenario(scen,
					model_files[j], false, instantiate_component);
			std::string filename_wo_path = strip_path(model_files[j]);

			if (list_components) {
				// Provide a detailed list
				for (unsigned int i = 0; i < list_automata.size(); ++i) {
					std::string aut_name = list_automata[i]->get_name();
					hybrid_automata::hybrid_automaton::ptr
							h =
									hybrid_automata::hybrid_automaton_cache::get_automaton(
											aut_name);
					std::cout << "[" << aut_name << "]" << std::endl;

					int show_subcomponents = 0;
					if (!list_components_context.empty())
						show_subcomponents = 1;
					list_automaton_info(std::cout, aut_name, show_subcomponents);

					std::cout << std::endl;
				}

				// For debugging:
				//hybrid_automata::hybrid_automaton_cache::print(std::cout);
			} else {
				// Just produce a basic list of automata
				std::stringstream os;
				os << "Read file " << filename_wo_path
						<< ", defined automata ";
				for (unsigned int i = 0; i < list_automata.size(); ++i) {
					if (i > 0)
						os << ", ";
					os << list_automata[i]->get_name();
				}
				os << "." << std::endl << std::flush;
				LOGGER_OS(HIGH, "apply_input_options") << os.str();
			}
		}
	}

	return !list_components; // if components were listed, suspend execution
}

}
