/*
 * parse_with_scenario.cpp
 *
 *  Created on: Sep 23, 2009
 *      Author: frehse
 */

#include "core/scenarios/parse_with_scenario.h"

#include <vector>
#include "core/scenarios/reachability_scenario.h"
#include "io/common_input/parse_type_chooser.h"
#include "io/common_input/state_parser.h"
#include "core/continuous/convert_predicate.h"
#include "core/hybrid_automata/adapt_automaton_visitor.h"
#include "core/hybrid_automata/automaton_cache.h"
#include "core/hybrid_automata/hybrid_automaton.h"
#include "core/hybrid_automata/hybrid_automaton_utility.h"
#include "core/symbolic_states/symbolic_state_collection.h"
#include "io/SX_format/SX_parser.h"
#include "io/CIF_format/CIF_parser.h"
#include "io/common_input/parse_type_chooser.h"
#include "io/common_input/symbol_table_cache.h"

#include "core/predicates/valuation_function_tree_nodes.h" // for locking context name lookup
#include "core/hybrid_automata/location_eq_node.h"

namespace hybrid_automata {

symbolic_state_collection::ptr parse_symbolic_state_collection_with_scenario(
		const reachability_scenario& scen, const std::string& s, const std::string& context,
		const parser::parse_policy& ppol) {
	adapt_automaton_visitor_ptr adaptor = scen.create_adapt_automaton_visitor();
	adaptor->init();
	parse_type_chooser::set_bool(adaptor->get_bool_type());
	parse_type_chooser::set_number(adaptor->get_number_type());

	// @todo This needs to be cleaned up. Locking should not be necessary thanks to parse policys. Include
	//       necessary info in parse policy.
	// @todo the lookup in the symbol table cache could be done inside predicate_parser::parse_state

	// lock the variable creation so that context dependent lookup is done
	bool oldlock=valuation_functions::variable_node_creator::is_locked();
	bool oldlock2=location_node_creator::is_locked();
	valuation_functions::variable_node_creator::set_locked(true);
	location_node_creator::set_locked(true);

	parser::symbol_table s_table = parser::symbol_table();
	if (ppol.use_symbol_table_cache) {
		s_table = parser::symbol_table_cache::get_symbol_table(context);
//std::cout << "parsing state collection using symbol stable:" << std::endl << s_table << std::endl;
	}
	tree::node::ptr p = predicate_parser::parse_state(s, context, s_table, ppol);

	valuation_functions::variable_node_creator::set_locked(oldlock);
	location_node_creator::set_locked(oldlock2);

	symbolic_state_collection::ptr sstates = scen.create_symbolic_state_collection();
	convert_and_add_predicate_to_symbolic_state_collection(p, *sstates);
	// Convert the predicates to continuous sets according to the scenario
	try {
		bool success = adapt(*adaptor, *sstates);
		if (!success) {
			throw basic_exception(
					"The following state predicate is not compatible with scenario "
							+ scen.get_name() + ": " + s);
		}
	} catch (std::exception& e) {
		throw basic_exception(
				"The following state predicate is not compatible with scenario "
						+ scen.get_name() + ".",e);
	}
	return sstates;
}

std::vector<hybrid_automata::hybrid_automaton::ptr> parse_models_in_file_with_scenario(const reachability_scenario& scen,
		const std::string& file_name, bool adapt_all, const std::string& sys_name) {

	// Get automata from XML file
	std::vector<hybrid_automata::hybrid_automaton::ptr> list_automata;
	adapt_automaton_visitor_ptr adaptor = scen.create_adapt_automaton_visitor();
	parse_type_chooser::set_bool(adaptor->get_bool_type());
	parse_type_chooser::set_number(adaptor->get_number_type());
	//automaton_parser::parse_automaton(file_name, list_automata);

	try {
		if (get_file_extension(file_name)=="cif")
			CIF_predicate_parser::parse_CIF(file_name, list_automata, sys_name);
		else
			parser::automaton_parser::parse_SX(file_name, list_automata, sys_name);
	} catch ( std::exception& e ) {
		throw basic_exception("Failed to parse model file "+strip_path(file_name) +".",e);
	}

	//hybrid_automaton_cache::print(cout);
	//	cout << std::flush;

	//std::cout << "parsed " << list_automatons.size() << " automata." << std::endl;

	if (adapt_all) {
		adaptor->init();
		for (unsigned int i = 0; i < list_automata.size(); ++i) {
			try {
				list_automata[i]->accept(*adaptor);
				if (!adaptor->get_success()) {
					throw basic_exception("Failed to adapt automaton "
							+ list_automata[i]->get_name() + " to scenario "
							+ scen.get_name());
				}
			} catch (std::exception& e) {
				throw basic_exception("Failed to adapt automaton "
						+ list_automata[i]->get_name() + " to scenario "
						+ scen.get_name(), e);
			}
			assert(adaptor->get_success());
		}
	}
	return list_automata;
}

hybrid_automaton::ptr parse_model_file_with_scenario(const reachability_scenario& scen,
		const std::string& file_name, const std::string& sys_name) {

	parse_models_in_file_with_scenario(scen, file_name, false, sys_name);

	// Get automata from XML file
	adapt_automaton_visitor_ptr adaptor = scen.create_adapt_automaton_visitor();
	parse_type_chooser::set_bool(adaptor->get_bool_type());
	parse_type_chooser::set_number(adaptor->get_number_type());

	adaptor->init();

	// Retrieve the network
	hybrid_automaton::ptr aut_net = hybrid_automaton_cache::get_automaton(sys_name);
	if (!aut_net)
		throw basic_exception("Could not find automaton "+sys_name+" in file " + file_name + ".");
	try {
		aut_net->accept(*adaptor);
		if (!adaptor->get_success()) {
			throw basic_exception("Failed to adapt system in file "
					+ file_name + " to current scenario.");
		}
	} catch (std::exception& e) {
		throw basic_exception("Failed to adapt system in file " + file_name
				+ " to current scenario.", e);
	}
	return aut_net;
}

}
