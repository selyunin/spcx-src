#include "SX_parser.h"

#include <iostream>
#include <map>
#include <string>
#include <boost/variant.hpp>

#include "extern/tinyxml/tinyxml.h"
#include "io/common_input/parse_exception.h"
#include "io/common_input/parser_basics.h"
#include "io/common_input/symbol_table.h"
#include "io/common_input/symbol_table_operators.h"
#include "io/common_input/network_creation.h"
#include "core/hybrid_automata/automaton_cache.h"
#include "core/hybrid_automata/hybrid_automaton.h"
#include "core/hybrid_automata/location.h"
#include "core/hybrid_automata/transition.h"
#include "core/continuous/continuous_set_transforms/continuous_set_transform.h"

#include "io/common_input/automaton_creation.h"
#include "io/common_input/symbol_table_cache.h"
#include "core/hybrid_automata/hybrid_automaton_network.h"
#include "core/predicates/dot_context_lookup.h"


using namespace dot_context;

namespace parser {
namespace automaton_parser {

/** Create a symbol or value from a TiXmlElement* param.
 *
 * The symbol table is used to look up the dimensions, which can
 * be symbols.
 * The value of the symbol is yet undetermined.
 * */
symbol get_unvalued_symbol(TiXmlElement* param, const symbol_table& s_table) {
	symbol res;
	std::string symbol_name(param->Attribute("name"));
	bool local = !(strcmp(param->Attribute("local"), "false")==0);
	bool is_label = (strcmp(param->Attribute("type"), "label")==0);

	symbol::symbol_type s_type = symbol::LABEL;
	symbol::dynamics_type dyn_type = symbol::ANY;
	symbol::data_type dat_type = symbol::REAL;
	bool controlled = true;
	unsigned int dim1 = 0;
	unsigned int dim2 = 0;

	if (!is_label) {
		/**
		 * Determinate dynamics_type, data_type and dimensions of parameter, if it is controlled.
		 */
		// If it's not a label, then it can only be a variable.
		// (constants are defined in map tags, not in param tags)
		s_type = symbol::VARIABLE;

		dyn_type = symbol::ANY;
		if (strcmp(param->Attribute("dynamics"), "const") == 0)
			dyn_type = symbol::CONSTANT;
		if (strcmp(param->Attribute("dynamics"), "explicit") == 0)
			dyn_type = symbol::EXPLICIT;

		if (strcmp(param->Attribute("type"), "int") == 0)
			dat_type = symbol::INT;

		dim1 = instantiate_dimension(param->Attribute("d1"),
				s_table);
		dim2 = instantiate_dimension(param->Attribute("d2"),
				s_table);

		if (param->Attribute("controlled") && strcmp(param->Attribute(
				"controlled"), "false") == 0)
			controlled = false;
	}

	res = symbol(symbol_name, s_type, dat_type, dyn_type, dim1, dim2, math::matrix<boost::any>(),
			local, controlled);

	return res;
};

/*!
 * Turns an xml element into an automaton object (with eventually parameters in a map)
 */
hybrid_automata::hybrid_automaton::ptr create_sspaceex_automaton(
		TiXmlElement* Temp, std::string name, bool bind,
		symbol_table s_table) {

	parse_policy ppol=parse_policy::SX_policy();

	hybrid_automata::hybrid_automaton::ptr aut;

	try {

		/**
		 * Create the symbol table of the component.
		 */
		std::map<std::string, std::string> id_to_name; // location ids
		aut = create_automaton(name);

//std::cout << "creating automaton " << name << std::endl;
//std::cout << "predefined symbol table:" << s_table << std::endl;

		/*!
		 * Read the template element, and for each child element:
		 */
		for (TiXmlElement* TempChild = Temp->FirstChildElement(); TempChild; TempChild
				= TempChild->NextSiblingElement()) {

			/*!
			 * If it's a param, instantiate the symbol table of the component and add variable to automaton
			 */
			if (strcmp(TempChild->Value(), "param") == 0) {
				symbol unv_symbol = get_unvalued_symbol(TempChild,s_table);
				symbol res_symbol = instantiate_symbol(unv_symbol,s_table,bind);
//std::cout << "adding symbol " << res_symbol << std::endl;
				add_symbol(aut,res_symbol);
			}
			else {
				// lock the symbol table
				s_table.set_locked();

//std::cout << "final symbol table:" << s_table << std::endl;
			/*!
			 * If it's a location:
			 */
			 if (strcmp(TempChild->Value(), "location") == 0) {
				std::string loc_id(TempChild->Attribute("id"));
				std::string loc_name(TempChild->Attribute("name"));
				id_to_name[loc_id] = loc_name;
				std::string invariant("");
				std::string flow("");

				for (TiXmlElement* CurrentElement =
						TempChild->FirstChildElement(); CurrentElement; CurrentElement
						= CurrentElement->NextSiblingElement()) {
					/*!
					 * Next, extract the invariant and flow, if they exist
					 */
					if (strcmp(CurrentElement->Value(), "invariant") == 0 && CurrentElement->FirstChild()) {
						invariant = CurrentElement->FirstChild()->Value();
					}
					if (strcmp(CurrentElement->Value(), "flow") == 0 && CurrentElement->FirstChild()) {
						flow = CurrentElement->FirstChild()->Value();
					}
				}
				/*!
				 * In the end, create the location and add to the automaton
				 */
				hybrid_automata::location::ptr l = create_location_from_string(
						loc_name, invariant, flow, aut->get_const_variables(), s_table, ppol);
				// Check whether variables are present
				throw_if_not_contains_context_free(
						aut->get_variable_ids(),
						l->get_time_constraints().get_invariant()->get_variable_ids(),
						"Unknown variable(s) ", " in invariant of location "
								+ l->get_name());
				throw_if_not_contains_context_free(
						aut->get_variable_ids(),
						get_unprimed_variables(
								l->get_time_constraints().get_dynamics()->get_variable_ids()),
						"Unknown variable(s) ", " in dynamics of location "
								+ l->get_name());

				aut->add_location(l);
			}

			/*!
			 * If it's a transition:
			 */
			else if (strcmp(TempChild->Value(), "transition") == 0) {
				/*!
				 * Copy the names of the source location and of the target location, by using the map id_to_name
				 */
				std::string source_loc_name(id_to_name[TempChild->Attribute(
						"source")]);
				std::string target_loc_name(id_to_name[TempChild->Attribute(
						"target")]);
				std::string label=hybrid_automata::named_label::silent_name();
				std::string guard(""), assignment("");

				/*!
				 * Extract the label, the guard and the assignment, if they exist
				 */
				for (TiXmlElement* CurrentElement =
						TempChild->FirstChildElement(); CurrentElement; CurrentElement
						= CurrentElement->NextSiblingElement()) {
					if (strcmp(CurrentElement->Value(), "label") == 0 && CurrentElement->FirstChild())
						label = CurrentElement->FirstChild()->Value();

					else if (strcmp(CurrentElement->Value(), "assignment") == 0 && CurrentElement->FirstChild())
						assignment = CurrentElement->FirstChild()->Value();

					else if (strcmp(CurrentElement->Value(), "guard") == 0 && CurrentElement->FirstChild())
						guard = CurrentElement->FirstChild()->Value();
				}

				/*!
				 * Add the transition to the automaton_template with automaton_template::add_transition_to_list
				 */
				hybrid_automata::location_id sloc = aut->get_location_id(
						source_loc_name);
				hybrid_automata::location_id tloc = aut->get_location_id(
						target_loc_name);

				std::string trans_message = "transition from location "
								+ source_loc_name + " to location "
								+ target_loc_name;

				if (label != hybrid_automata::named_label::silent_name()) {
						if (s_table.is_symbol(label)) {
							label = s_table.get_symbol(label).my_name;
						} else {
							throw parse_exception(
									"Unknown label " + label + ".");
						}
					} else {
						// add the silent label to the alphabet
						aut->add_label(
								hybrid_automata::named_label::silent_name());
					}

				hybrid_automata::transition::ptr t =
						create_transition_from_string(sloc, label, guard,
								assignment, tloc,
								aut->get_controlled_variables(), aut->get_input_variables(), aut->get_const_variables(), s_table, ppol);

				// check if variables exist
				throw_if_not_contains_context_free(
						aut->get_variable_ids(),
						t->get_jump_constraints().get_guard()->get_variable_ids(),
						"Unknown variable(s) ",
						" in guard of "+trans_message+".");
				variable_id_set used_vars;
				variable_id_set modif_vars;
				t->get_jump_constraints().get_transform()->get_used_and_modif_variables(
						used_vars, modif_vars);

				throw_if_not_contains_context_free(aut->get_variable_ids(),
						get_unprimed_variables(used_vars),
						"Unknown variable(s) ",
						" used in "+trans_message+".");
				throw_if_not_contains_context_free(aut->get_variable_ids(),
						get_unprimed_variables(modif_vars),
						"Unknown variable(s) ",
						" assigned in "+trans_message+".");

				aut->add_transition(t,false);
			}
			}
		}

		/** According to the policy, add symbol table to the cache. */
		if (ppol.use_symbol_table_cache) {
			symbol_table_cache::get_symbol_table(name)=s_table;
		}
	} catch (std::exception& e) {
		throw parse_exception("Could not parse base component " + name + ".", e);
	}

//	// Debug: print all automata in the cache
//	std::cout << std::endl << "After creating automaton " << name << std::endl;
//	hybrid_automata::hybrid_automaton_cache::print_all_automata(std::cout);

	return aut;
}

/*!
 * Turns an xml element into an automaton object (with eventually parameters in a map)
 */
hybrid_automata::hybrid_automaton_network::ptr create_sspaceex_network(
		TiXmlElement* Temp,
		std::string name,
		bool bind,
		std::map<std::string, std::pair<TiXmlElement*, bool> > list_automatons_node, symbol_table s_table) {

	parse_policy ppol=parse_policy::SX_policy();
	hybrid_automata::hybrid_automaton_network::ptr aut_net;

	try {

		/**
		 * Create the symbol table of the component.
		 */
		aut_net = create_automaton_network(name);
		s_table.set_context(name);

		/*!
		 * Read the template element, and for each child element:
		 */
		for (TiXmlElement* TempChild = Temp->FirstChildElement(); TempChild; TempChild
				= TempChild->NextSiblingElement()) {

			/*!
			 * If it's a param, instantiate the symbol table of the component and add variable to automaton
			 */
			if (strcmp(TempChild->Value(), "param") == 0) {
				symbol unv_symbol = get_unvalued_symbol(TempChild,s_table);
				symbol res_symbol = instantiate_symbol(unv_symbol,s_table,bind);
			}
			/*!
			 * If it's a bind:
			 */
			if (strcmp(TempChild->Value(), "bind") == 0) {
				std::string bind_name = TempChild->Attribute("as");
				std::string aut_name = name + "." + bind_name;

				/** New symbol table for the instantiation. */
				symbol_table aut_s_table(aut_name);

				for (TiXmlElement* CurrentElement =
						TempChild->FirstChildElement(); CurrentElement; CurrentElement
						= CurrentElement->NextSiblingElement()) {
					if (strcmp(CurrentElement->Value(), "map") != 0) {
						// ignore
						// throw parse_exception("Incorrect bind in component "+aut_name+", could not find map, got " + std::string(CurrentElement->Value()));
					} else {

						/**
						 * If the value of map is a number, add in key_map the pair "key, value"
						 * Else (it's a symbol), add in key_map the pair "key, s_table.get_symbol(value)"
						 */

						/** Check if key and value exist */
						if (!CurrentElement->Attribute("key"))
							throw parse_exception(
									"Incorrect bind in component " + aut_name
											+ ", could not find key.");
						std::string key = CurrentElement->Attribute("key");
						if (!(CurrentElement->FirstChild()
								&& CurrentElement->FirstChild()->Value())) {
							throw parse_exception(
									"Incorrect bind in component " + aut_name
											+ ", could not find value or mapping of parameter "
											+ key
											+ ".\n Note that every non-local parameter must be bound to either a constant value or to a parameter of the network component.");
						}

						std::string value =
								CurrentElement->FirstChild()->Value();
						if (is_number(value)) {
							// add the constant to the symbol table
							symbol const_symb(key, symbol::CONST_VALUE,
									symbol::UNINTERPR_STRING, symbol::CONSTANT,
									value);
							aut_s_table.add_symbol(const_symb);
						} else {
							aut_s_table.add_symbol(key,
									s_table.get_symbol(value));
						}
					}
				}

				/**
				 * If the component is a automaton template, create automaton with create_sspaceex_automaton and add it to network
				 * If the component is a network, create automaton with create_sspaceex_network and add it to network
				 */
				std::string component_name = TempChild->Attribute("component");
				TiXmlElement * component_node =
						list_automatons_node[component_name].first;
				bool is_network = list_automatons_node[component_name].second;
				if (component_node) {
					hybrid_automata::hybrid_automaton::ptr new_aut;
//std::cout << "bind with symbol table " << aut_s_table << std::endl;
					if (is_network)
						new_aut = create_sspaceex_automaton(component_node,
								aut_name, true, aut_s_table);
					else
						new_aut = create_sspaceex_network(component_node,
								aut_name, true, list_automatons_node, aut_s_table);
					aut_net = aut_net->compute_or_assign_composition(new_aut);
				} else
					throw parse_exception("Component " + component_name
							+ " is not declared.");
			}
		}

		/** According to the policy, add symbol table to the cache. */
		if (ppol.use_symbol_table_cache) {
			symbol_table_cache::get_symbol_table(name)=s_table;
		}
	} catch (std::exception& e) {
		throw parse_exception(
				"Could not parse network component " + name + ".", e);
	}

	return aut_net;
}

void parse_SX(std::string source, std::vector<
		hybrid_automata::hybrid_automaton_ptr>& list_automatons, std::string component_to_create) {

	/*!
	 * First, load the xml file in a  TiXmlDocument
	 */
	TiXmlBase::SetCondenseWhiteSpace(true);
	TiXmlDocument doc(source.c_str());
	bool load = doc.LoadFile();

	if (!load) {
		std::stringstream ss;
		ss << "Error loading XML file \"" << source << "\""<< std::endl;
		ss << "TinyXML error #" << doc.ErrorId() << " : "
				<< doc.ErrorDesc();
		throw parse_exception(ss.str());
		return;
	}

	TiXmlElement* root = doc.RootElement();
	/*!
	 * Next, verify that the root of the xml file is sspaceex
	 */
	if (strcmp(root->Value(), "sspaceex") != 0) {
		throw parse_exception("Error in XML file: XML root incorrect");
	}

	/**
	 * Create automaton list, with name of template, TinyXml node and a boolean.
	 * This boolean is true for template automaton and false for network.
	 */
	std::map<std::string, std::pair<TiXmlElement*, bool> > list_automatons_node;
	TiXmlElement* DocChild = root->FirstChildElement();

	/*!
	 * Read the xml file to extract elements.
	 * Automaton template elements are turns into automaton, then adds them to the ListAutomatons list.
	 * When element is not a automaton template, turn it into network with create_sspaceex_network.
	 */
	for (; DocChild; DocChild = DocChild->NextSiblingElement()) {

		if (strcmp(DocChild->Value(), "component") != 0
				|| !DocChild->Attribute("id")) {
			throw parse_exception("Incorrect syntax in file");
		}

		std::string component_name = std::string(DocChild->Attribute("id"));

		/**
		 * Check if the component is template or network.
		 * To do that, search a location. If its exist, the component is a template, else is a network.
		 */
		bool is_template = false;
		for (TiXmlElement* search_loc = DocChild->FirstChildElement(); search_loc; search_loc
				= search_loc->NextSiblingElement())
			if (strcmp(search_loc->Value(), "location") == 0)
				is_template = true;

		/** Only create automaton if necessary. */
		bool instantiate_component = (component_to_create.empty()
				|| component_to_create == component_name);

		symbol_table s_table(component_name);

		if (is_template) {
			/**
			 * Add node to list_automatons_node.
			 * If it is a template, create an automaton with create_sspaceex_automaton.
			 * Add this automaton to list_automatons.
			 */
			list_automatons_node[component_name] = std::make_pair<
					TiXmlElement*, bool>(DocChild, true);

			if (instantiate_component)
				list_automatons.push_back(create_sspaceex_automaton(DocChild,
						component_name, false,s_table));
		}

		else { //is a network, with bind
			/**
			 * Add node to list_automatons_node.
			 * If it is a network, create an automaton with create_sspaceex_automaton.
			 * Add this automaton to list_automatons.
			 */
			list_automatons_node[component_name] = std::make_pair<
					TiXmlElement*, bool>(DocChild, false);

			if (instantiate_component)
				list_automatons.push_back(create_sspaceex_network(DocChild,
						component_name, false, list_automatons_node, s_table));
		}
	}
}

}
}
