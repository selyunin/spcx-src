/*
 * specification_parser.cpp
 *
 *  Created on: Sep 27, 2009
 *      Author: gvincent
 */


#include <fstream>
#include "state_parser.h"
#include "specification_parser.h"
#include "specification_grammar.h"
#include "utility/tree_node.h"
#include "core/symbolic_states/symbolic_state_collection_stl_list.h"
#include "core/predicates/node_print_visitor.h"

using namespace parser;

namespace predicate_parser {

hybrid_automata::symbolic_state_collection_stl_list create_symbolic_state_collection(
		const tree::node::ptr& p) {
	//std::cout<<"\n in tree: "<<p<<"\n";
	hybrid_automata::symbolic_state_collection_stl_list stc;
	return stc;
}

void execute_reach(std::string const& s, hybrid_automata::symbolic_state_collection_stl_list const& collec){
	//std::cout<<"\n name: "<<s<<"\n";
}

void execute_commands(std::vector<command> const& command_list, symbol_table symbol_table = symbol_table()) {
	for (std::vector<command>::size_type i = 0; i < command_list.size(); ++i){

		if(command_list[i].my_name == "reach"){
			if(command_list[i].my_arguments.size() != 2)
				throw std::runtime_error("Bad number of parameters for command " + command_list[i].my_name);

			std::string aut(boost::get<std::string>(command_list[i].my_arguments[0].my_arg));
			state_ast ast(boost::get<state_ast>(command_list[i].my_arguments[1].my_arg));
			/*if(aut == NULL || ast == NULL)
				throw std::runtime_error("Bad type of parameters (type must be string and state predicate)");*/

			execute_reach(aut, create_symbolic_state_collection(make_state_tree(ast, symbol_table)));
		}
		else //only command reach are implemented
			throw std::runtime_error(command_list[i].my_name +" : Unknown command");

	}
}

void parse_specification(std::string const& file, symbol_table symbol_table) {
	std::ifstream command_file(file.c_str());
	std::string commands = "";

	if (command_file) {
		std::string line;
		while (std::getline(command_file, line)) {
			commands = commands + line;
		}
		command_file.close();
	} else
		throw std::runtime_error("Error loading");

	if (commands != "") {
		specification_grammar g;
				std::vector<command> command_list;
				std::string::const_iterator iter = commands.begin();
				std::string::const_iterator end = commands.end();
				bool r = boost::spirit::qi::phrase_parse(iter, end, g, boost::spirit::ascii::space, command_list);

				if (r && iter == end) {
					execute_commands(command_list, symbol_table);
				} else {
					/**
					* error handle : precise as to when the parser fails
					* (the command fails)
					*/
					throw std::runtime_error("Could not parse commands '" + commands + "' \nfrom '" + std::string(iter, end) + "'.");
				}
	}

}

}
