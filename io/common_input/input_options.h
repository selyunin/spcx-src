/*
 * input_options.h
 *
 *  Created on: Sep 28, 2009
 *      Author: frehse
 */

#ifndef input_options_H_
#define input_options_H_

#include "application/options.h"

/** Forward declarations */
namespace hybrid_automata {
class hybrid_automaton;
typedef boost::shared_ptr<hybrid_automaton> hybrid_automaton_ptr;
class symbolic_state_collection;
typedef boost::shared_ptr<symbolic_state_collection> symbolic_state_collection_ptr;
}

namespace options {

std::string get_system_name(options::options_processor::variables_map& vmap);

hybrid_automata::hybrid_automaton_ptr get_system(options::options_processor::variables_map& vmap);

hybrid_automata::hybrid_automaton_ptr adapt_system(
		options::options_processor::variables_map& vmap);


hybrid_automata::symbolic_state_collection_ptr define_initial_states(
		options::options_processor::variables_map& vmap,
		hybrid_automata::hybrid_automaton_ptr sys);

void add_input_options();
bool check_input_options(options::options_processor::variables_map& vmap);
bool apply_input_options(options::options_processor::variables_map& vmap);


/**
 * Remove the excess names to shorten identifier
 *
 * @param ident Name of the Identifier
 * @param h Name of the context (Component)
 * @return Shortened Identifier
 */
std::string format_ident(const std::string& ident, const std::string& h);

/**
 * Remove the excess names to shorten identifier
 *
 * @param ident Name of the automaton
 * @param h Name of the context (Component)
 * @return Shortened Identifier
 */
std::string format_aut(const std::string& ident, const std::string& h);

/**
 * Remove the excess names to shorten identifier
 *
 * @param ident Name of the Variable
 * @param h Name of the context (Component)
 * @return Shortened Identifier
 */
std::string format_var(const std::string& ident, const std::string& h);


}

#endif /* input_options_H_ */
