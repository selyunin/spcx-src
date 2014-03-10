#ifndef SX_PARSER_H_
#define SX_PARSER_H_

#include <string>
#include <vector>
#include <map>
#include <boost/shared_ptr.hpp>

/** Forward declarations */
namespace hybrid_automata {
class hybrid_automaton;
typedef boost::shared_ptr<hybrid_automaton> hybrid_automaton_ptr;
class hybrid_automaton_network;
typedef boost::shared_ptr<hybrid_automaton_network> hybrid_automaton_network_ptr;
}

namespace parser {
namespace automaton_parser {

/*!
 * \fn parse_SX(string source, vector<hybrid_automaton_ptr>& list_automatons)
 * \brief The main Parser function, which parse a xml file and return a vector of automaton
 *
 * \param source the address of the xml file
 * \param list_automatons a list of where to put the automata extracted from the xml file
 * \param component_to_create is the name of the component to be instantiated as an automaton.
 * If component_to_create is omitted or empty, all components are instantiated.
 */
void parse_SX(std::string source, std::vector<
		hybrid_automata::hybrid_automaton_ptr>& list_automatons,
		std::string component_to_create = "");

}
}

#endif

