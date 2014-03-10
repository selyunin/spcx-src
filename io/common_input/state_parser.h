/*
 * state_parser.h
 *
 *  Created on: Jul 10, 2009
 *      Author: gvincent
 */

#ifndef STATE_PARSER_H_
#define STATE_PARSER_H_

#include "boost/shared_ptr.hpp"
#include "io/common_input/symbol_table.h"
#include "parse_type_chooser.h" // this is included for convenience (avoids include in user headers)
#include "parse_policy.h"

/** Forward declarations */
namespace tree {
class node;
typedef boost::shared_ptr<node> node_ptr;
}

namespace predicate_parser {

tree::node_ptr parse_state(std::string expr,
		std::string const& context = "",
		parser::symbol_table symbol_table = parser::symbol_table(),
		const parser::parse_policy& ppol = parser::parse_policy());

struct state_ast;
tree::node_ptr make_state_tree(state_ast const& a,
		parser::symbol_table symbol_table = parser::symbol_table(),
		const parser::parse_policy& ppol = parser::parse_policy());

}
#endif /* STATE_PARSER_H_ */
