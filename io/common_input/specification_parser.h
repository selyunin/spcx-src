/*
 * specification_parser.h
 *
 *  Created on: Sep 7, 2009
 *      Author: gvincent
 */

#ifndef SPECIFICATION_PARSER_H_
#define SPECIFICATION_PARSER_H_

#include "parse_type_chooser.h" // this is included for convenience (avoids include in user headers)
#include "io/common_input/symbol_table.h"

/** Forward declarations */
namespace tree {
class node;
typedef boost::shared_ptr<node> node_ptr;
}

namespace predicate_parser {

void parse_specification(std::string const& file,
		parser::symbol_table symbol_table = parser::symbol_table());

}

#endif /* SPECIFICATION_PARSER_H_ */
