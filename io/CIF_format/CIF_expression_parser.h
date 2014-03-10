/*
 * expression_parser.h
 *
 *  Created on: Jul 1, 2009
 *      Author: gvincent
 */

#ifndef EXPRESSION_PARSER_H_

#include "boost/shared_ptr.hpp"
#include "io/common_input/symbol_table.h"
#include "io/common_input/parse_type_chooser.h" // this is included for convenience (avoids include in user headers)
#include "io/common_input/parse_policy.h"

/** Forward declarations */
namespace tree {
class node;
typedef boost::shared_ptr<node> node_ptr;
}

namespace CIF_predicate_parser {

/*!
 *  Parse the string expr and return a arithmetic tree::node::ptr, using boost::spirit).
 */
tree::node_ptr parse_expression(std::string expr,
		parser::symbol_table symbol_table = parser::symbol_table(), const parser::parse_policy& ppol = parser::parse_policy());

struct CIF_expression_ast;
tree::node_ptr make_expr_tree(CIF_expression_ast const& a,
		parser::symbol_table symbol_table = parser::symbol_table(),
		const parser::parse_policy& ppol = parser::parse_policy());

}
#endif /* EXPRESSION_PARSER_H_ */
