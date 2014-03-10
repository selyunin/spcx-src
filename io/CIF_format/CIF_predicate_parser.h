#ifndef PREDICATE_PARSER_H_
#define PREDICATE_PARSER_H_

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

/** Parse the string expr and return a tree::ptr, using boost::spirit.
 */
tree::node_ptr parse_predicate(std::string expr,
		const parser::symbol_table& symbol_table = parser::symbol_table(),
		const parser::parse_policy& ppol = parser::parse_policy());

/** Parse the string expr and return a tree::ptr, looking up symbols within context.
 */
tree::node_ptr parse_predicate(std::string expr, std::string context,
		const parser::parse_policy& ppol = parser::parse_policy());

/** Parse string and return a tree with types bool and Rational. */
tree::node_ptr parse_predicate_bool_rational(std::string expr,
		const parser::symbol_table& symbol_table = parser::symbol_table(),
		const parser::parse_policy& ppol = parser::parse_policy());

/** Parse string and return a tree with type calc_string. */
tree::node_ptr parse_predicate_calc_string(std::string expr,
		const parser::symbol_table& symbol_table = parser::symbol_table(),
		const parser::parse_policy& ppol = parser::parse_policy());

/*!
 * \note The construction makes it illegal to use a == b to compare two boolean
 * expressions a and b for equality.
 * However, there is an equivalent expression that can be parsed:
 * (a && b) || !(a || b).
 */
struct predicate_ast;
tree::node_ptr make_predicate_tree(predicate_ast const& a,
		const parser::symbol_table& symbol_table = parser::symbol_table(),
		const parser::parse_policy& ppol = parser::parse_policy());

}
#endif
