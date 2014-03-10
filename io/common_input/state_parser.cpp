/*
 * state_parser.cpp
 *
 *  Created on: Sep 2, 2009
 *      Author: frehse
 */

#include "state_parser.h"

#include "utility/basic_exception.h"
#include "core/predicates/valuation_function_tree_nodes.h"
#include "core/hybrid_automata/location_eq_node.h"
#include "state_grammar.h"
#include "predicate_parser.h"

using namespace parser;

namespace predicate_parser {

//  Walk the tree
struct state_tree
{
	typedef tree::node_ptr result_type;
	symbol_table my_symbol_table;
	const parser::parse_policy& my_parse_policy;

	state_tree(symbol_table symb_table, const parser::parse_policy& ppol =
			parser::parse_policy()) :
		my_parse_policy(ppol) {
		my_symbol_table = symb_table;
	}

	tree::node_ptr operator()(nil) const {
		return tree::node::null_node();
	}

    tree::node_ptr operator()(predicate_ast const& ast) const
    {
		return make_predicate_tree(ast, my_symbol_table, my_parse_policy);
	}

    tree::node_ptr operator()(state_ast const& ast) const
        {
    	return boost::apply_visitor(*this, ast.my_pred);
        }

    tree::node_ptr operator()(state_binary_op const& pred) const
    {
    	return tree::node_ptr(new valuation_functions::boolean_node(pred.op,
    			boost::apply_visitor(*this, pred.left.my_pred), boost::apply_visitor(*this, pred.right.my_pred)));
    }

    tree::node_ptr operator()(state_loc const& pred) const {
		std::string autom(pred.aut);
		std::string context=my_symbol_table.get_context();
		if (autom == "") {
			autom = context;
			context = "";
		}

		return tree::node_ptr(hybrid_automata::location_node_creator::create(
				autom, context, pred.loc, pred.eq));
	}

};

tree::node_ptr make_state_tree(state_ast const& a, symbol_table sym_table,
		const parser::parse_policy& ppol){
	state_tree make_tree(sym_table,ppol);
	return make_tree(a);
}

tree::node::ptr parse_state(std::string expr, std::string const& context, symbol_table symbol_table,
		const parser::parse_policy& ppol) {
	if (expr == "") {
		tree::node::ptr p = tree::node::ptr(tree::node::null_node());
		return p;
	} else {
		symbol_table.set_context(context);
		state_grammar g;
		state_ast my_ast;
		std::string::const_iterator iter = expr.begin();
		std::string::const_iterator end = expr.end();

		bool r = boost::spirit::qi::phrase_parse(iter, end, g, boost::spirit::ascii::space, my_ast);

		if (r && iter == end) {
			return make_state_tree(my_ast, symbol_table,ppol);
		} else {
			/**
			* error handle : precise as to when the parser fails
			* (part of the predicate fails)
			*/
			throw basic_exception("Could not parse predicate:\n" + expr + "\n" + std::string(iter-expr.begin(),' ') + "^"+ std::string(iter, end));
		}
	}
}

}
