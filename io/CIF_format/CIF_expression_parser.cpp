/*
 * expression_parser.cpp
 *
 *  Created on: Sep 4, 2009
 *      Author: frehse
 */

#include <ctype.h>
//#include <boost/spirit/include/phoenix_stl.hpp>
#include "io/common_input/symbol_table.h"
#include "CIF_expression_parser.h"
#include "io/common_input/parse_type_chooser.h" // this is included for convenience (avoids include in user headers)
#include "core/predicates/valuation_function_tree_nodes.h"
#include "io/common_input/scalar_node_creation.h"
#include "CIF_expression_grammar.h"
#include "utility/basic_exception.h"

using namespace parser;

namespace CIF_predicate_parser {

/** Translate matrix or vector index notation to pair of indices.
 *
 * The indices start at 1.
 */
bool string_to_index_pair(const std::string& var_matrix_index, unsigned int& i,
		unsigned int& j) {
	// It's a matrix element
	i=1;
	j=1;
	std::string::const_iterator iter = var_matrix_index.begin();
	std::string::const_iterator end = var_matrix_index.end();
	bool r;
	if (var_matrix_index.find(",", 0) != std::string::npos) {
		r = qi::phrase_parse(iter, end, ('('
				>> qi::uint_[boost::phoenix::ref(i) = qi::_1] >> ','
				>> qi::uint_[boost::phoenix::ref(j) = qi::_1] >> ')'),
				ascii::space);
	} else {
		r = qi::phrase_parse(iter, end, ('('
				>> qi::uint_[boost::phoenix::ref(i) = qi::_1] >> ')'),
				ascii::space);
	}
	return !(iter != end || !r); // fail if we did not get a full match
}

//  Walk the tree
struct expression_tree {
	expression_tree(symbol_table sym_table, const parse_policy& ppol =
			parse_policy()) :
		my_parse_policy(ppol) {
		my_symbol_table = sym_table;
	}

	typedef tree::node::ptr result_type;

	tree::node_ptr operator()(nil) {
		return tree::node::null_node();
	}

	tree::node::ptr operator()(std::string s) {
		if (isdigit(s[0]) || s[0] == '-') {
			//it's a constant because a variable can't begin by a number or minus
			return predicate_parser::create_scalar_const_node(s);
		} else {
			/**Add context  */
			std::string context = my_symbol_table.get_context();

			//    		if(symbol_table.get_context() != "")
			//    		    s = symbol_table.get_context() + "." + s;

			/**
			 * If variable is primed, search the unprimed variable (and without index if it's a matrix) in symbol table.
			 * Next, add primes to the value of this variable (in variable node creation)
			 */
			std::string primes("");
			std::string var_matrix_index("");
			std::string unprimed_var(s);
			unsigned int prime_count = 0;
			variable::get_unprimed_name_and_add_prime_count(unprimed_var,
					prime_count);
			for (int i = 0; i < prime_count; ++i)
				primes = primes + "'";

			if (s.find("(") != std::string::npos) {
				std::string temp_s = unprimed_var.substr(0, unprimed_var.find(
						"("));
				var_matrix_index = unprimed_var.substr(unprimed_var.find("("),
						unprimed_var.size() - unprimed_var.find("("));
				unprimed_var = temp_s;
			}

			std::string
					unpr_var_in_context =
							valuation_functions::variable_node_creator::name_dot_context(
									unprimed_var, context);
			//std::cout << "treating token " << s << " unprimed: " << unprimed_var;
			if (my_symbol_table.is_symbol(unprimed_var)) {
				const symbol& sym = my_symbol_table.get_symbol(unprimed_var);

//std::cout << ", found symbol " << sym << std::endl;
				unsigned int i, j;
				if (var_matrix_index != "") {
					bool ok = string_to_index_pair(var_matrix_index, i, j);
					if (!ok) {
						throw std::runtime_error("Bad matrix coordinates in '"
								+ unprimed_var + var_matrix_index + "'.");
					}
					if (i > sym.dim1 || j > sym.dim2) {
						throw std::runtime_error(
								"Index exceeds dimensions in '" + unprimed_var
										+ var_matrix_index + "'.");
					}
				}

				// Resolve value if it's a constant
				if (sym.my_symbol_type == symbol::CONST_VALUE) {
					std::string value;
					/** If variable is a matrix element
					 */
					if (sym.dim1 != 1 || sym.dim2 != 1) {
						if (var_matrix_index != "") {
							value = boost::any_cast<std::string>(sym.my_value(i
									- 1, j - 1));
							return predicate_parser::create_scalar_const_node(value);
						} else {
							/** If symbol is a (entire) matrix of constant values */
							return predicate_parser::create_matrix_scalar_const_node(sym);
						}
					} else {
						// It's a constant scalar
						value
								= boost::any_cast<std::string>(sym.my_value(0,
										0));
						return predicate_parser::create_scalar_const_node(value);
					}
				} else {
					/** It's a variable */
					if (var_matrix_index != "") {
						// it's an element of a vector or matrix
						//std::cout << "instantiating scalar variable node " << sym.my_name + var_matrix_index + primes << std::endl;
						return tree::node::ptr(
								valuation_functions::variable_node_creator::create(
										sym.my_name + "(" + to_string((i-1)*sym.dim2+j)
												+ ")" + primes, "", 0,
										my_parse_policy));
					} else {
						// it's a vector or matrix
						//std::cout << "instantiating matrix variable node " << sym.my_name + primes << std::endl;
						return tree::node::ptr(
								valuation_functions::variable_node_creator::create(
										sym.my_name + primes, "", sym.dim1,
										sym.dim2, my_parse_policy));
					}
				}
			} else {
				/** It's not a known symbol. */
				if (!my_parse_policy.add_unknown_vars) {
					throw basic_exception("Can't find variable "
							+ unpr_var_in_context + " in component " + context
							+ ".");
					return tree::node::ptr();
				} else {
					return tree::node::ptr(
							valuation_functions::variable_node_creator::create(
									s, context, 0, my_parse_policy));
				}
			}

		}
	}

	tree::node_ptr operator()(CIF_expression_ast const& ast) {
		return boost::apply_visitor(*this, ast.expr);
	}

	tree::node_ptr operator()(binary_op const& expr) {
		return tree::node_ptr(new valuation_functions::arithmetic_node(expr.op,
				boost::apply_visitor(*this, expr.left.expr),
				boost::apply_visitor(*this, expr.right.expr)));
	}

	symbol_table my_symbol_table;
	const parse_policy& my_parse_policy;
};

tree::node_ptr make_expr_tree(CIF_expression_ast const& a,
		symbol_table symbol_table, const parser::parse_policy& ppol) {
	expression_tree make_tree(symbol_table, ppol);
	return make_tree(a);
}

tree::node::ptr parse_expression(std::string expr, symbol_table symbol_table, const parser::parse_policy& ppol) {
	if (expr == "") {
		tree::node::ptr p = tree::node::ptr(tree::node::null_node());
		return p;
	} else {
		expr_grammar g;
		CIF_expression_ast my_ast;
		std::string::const_iterator iter = expr.begin();
		std::string::const_iterator end = expr.end();

		bool r = boost::spirit::qi::phrase_parse(iter, end, g,
				boost::spirit::ascii::space, my_ast);

		if (r && iter == end) {
			return make_expr_tree(my_ast, symbol_table, ppol);
		} else {
			/**
			 * error handle : precise as to when the parser fails
			 * (part of the expression fails)
			 */
			throw basic_exception("Could not parse expression '" + expr
					+ "' \nfrom '" + std::string(iter, end) + "'.");
		}
	}
}
}
