/*
 * predicate_parser.cpp
 *
 *  Created on: Sep 4, 2009
 *      Author: frehse
 */

#include "CIF_predicate_parser.h"
#include "io/common_input/bool_node_creation.h"
#include "CIF_expression_parser.h"
#include "io/common_input/parse_type_chooser.h" // this is included for convenience (avoids include in user headers)
#include "CIF_predicate_grammar.h"
#include "core/predicates/valuation_function_tree_nodes.h"
#include "core/predicates/node_print_visitor.h"


namespace CIF_predicate_parser {

tree::node_ptr parse_predicate_bool_rational(std::string expr, const parser::symbol_table& symbol_table, const parser::parse_policy& ppol) {
	parse_type_chooser::set_bool(global_types::STD_BOOL);
	parse_type_chooser::set_number(global_types::GMP_RATIONAL);
	return parse_predicate(expr, symbol_table, ppol);
}

tree::node_ptr parse_predicate_calc_string(std::string expr, const parser::symbol_table& symbol_table, const parser::parse_policy& ppol) {
	parse_type_chooser::set_bool(global_types::CALC_STR);
	parse_type_chooser::set_number(global_types::CALC_STR);
	return parse_predicate(expr, symbol_table, ppol);
}


//  Walk the tree
struct predicate_tree
{
	predicate_tree(parser::symbol_table symb_table, const parser::parse_policy& ppol) :  my_policy(ppol) {
		my_symbol_table=symb_table;
	}

	typedef tree::node_ptr result_type;
	parser::symbol_table my_symbol_table;
	parser::parse_policy my_policy;

    tree::node_ptr operator()(nil)  {
    	return tree::node::null_node();
    }

    tree::node_ptr operator()(predicate_ast const& ast)
    {
        return boost::apply_visitor(*this, ast.my_pred);
    }

    tree::node_ptr operator()(CIF_expression_ast const& ast) {

		return make_expr_tree(ast.expr, my_symbol_table, my_policy);
	}

    tree::node_ptr operator()(predicate_binary_op const& pred)
    {
    	return tree::node_ptr(new valuation_functions::boolean_node(pred.op,
    			boost::apply_visitor(*this, pred.left.my_pred), boost::apply_visitor(*this, pred.right.my_pred)));
    }

    tree::node_ptr operator()(comparison_op const& pred)
       {
    	tree::node_ptr left = boost::apply_visitor(*this, pred.left.my_pred);
    	if (valuation_functions::comparison_node* c= dynamic_cast<valuation_functions::comparison_node*>(left.get()))
    	{
    		//check that is a correct ternary node
    		if(((c->my_op == GT || c->my_op == GE) && (pred.op != GT && pred.op != GE))
    				|| ((c->my_op == LT || c->my_op == LE) && (pred.op != LT && pred.op != LE))) {
    			std::stringstream s;
    			s << left << " " << sign_string(pred.op) << " ... ";
    			throw basic_exception("Cannot handle the following ternary expression:\n"+s.str());
    		}
    		tree::node_ptr right = tree::node_ptr(new valuation_functions::comparison_node(pred.op, c->child2, boost::apply_visitor(*this, pred.right.my_pred)));
    		return tree::node_ptr(new valuation_functions::boolean_node(AND, left, right));
    	}
    	else
    		return tree::node_ptr(new valuation_functions::comparison_node(pred.op, left, boost::apply_visitor(*this, pred.right.my_pred)));
       }

    tree::node_ptr operator()(assignment_op const& pred)
    {
        	tree::node_ptr var = boost::apply_visitor(*this, pred.left.my_pred);
        	if (valuation_functions::variable_node* var_left = dynamic_cast<valuation_functions::variable_node*>(var.get())){
        		tree::node_ptr var_primed = tree::node_ptr(new valuation_functions::variable_node(variable::get_id_primedness_increased(var_left->my_id)));
        		return tree::node_ptr(new valuation_functions::comparison_node(EQ, var_primed, boost::apply_visitor(*this, pred.right.my_pred)));
        	}
        	else {
    			std::stringstream s;
    			s << var;
    			throw basic_exception("Left side of an assignment needs to be a variable. Cannot handle assignment to "+s.str()+".");
        	}
    }

    tree::node_ptr operator()(assignment_vector const& pred)
    {
    	if (pred.left.size() != pred.right.size()) {
    		throw basic_exception(
					"Cannot handle vector assignment with vectors of different length ("
							+ to_string(pred.left.size()) + " and "
							+ to_string(pred.right.size()) + ")");
		}

    	predicate_tree make_assignments(my_symbol_table,my_policy);
    	predicate_ast var(pred.left[0]);
    	tree::node_ptr new_node = make_assignments(assignment_op(var, pred.right[0]));
    	for (std::vector<CIF_expression_ast>::size_type i = 1; i < pred.left.size(); ++i){
    		predicate_ast var2(pred.left[i]);
    		tree::node_ptr new_node2 =  make_assignments(assignment_op(var2, pred.right[i]));
    		new_node = tree::node_ptr(new valuation_functions::boolean_node(AND, new_node, new_node2));
    	}
    	return new_node;
   }

    tree::node_ptr operator()(std::string s)  {
        	if(s == "true" || s == "false")
        		return predicate_parser::create_bool_const_node(s);
        	else {
        		throw basic_exception("Cannot parse "+s+". Note that boolean variables are not handled.");
        	}
        }
};

tree::node_ptr make_predicate_tree(predicate_ast const& a, const parser::symbol_table& symbol_table, const parser::parse_policy& ppol){
	predicate_tree make_tree(symbol_table,ppol);
	return make_tree(a);
}

tree::node::ptr parse_predicate(std::string expr, const parser::symbol_table& symbol_table, const parser::parse_policy& ppol) {
	if (expr == "") {
		tree::node::ptr p = tree::node::ptr(tree::node::null_node());
		return p;
	} else {
		predicate_grammar g;
		predicate_ast my_ast;
		std::string::const_iterator iter = expr.begin();
		std::string::const_iterator end = expr.end();

		bool r = boost::spirit::qi::phrase_parse(iter, end, g, boost::spirit::ascii::space, my_ast);

		if (r && iter == end) {
			return make_predicate_tree(my_ast, symbol_table, ppol);
		} else {
			/**
			* error handle : precise as to when the parser fails
			* (part of the predicate fails)
			*/
			throw basic_exception("Could not parse predicate '" + expr + "' \nfrom '" + std::string(iter, end) + "'.");
		}
	}
}

tree::node::ptr parse_predicate(std::string expr, std::string context, const parser::parse_policy& ppol) {
	parser::symbol_table symbol_table;
	symbol_table.set_context(context);
	return parse_predicate(expr,symbol_table, ppol);
}

}
