/*
 * predicate_grammar.h
 *
 *  Created on: Sep 4, 2009
 *      Author: frehse
 */

#ifndef PREDICATE_GRAMMAR_H_
#define PREDICATE_GRAMMAR_H_

#include "expression_grammar.h"

namespace predicate_parser {

namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;

    //  Our AST
	struct predicate_binary_op;
	struct comparison_op;
    struct assignment_op;
    struct assignment_vector;

    struct predicate_ast
    {
        typedef
            boost::variant<
                nil // can't happen!
              , std::string
              , boost::recursive_wrapper<predicate_ast>
			  , boost::recursive_wrapper<expression_ast>
              , boost::recursive_wrapper<predicate_binary_op>
			  , boost::recursive_wrapper<comparison_op>
			  , boost::recursive_wrapper<assignment_op>
			  , boost::recursive_wrapper<assignment_vector>
            >
        type;

        predicate_ast()
          : my_pred(nil()) {}

        template <typename expr_type>
        predicate_ast(expr_type const& pred)
          : my_pred(pred) {}

        //boolean operators
        predicate_ast& operator&&(predicate_ast const& pred);
        predicate_ast& operator||(predicate_ast const& pred);
        predicate_ast& operator!();

        //comparison operators
        predicate_ast& operator==(expression_ast const& exp);
        predicate_ast& operator>(expression_ast const& exp);
        predicate_ast& operator<=(expression_ast const& exp);
        predicate_ast& operator>=(expression_ast const& exp);
        predicate_ast& operator<(expression_ast const& exp);

        predicate_ast& operator*=(expression_ast const& exp); //assignment operator

        type my_pred;
    };

    struct predicate_binary_op
    {
    	predicate_binary_op(boolean_operator ope
          , predicate_ast const& c1
          , predicate_ast const& c2)
        : op(ope), left(c1), right(c2) {}

    	predicate_binary_op(boolean_operator ope
    	          , predicate_ast const& c1)
    	        : op(ope), left(c1), right() {}

        boolean_operator op;
        predicate_ast left;
        predicate_ast right;
    };

    struct comparison_op
        {
    	comparison_op(comparison_operator ope
              , predicate_ast const& c1
              , expression_ast const& c2)
            : op(ope), left(c1), right(c2) {}

            comparison_operator op;
            predicate_ast left;
            predicate_ast right;
        };

    struct assignment_op
    {
    	assignment_op(predicate_ast const& c1
                  , expression_ast const& c2)
                : left(c1), right(c2) {}

                predicate_ast left;
                predicate_ast right;
     };

    struct assignment_vector
        {
			assignment_vector(std::vector<expression_ast> const& c1
                      , std::vector<expression_ast> const& c2)
                    : left(c1), right(c2) {}

			assignment_vector(std::vector<expression_ast> const& c1)
			                    : left(c1) {}
			assignment_vector(){}

                   std::vector<expression_ast> left;
                   std::vector<expression_ast> right;

                   assignment_vector& operator*=(std::vector<expression_ast> const& exp); //assignment vector operator (define right member)
         };

    //  predicate grammar
    struct predicate_definition : expr_definition
    {
    	predicate_definition();

    	qi::rule<grammar_type, predicate_ast(), ascii::space_type > and_equa, or_equa, factor_equa, equa, assignment_expression, tern_equa;
    	qi::rule<grammar_type, std::string(), ascii::space_type > bool_const;
    	qi::rule<grammar_type, assignment_vector(), ascii::space_type > assignment_vec;
    	qi::rule<grammar_type, std::vector<expression_ast>(), ascii::space_type > expression_vector;
    };

    struct predicate_grammar : qi::grammar<grammar_type, predicate_ast(), ascii::space_type>, predicate_definition
        {
        	predicate_grammar() : predicate_grammar::base_type(or_equa){};
        };


}

#endif /* PREDICATE_GRAMMAR_H_ */
