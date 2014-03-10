/*
 * expression_grammar.h
 *
 *  Created on: Sep 4, 2009
 *      Author: frehse
 */

#ifndef EXPRESSION_GRAMMAR_H_
#define EXPRESSION_GRAMMAR_H_

#include <string>
#include <ctype.h>
#include <iostream>
#include <boost/algorithm/string.hpp>


#include <boost/config/warning_disable.hpp>
#include <boost/variant.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_function.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/fusion/include/io.hpp>
#include "utility/operator_enums.h"
#include "io/common_input/scalar_node_creation.h"

namespace predicate_parser {

namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;

typedef std::string::const_iterator grammar_type;

    //  Our AST
    struct binary_op;
    struct nil {};

    struct expression_ast
    {
        typedef
            boost::variant<
                nil // can't happen!
              , std::string
              , boost::recursive_wrapper<expression_ast>
              , boost::recursive_wrapper<binary_op>
            >
        type;

        expression_ast()
          : expr(nil()) {}

        template <typename expr_type>
        expression_ast(expr_type const& expr)
          : expr(expr) {}

        expression_ast& operator+=(expression_ast const& exp);
        expression_ast& operator-=(expression_ast const& exp);
        expression_ast& operator*=(expression_ast const& exp);
        expression_ast& operator/=(expression_ast const& exp);
        expression_ast& operator&=(expression_ast const& exp);


        type expr;
    };

    struct binary_op
    {
        binary_op(arithmetic_operator ope
          , expression_ast const& c1
          , expression_ast const& c2)
        : op(ope), left(c1), right(c2) {}

        binary_op(arithmetic_operator ope
                 , expression_ast const& c1)
               : op(ope), left(c1), right() {}

        arithmetic_operator op;
        expression_ast left;
        expression_ast right;
    };


    // We should be using expression_ast::operator-. There's a bug
    // in phoenix type deduction mechanism that prevents us from
    // doing so. Phoenix will be switching to BOOST_TYPEOF. In the
    // meantime, we will use a phoenix::function below:
    struct negate_expr
    {
        template <typename T>
        struct result { typedef T type; };

        expression_ast operator()(expression_ast const& expr) const
        {
            return expression_ast(binary_op(NEG, expr));
        }
    };

    struct sqroot_expr
          {
              template <typename T>
              struct result { typedef T type; };

             expression_ast operator()(expression_ast const& expr) const
              {
            	 return expression_ast(binary_op(SQRT, expr));
              }
          };

        struct sin_expr
        {
             template <typename T>
             struct result { typedef T type; };

             expression_ast operator()(expression_ast const& expr) const
             {
                return expression_ast(binary_op(SIN, expr));
              }
        };

        struct log_expr
        {
              template <typename T>
              struct result { typedef T type; };

              expression_ast operator()(expression_ast const& expr) const
              {
                 return expression_ast(binary_op(LOG, expr));
               }
        };

        struct exp_expr
        {
        	template <typename T>
        	struct result { typedef T type; };

           expression_ast operator()(expression_ast const& expr) const
           {
              return expression_ast(binary_op(EXP, expr));
           }
        };

        struct cos_expr
        {
         	template <typename T>
         	struct result { typedef T type; };

            expression_ast operator()(expression_ast const& expr) const
            {
            	return expression_ast(binary_op(COS, expr));
            }
        };

        struct tan_expr
        {
          	template <typename T>
           	struct result { typedef T type; };

            expression_ast operator()(expression_ast const& expr) const
            {
            	return expression_ast(binary_op(TAN, expr));
            }
        };


    //  expression grammar

    struct expr_definition
    {
    	expr_definition();

        qi::rule<grammar_type, expression_ast(), ascii::space_type>
        expression, term, factor, pow,non_lin_factor;
        qi::rule<grammar_type, std::string(), ascii::space_type> token, mat_token, var_token,constant, sin_token,
        log_token, sqrt_token, exp_token, cos_token, tan_token;
    };

    struct expr_grammar : qi::grammar<grammar_type, expression_ast(), ascii::space_type>, expr_definition
        {
        	expr_grammar() : expr_grammar::base_type(expression){};

        };

}

#endif /* EXPRESSION_GRAMMAR_H_ */
