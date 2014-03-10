/*
 * expression_grammar.cpp
 *
 *  Created on: Mar 28, 2010
 *      Author: gvincent
 */

#include <boost/spirit/include/phoenix_object.hpp>
#include "CIF_expression_grammar.h"

namespace CIF_predicate_parser {

CIF_expression_ast& CIF_expression_ast::operator+=(CIF_expression_ast const& exp)
    {
        expr = binary_op(ADD, expr, exp);
        return *this;
    }

CIF_expression_ast& CIF_expression_ast::operator-=(CIF_expression_ast const& exp)
    {
        expr = binary_op(SUB, expr, exp);
        return *this;
    }

CIF_expression_ast& CIF_expression_ast::operator*=(CIF_expression_ast const& exp)
    {
        expr = binary_op(MUL, expr, exp);
        return *this;
    }

CIF_expression_ast& CIF_expression_ast::operator/=(CIF_expression_ast const& exp)
    {
        expr = binary_op(DIV, expr, exp);
        return *this;
    }


		boost::phoenix::function<negate_expr> neg;
		boost::phoenix::function<sqroot_expr> sqr;
		boost::phoenix::function<sin_expr> sin;
		boost::phoenix::function<log_expr> log;
		boost::phoenix::function<exp_expr> exp;
		boost::phoenix::function<cos_expr> cos;
		boost::phoenix::function<tan_expr> tan;


    expr_definition::expr_definition()
    {
                using qi::_val;
                using qi::_1;
                using qi::lexeme;
                using ascii::char_;
                using qi::lit;

                expression =
                    term                            [_val = _1]
                    >> *(   ('+' >> term            [_val += _1])
                        |   ('-' >> term            [_val -= _1])
                        )
                    ;

                term =
                    factor                          [_val = _1]
                    >> *(   ('*' >> factor          [_val *= _1])
                        |   ('/' >> factor          [_val /= _1])
                        )
                    ;

                factor =
                	constant                    	[_val = _1]
                	|   ( '(' >> expression [_val = _1] >> ')')
                    |   ('-' >> factor              [_val = neg(_1)])
                    |	non_lin_factor				[_val = _1]
                    | token                         [_val = _1]
                    ;


                sin_token = lexeme[qi::eps[_val = ""]
                             >> (char_('s')[_val += _1]
                             >> char_('i')[_val += _1] >> char_('n')[_val += _1])];

                cos_token = lexeme[qi::eps[_val = ""]
                            >> (char_('c')[_val += _1]
                            >> char_('o')[_val += _1] >> char_('s')[_val += _1])];

                tan_token = lexeme[qi::eps[_val = ""]
                            >> (char_('t')[_val += _1]
                            >> char_('a')[_val += _1] >> char_('n')[_val += _1])];

                log_token = lexeme[qi::eps[_val = ""]
                            >> (char_('l')[_val += _1]
                            >> char_('o')[_val += _1] >> char_('g')[_val += _1])];

                sqrt_token = lexeme[qi::eps[_val = ""]
                             >> (char_('s')[_val += _1]
                             >> char_('q')[_val += _1] >> char_('r')[_val += _1]) >> char_('t')[_val += _1]];

                exp_token = lexeme[qi::eps[_val = ""]
                            >> (char_('e')[_val += _1]
                            >> char_('x')[_val += _1] >> char_('p')[_val += _1])];

                non_lin_factor = (sin_token >>  '(' >> expression [_val = sin(_1)] >> ')')
                               | (sqrt_token >>  '(' >> expression [_val = sqr(_1)] >> ')')
                               | (log_token >>  '(' >> expression [_val = log(_1)] >> ')')
                               | (exp_token >>  '(' >> expression [_val = exp(_1)] >> ')')
                               | (cos_token >>  '(' >> expression [_val = cos(_1)] >> ')')
                               | (tan_token >>  '(' >> expression [_val = tan(_1)] >> ')');


                /** variable token can either be a matrix or normal token */
                var_token = token					[_val = _1]
                		 | mat_token				[_val = _1]
                		 ;

                token = lexeme[qi::eps[_val = ""]
                     >> (ascii::alpha[_val += _1] || char_('_')[_val += _1])
                     >> *(ascii::alnum[_val += _1] || char_('_')[_val += _1] || char_('.')[_val += _1])
                     >> *(char_("'")[_val += _1])];

                /** For the tokens of type matrix */
                mat_token = lexeme[qi::eps[_val = ""]
                        >> (ascii::alpha[_val += _1] || char_('_')[_val += _1])
                        >> *(ascii::alnum[_val += _1] || char_('_')[_val += _1] || char_('.')[_val += _1])
                        >> -(char_('(')[_val += _1] >> +ascii::digit[_val += _1]
                        >> -(char_(',')[_val += _1] >> +ascii::digit[_val += _1]) >> char_(')')[_val += _1] )
                        >> *(char_("'")[_val += _1])];

                constant = lexeme[qi::eps[_val = ""] >> -(char_('-')[_val += _1]) >> +ascii::digit[_val += _1] >> -(char_('.')[_val += _1] >> +ascii::digit[_val += _1])];
     }

}
