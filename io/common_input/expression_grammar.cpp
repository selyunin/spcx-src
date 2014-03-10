/*
 * expression_grammar.cpp
 *
 *  Created on: Mar 28, 2010
 *      Author: gvincent
 */

#include "expression_grammar.h"

#include <boost/spirit/include/phoenix_object.hpp>
#include "parse_exception.h"

namespace predicate_parser {

 expression_ast& expression_ast::operator+=(expression_ast const& exp)
    {
        expr = binary_op(ADD, expr, exp);
        return *this;
    }

    expression_ast& expression_ast::operator-=(expression_ast const& exp)
    {
        expr = binary_op(SUB, expr, exp);
        return *this;
    }

    expression_ast& expression_ast::operator*=(expression_ast const& exp)
    {
        expr = binary_op(MUL, expr, exp);
        return *this;
    }

    expression_ast& expression_ast::operator/=(expression_ast const& exp)
    {
        expr = binary_op(DIV, expr, exp);
        return *this;
    }

    expression_ast& expression_ast::operator&=(expression_ast const& exp)
    {
         expr = binary_op(POW, expr, exp);
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
	                using ascii::char_;
	                using ascii::digit;
	                using ascii::alpha;
	                using ascii::alnum;

                    expression =
                        term                            [_val = _1]
                        >> *(   ('+' >> term            [_val += _1])
                            |   ('-' >> term            [_val -= _1])
                            )
                        ;

                    term =
                        pow                          [_val = _1]
                        >> *(   ('*' >> pow          [_val *= _1])
                            |   ('/' >> pow          [_val /= _1])
                            )
                        ;

	                pow =
	                    factor                          [_val = _1]
	                     >> *('^' >> factor             [_val &= _1])
                        |   ('-' >> pow              	[_val = neg(_1)])
	                    ;

                    factor =
                    	constant                    	[_val = _1]
                       	|   ( '(' >> expression [_val = _1] >> ')')
                        |	non_lin_factor				[_val = _1]
                        | token                         [_val = _1]
                        ;

                    //PI_constant = PI_token >> expression[_val = pi(_1)];


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
                            >> (char_('(')[_val += _1] >> +ascii::digit[_val += _1]
                            >> -(char_(',')[_val += _1] >> +ascii::digit[_val += _1]) >> char_(')')[_val += _1] )
                            >> *(char_("'")[_val += _1])];

				constant %=  lexeme[qi::raw[-char_('-') >> (
											(&char_('.') >> char_('.') >> *digit) |
											(&digit >> +digit >> -char_('.') >> *digit)
										   )
							>> -(
									(char_('e')| char_('E')) >
									(+digit | (char_('+') > +digit) | (char_('-') > +digit))
								)]];
                // for error handling
                using qi::fail;
                using qi::on_error;
                using namespace qi::labels;
				boost::phoenix::function<parser::parse_exception_handler_impl> parse_exception_handler;

				on_error<fail>( constant,parse_exception_handler("number",_1,_2,_3,_4));
                on_error<fail>( token,parse_exception_handler("name",_1,_2,_3,_4));
                on_error<fail>( expression,parse_exception_handler("expression",_1,_2,_3,_4));
     }
}
