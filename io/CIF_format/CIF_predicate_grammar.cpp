/*
 * predicate_grammar.cpp
 *
 *  Created on: Mar 28, 2010
 *      Author: gvincent
 */

#include "CIF_predicate_grammar.h"
#include <boost/spirit/include/phoenix_stl.hpp>

namespace CIF_predicate_parser {

predicate_ast& predicate_ast::operator&&(predicate_ast const& pred)
    {
	my_pred = predicate_binary_op(AND, my_pred, pred);
        return *this;
    }

predicate_ast& predicate_ast::operator||(predicate_ast const& pred)
    {
    my_pred = predicate_binary_op(OR, my_pred, pred);
        return *this;
    }

predicate_ast& predicate_ast::operator!()
    {
	my_pred = predicate_binary_op(NOT, my_pred);
        return *this;
    }

//comparison operators
predicate_ast& predicate_ast::operator==(CIF_expression_ast const& exp)
   {
	my_pred = comparison_op(EQ, my_pred, exp);
	    return *this;
   }

predicate_ast& predicate_ast::operator>(CIF_expression_ast const& exp)
   {
    my_pred = comparison_op(GT, my_pred, exp);
       return *this;
   }

predicate_ast& predicate_ast::operator<=(CIF_expression_ast const& exp)
   {
	my_pred = comparison_op(LE, my_pred, exp);
		    return *this;
   }

predicate_ast& predicate_ast::operator>=(CIF_expression_ast const& exp)
   {
	my_pred = comparison_op(GE, my_pred, exp);
		    return *this;
   }

predicate_ast& predicate_ast::operator<(CIF_expression_ast const& exp)
   {
	my_pred = comparison_op(LT, my_pred, exp);
		    return *this;
   }


predicate_ast& predicate_ast::operator*=(CIF_expression_ast const& exp)
   {
	my_pred = assignment_op(my_pred, exp);
		    return *this;
   }

assignment_vector& assignment_vector::operator*=(std::vector<CIF_expression_ast> const& exp)
   {
			right = exp;
		    return *this;
   }

predicate_definition::predicate_definition()
    {
                using qi::_val;
                using qi::_1;
                using qi::lit;
                using qi::char_;
                using ascii::string;
                using boost::phoenix::push_back;


                equa = expression[_val = _1] >>
						( (lit(">") >> expression [_val > _1])
                		| (lit(">=") >> expression [_val >= _1])
                		| (lit("=") >> expression [_val == _1])
                		| (lit("<") >> expression [_val < _1])
						| (lit("<=") >> expression [_val <= _1]));

                tern_equa = equa[_val = _1] >>
                						( (lit(">") >> expression [_val > _1])
                                		| (lit(">=") >> expression [_val >= _1])
                                		| (lit("<") >> expression [_val < _1])
                						| (lit("<=") >> expression [_val <= _1]));

                factor_equa =
                		tern_equa [_val = _1] | equa [_val = _1] | assignment_expression[_val = _1] | bool_const[_val = _1] | assignment_vec[_val = _1]
                                    |   '(' >> or_equa           [_val = _1] >> ')'
                                    |   ('!' >> factor_equa              [_val = !_1]);


           		and_equa = factor_equa  [_val = _1]
           		                         >> *(  ( lit("&&") >> factor_equa [_val && _1])
           		                        		 | (lit("&") >> factor_equa [_val && _1]));

           		or_equa = and_equa [_val = _1]
       		                         >> *(  ( lit("||") >> and_equa [_val || _1])
       		                        		 | (lit("|") >> and_equa [_val || _1]));

           		assignment_expression = expression[_val = _1] >> ( (lit(":=") >> expression [_val *= _1]));

           		bool_const = string("true")[_val = _1] | string("false")[_val = _1];

           		assignment_vec = expression_vector[_val = _1] >> (  ( lit("=") >> expression_vector [_val *= _1])
						| (lit(":=") >> expression_vector [_val *= _1]));

           		expression_vector = lit('(') >> (expression[push_back(_val, _1)] % ',') >> lit(')');
     }

}
