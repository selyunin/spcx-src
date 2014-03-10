/*
 * state_grammar.cpp
 *
 *  Created on: Mar 28, 2010
 *      Author: gvincent
 */

#include "state_grammar.h"

#include "parse_exception.h"

namespace predicate_parser {

state_ast& state_ast::operator&&(state_ast const& pred)
    {
	my_pred = state_binary_op(AND, my_pred, pred);
        return *this;
    }

state_ast& state_ast::operator||(state_ast const& pred)
    {
    my_pred = state_binary_op(OR, my_pred, pred);
        return *this;
    }

state_ast& state_ast::operator!()
    {
	my_pred = state_binary_op(NOT, my_pred);
        return *this;
    }


state_loc& state_loc::operator==(std::string const& exp)
   {
		loc = exp;
		eq = true;
	    return *this;
   }

state_loc& state_loc::operator!=(std::string const& exp)
   {
		loc = exp;
		eq = false;
        return *this;
   }

state_loc::state_loc(){
	aut = "";
	loc = "";
}

state_definition::state_definition()
    {
                using qi::_val;
                using qi::_1;
                using qi::lit;
                using qi::char_;
                using ascii::string;

                loc_equa = lit("loc(") > -(token[_val = _1]) > lit(")") >
						( (lit("==") > token [_val == _1])
                		| (lit("!=") > token [_val != _1]));

                factor_loc = loc_equa [_val = _1]
                		     | ( '(' >> or_loc [_val = _1] >> ')')
                		     | factor_equa[_val = _1]
                		     | ('!' >> factor_loc[_val = !_1]);


           		and_loc = factor_loc  [_val = _1]
           		                         >> *(  ( lit("&&") >> factor_loc [_val && _1])
           		                        		 | (lit("&") >> factor_loc [_val && _1]));

           		or_loc = and_loc [_val = _1]
       		                         >> *(  ( lit("||") >> and_loc [_val || _1])
       		                        		 | (lit("|") >> and_loc [_val || _1]));

                // for error handling
                using qi::fail;
                using qi::on_error;
                using namespace qi::labels;
				boost::phoenix::function<parser::parse_exception_handler_impl> parse_exception_handler;

                on_error<fail>( loc_equa,parse_exception_handler("location constraint",_1,_2,_3,_4));
                on_error<fail>( factor_loc,parse_exception_handler("state predicate",_1,_2,_3,_4));
                on_error<fail>( and_loc,parse_exception_handler("state predicate",_1,_2,_3,_4));
                on_error<fail>( or_loc,parse_exception_handler("state predicate",_1,_2,_3,_4));
     }

}
