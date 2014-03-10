/*
 * state_grammar.h
 *
 *  Created on: Sep 27, 2009
 *      Author: gvincent
 */

#ifndef STATE_GRAMMAR_H_
#define STATE_GRAMMAR_H_

#include "predicate_grammar.h"

namespace predicate_parser {

namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;

    //  Our AST
    struct state_binary_op;
    struct state_loc;
    struct state_ast
    {
        typedef
            boost::variant<
                nil // can't happen!
			  , boost::recursive_wrapper<state_ast>
			  , boost::recursive_wrapper<predicate_ast>
              , boost::recursive_wrapper<state_binary_op>
			  , boost::recursive_wrapper<state_loc>
            >
        type;

        state_ast()
          : my_pred(nil()) {}

        template <typename expr_type>
        state_ast(expr_type const& pred)
          : my_pred(pred) {}

        //boolean operators
        state_ast& operator&&(state_ast const& pred);
        state_ast& operator||(state_ast const& pred);
        state_ast& operator!();

        type my_pred;
    };

    struct state_binary_op
    {
    	state_binary_op(
    			boolean_operator ope
          , state_ast const& c1
          , state_ast const& c2)
        : op(ope), left(c1), right(c2) {}

    	state_binary_op(boolean_operator ope, state_ast const& c1) : op(ope), left(c1), right() {}

        boolean_operator op;
        state_ast left;
        state_ast right;
    };

    struct state_loc
        {
    	state_loc(bool ope
              , std::string const& c1
              , std::string const& c2)
            : eq(ope), aut(c1), loc(c2) {}

    	state_loc(std::string const& c1) : aut(c1){}
    	state_loc();

            bool eq;
            std::string aut;
            std::string loc;

            //comparison operators
            state_loc& operator==(std::string const& exp);
            state_loc& operator!=(std::string const& exp);
        };


    //  state grammar
    struct state_definition : predicate_definition
    {
    	state_definition();

    	qi::rule<grammar_type, state_ast(), ascii::space_type > and_loc, or_loc, factor_loc;
    	qi::rule<grammar_type, state_loc(), ascii::space_type > loc_equa;
    };

    struct state_grammar : qi::grammar<grammar_type, state_ast(), ascii::space_type>, state_definition
        {
    	state_grammar() : state_grammar::base_type(or_loc){};
        };


}

#endif /* STATE_GRAMMAR_H_ */
