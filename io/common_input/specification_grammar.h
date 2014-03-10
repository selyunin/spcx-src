/*
 * specification_grammar.h
 *
 *  Created on: Sep 27, 2009
 *      Author: gvincent
 */

#ifndef SPECIFICATION_GRAMMAR_H_
#define SPECIFICATION_GRAMMAR_H_

#include "state_grammar.h"

namespace predicate_parser {

namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;

    struct argument
    {
        typedef
            boost::variant<
                nil // can't happen!
              , std::string
			  , boost::recursive_wrapper<state_ast>
            >
        type;

        argument()
          : my_arg(nil()) {}

        template <typename expr_type>
        argument(expr_type const& pred)
          : my_arg(pred) {}

        type my_arg;
    };

    struct command
    {
    	command(std::string const& name, std::vector<argument> const& arguments) : my_name(name), my_arguments(arguments){}
    	command(){}

        std::string my_name;
        std::vector<argument> my_arguments;
    };
}

BOOST_FUSION_ADAPT_STRUCT(
		predicate_parser::command,
    (std::string, my_name)
    (std::vector<predicate_parser::argument>, my_arguments)
)

namespace predicate_parser {


    //  specification grammar
    struct specification_definition : state_definition
    {
    	specification_definition();

    	qi::rule<grammar_type, state_ast(), ascii::space_type > symbolic_state;
    	qi::rule<grammar_type, argument(), ascii::space_type > arg;
    	qi::rule<grammar_type, std::vector<argument>(), ascii::space_type > arguments;
    	qi::rule<grammar_type, command(), ascii::space_type > com;
    	qi::rule<grammar_type, std::vector<command>(), ascii::space_type > command_list;
    };

    struct specification_grammar : qi::grammar<grammar_type, std::vector<command>(), ascii::space_type>, specification_definition
        {
    	specification_grammar() : specification_grammar::base_type(command_list){};
        };


}

#endif /* SPECIFICATION_GRAMMAR_H_ */
