/*
 * specification_grammar.cpp
 *
 *  Created on: Mar 28, 2010
 *      Author: gvincent
 */

#include "specification_grammar.h"
#include <boost/spirit/include/phoenix_stl.hpp>

namespace predicate_parser {

specification_definition::specification_definition()
    {
                using qi::_val;
                using qi::_1;
                using boost::phoenix::push_back;

                symbolic_state %= '{' >> or_loc >> '}';

                arg = token[_val = _1] | symbolic_state[_val = _1];

                arguments = '(' >> (arg[push_back(_val, _1)] % ',') >> ')';

                com %=  token >> arguments;

                command_list = (com[push_back(_val, _1)] % ';') >> ';';
     }

}
