/*
 * parser_basics.h
 *
 *  Created on: Mar 1, 2011
 *      Author: frehse
 */

#ifndef PARSER_BASICS_H_
#define PARSER_BASICS_H_

#include <string>

namespace parser {

/** Returns true if the string can be interpreted as a number. */
bool is_number(const std::string& value);

}

#endif /* PARSER_BASICS_H_ */
